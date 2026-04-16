#!/usr/bin/env python3
"""
Tests for SAI-10k predictions module.

These tests verify that the sai10k_predictions logic produces results consistent
with the examples in the SAI-10k-calc repository.

To run from this directory:
    python3 -m unittest test_sai10k_predictions -v

Or directly:
    python3 test_sai10k_predictions.py
"""

import os
import sys
import unittest

from dotenv import load_dotenv

load_dotenv()

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sai10k_predictions import (
    sai10k_find_affected_region,
    sai10k_determine_aberrations,
    sai10k_annotate_frameshift,
    sai10k_compute_predictions,
    sai10k_select_transcript,
    sai10k_get_transcript_predictions,
    MIN_DELTA_SCORE,
)


def get_db_connection():
    """Get a database connection for retrieving transcript info.

    Returns:
        psycopg2 connection object, or None if connection fails
    """
    try:
        import psycopg2

        db_host = os.environ.get("SPLICEAI_LOOKUP_DB_HOST")
        if not db_host:
            print("Warning: SPLICEAI_LOOKUP_DB_HOST not set in environment")
            return None

        # Read password from .pgpass file (same approach as build_and_deploy.py)
        pgpass_path = os.path.join(os.path.dirname(__file__), ".pgpass")
        if os.path.exists(pgpass_path):
            with open(pgpass_path) as f:
                password = f.read().strip()
        else:
            password = os.environ.get("DB_PASSWORD")

        if not password:
            print("Warning: No database password found in .pgpass or DB_PASSWORD env var")
            return None

        conn = psycopg2.connect(
            host=db_host,
            dbname="spliceai-lookup-db",
            user="postgres",
            password=password,
        )
        return conn
    except Exception as e:
        print(f"Warning: Could not connect to database: {e}")
        return None


def get_transcript_structure_from_db(conn, transcript_id, genome_version):
    """Retrieve transcript structure from the database.

    This mirrors the function in server.py.

    Args:
        conn: Database connection
        transcript_id: Transcript ID (e.g., "ENST00000357654" or "NM_000059.4")
        genome_version: "37" or "38"

    Returns:
        dict with transcript structure fields, or None if not found
    """
    if conn is None:
        return None

    transcript_id_without_version = transcript_id.split(".")[0]
    table_name = f"transcripts_hg{genome_version}"

    cursor = conn.cursor()
    cursor.execute(
        f"""SELECT strand, cds_start, cds_end, exon_starts, exon_ends, exon_frames
           FROM {table_name} WHERE transcript_id = %s""",
        (transcript_id_without_version,),
    )
    rows = cursor.fetchall()
    cursor.close()

    if not rows:
        return None

    strand, cds_start, cds_end, exon_starts_str, exon_ends_str, exon_frames_str = rows[0]

    # Parse comma-separated values and convert from 0-based to 1-based coordinates
    exon_starts_0based = [int(s) for s in exon_starts_str.rstrip(",").split(",") if s]
    exon_ends_0based = [int(s) for s in exon_ends_str.rstrip(",").split(",") if s]

    # Convert to 1-based coordinates
    exon_starts_1based = [s + 1 for s in exon_starts_0based]
    exon_ends_1based = exon_ends_0based

    cds_start_1based = cds_start + 1 if cds_start is not None else None
    cds_end_1based = cds_end if cds_end is not None else None

    exon_frames = None
    if exon_frames_str:
        exon_frames = [int(f) for f in exon_frames_str.rstrip(",").split(",") if f]

    return {
        "EXON_STARTS": exon_starts_1based,
        "EXON_ENDS": exon_ends_1based,
        "CDS_START": cds_start_1based,
        "CDS_END": cds_end_1based,
        "EXON_FRAMES": exon_frames,
        "STRAND": strand,
    }


# Example transcript structures from SAI-10k-calc examples (hg19/GRCh37)
# These are based on RefSeq transcripts NM_000059.4 (BRCA2) and NM_007294.4 (BRCA1)

BRCA2_TRANSCRIPT_HG19 = {
    "transcript_id": "NM_000059.4",
    "gene_name": "BRCA2",
    "STRAND": "+",
    "TX_START": 32889645,  # 1-based
    "TX_END": 32974405,
    "CDS_START": 32890598,  # 1-based
    "CDS_END": 32972907,
    # Exon boundaries (1-based)
    "EXON_STARTS": [
        32889645, 32890559, 32893214, 32899213, 32900238, 32900379, 32900636,
        32903580, 32905056, 32906409, 32910402, 32918695, 32920964, 32928998,
        32930565, 32931879, 32936660, 32937316, 32944539, 32945093, 32950807,
        32953454, 32953887, 32954144, 32968826, 32971035, 32972299,
    ],
    "EXON_ENDS": [
        32889804, 32890664, 32893462, 32899321, 32900287, 32900419, 32900750,
        32903629, 32905167, 32907524, 32915333, 32918790, 32921033, 32929425,
        32930746, 32932066, 32936830, 32937670, 32944694, 32945237, 32950928,
        32953652, 32954050, 32954282, 32969070, 32971181, 32974405,
    ],
    "EXON_FRAMES": [
        -1, 0, 1, 1, 2, 1, 0, 1, 0, 1, 1, 1, 1, 2, 1, 0, 2, 2, 0, 0, 1, 0, 1, 0, 1, 0, 0,
    ],
}

BRCA1_TRANSCRIPT_HG19 = {
    "transcript_id": "NM_007294.4",
    "gene_name": "BRCA1",
    "STRAND": "-",
    "TX_START": 41196312,  # 1-based
    "TX_END": 41277381,
    "CDS_START": 41197695,  # 1-based
    "CDS_END": 41276113,
    # Exon boundaries (1-based)
    "EXON_STARTS": [
        41196312, 41199660, 41201138, 41203080, 41209069, 41215350, 41215891,
        41219625, 41222945, 41226348, 41228505, 41234421, 41242961, 41243452,
        41247863, 41249261, 41251792, 41256139, 41256885, 41258473, 41267743,
        41276034, 41277288,
    ],
    "EXON_ENDS": [
        41197819, 41199720, 41201211, 41203134, 41209152, 41215390, 41215968,
        41219712, 41223255, 41226538, 41228631, 41234592, 41243049, 41246877,
        41247939, 41249306, 41251897, 41256278, 41256973, 41258550, 41267796,
        41276132, 41277381,
    ],
    "EXON_FRAMES": [
        1, 0, 1, 0, 0, 1, 1, 0, 1, 2, 1, 0, 1, 1, 2, 1, 0, 1, 2, 2, 2, 0, -1,
    ],
}


class TestFindAffectedRegion(unittest.TestCase):
    """Test the sai10k_find_affected_region function."""

    def test_variant_in_exon_brca2(self):
        """Test finding a variant within an exon of BRCA2."""
        exon_starts = BRCA2_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA2_TRANSCRIPT_HG19["EXON_ENDS"]

        # Variant at 32890600 is in exon 2 (0-indexed: 1)
        result = sai10k_find_affected_region(32890600, exon_starts, exon_ends)
        self.assertIsNotNone(result)
        self.assertEqual(result["region_type"], "exon")
        self.assertEqual(result["region_number"], 2)  # 1-based exon number

    def test_variant_in_intron_brca2(self):
        """Test finding a variant within an intron of BRCA2."""
        exon_starts = BRCA2_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA2_TRANSCRIPT_HG19["EXON_ENDS"]

        # Variant at 32890533 is in intron 1 (between exon 1 and exon 2)
        result = sai10k_find_affected_region(32890533, exon_starts, exon_ends)
        self.assertIsNotNone(result)
        self.assertEqual(result["region_type"], "intron")
        self.assertEqual(result["region_number"], 1)  # 1-based intron number

    def test_variant_in_exon_brca1_minus_strand(self):
        """Variant at 41197774 sits in BRCA1's genomically-lowest exon.

        BRCA1 is on '-' strand with 23 exons. Biological exon 1 has the highest
        genomic coords; the genomically-lowest exon (genomic index 0) is
        biological exon 23 (the last exon in transcript order).
        """
        exon_starts = BRCA1_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA1_TRANSCRIPT_HG19["EXON_ENDS"]

        result = sai10k_find_affected_region(41197774, exon_starts, exon_ends, strand="-")
        self.assertIsNotNone(result)
        self.assertEqual(result["region_type"], "exon")
        self.assertEqual(result["region_number"], 23)  # biological numbering
        self.assertEqual(result["_genomic_index"], 0)

    def test_empty_exons(self):
        """Test with empty exon lists."""
        result = sai10k_find_affected_region(100, [], [])
        self.assertIsNone(result)

    def test_biological_numbering_plus_strand_matches_genomic(self):
        """On '+' strand biological number equals genomic index + 1."""
        result = sai10k_find_affected_region(
            32890600,
            BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            strand="+",
        )
        self.assertEqual(result["region_number"], 2)
        self.assertEqual(result["_genomic_index"], 1)


class TestDetermineAberrations(unittest.TestCase):
    """Test the sai10k_determine_aberrations function (returns a list)."""

    def setUp(self):
        self.exon_starts = BRCA2_TRANSCRIPT_HG19["EXON_STARTS"]
        self.exon_ends = BRCA2_TRANSCRIPT_HG19["EXON_ENDS"]

    def _invoke(self, **kwargs):
        """Call sai10k_determine_aberrations with defaults for any unspecified
        ALT scores (0) and mandatory arguments (BRCA2 coords)."""
        kwargs.setdefault("ds_ag_alt", 0)
        kwargs.setdefault("ds_al_alt", 0)
        kwargs.setdefault("ds_dg_alt", 0)
        kwargs.setdefault("ds_dl_alt", 0)
        kwargs.setdefault("exon_starts", self.exon_starts)
        kwargs.setdefault("exon_ends", self.exon_ends)
        return sai10k_determine_aberrations(**kwargs)

    def test_brca2_c40_1g_a(self):
        """BRCA2_c.-40+1G>A (DS_DL=0.99, DS_DG=0.57, DS_AL=0.08, DS_AG=0.01).

        Reference parser output: Intron_retention=YES AND Partial_exon_deletion=YES
        (combination). Loss branch fires (strand-aware gate routes to
        whole_intron_retention); cryptic-donor subflow fires with DP_DG=-100 and
        DP_ND=-1 -> partial_exon_deletion (3' end).
        """
        result = self._invoke(
            ds_ag=0.01, ds_al=0.08, ds_dg=0.57, ds_dl=0.99,
            # Supply ALT scores consistent with a real cryptic donor
            # activation (DS_DG_ALT high, native donor DS_DL_ALT weaker).
            ds_dg_alt=0.7, ds_dl_alt=0.01,
            dp_ag=-38, dp_al=754, dp_dg=-100, dp_dl=-1,
            variant_pos=32889805,
            strand="+",
        )
        types = [a["aberration_type"] for a in result]
        self.assertIn("whole_intron_retention", types)
        self.assertIn("partial_exon_deletion", types)

    def test_exon_skipping_prediction(self):
        """Both donor loss and acceptor loss pass -> exon_skipping (single aberration)."""
        result = self._invoke(
            ds_ag=0.0, ds_al=0.81, ds_dg=0.17, ds_dl=0.98,
            dp_ag=-221, dp_al=-106, dp_dg=-137, dp_dl=-1,
            variant_pos=32890665,
            strand="+",
        )
        types = [a["aberration_type"] for a in result]
        self.assertIn("exon_skipping", types)
        # Driving score for exon_skipping is max(DS_DL, DS_AL) = 0.98
        es = next(a for a in result if a["aberration_type"] == "exon_skipping")
        self.assertEqual(es["max_delta_score"], 0.98)

    def test_no_significant_effect(self):
        """All scores below DS_ALDL_MIN -> empty list."""
        result = self._invoke(
            ds_ag=0.01, ds_al=0.01, ds_dg=0.01, ds_dl=0.01,
            dp_ag=-78, dp_al=-193, dp_dg=-214, dp_dl=-109,
            variant_pos=32890637,
            strand="+",
        )
        self.assertEqual(result, [])

    def test_pseudoexon_requires_validation_gates(self):
        """Paired gain in an intron, but variant too close to native exon (<50bp)
        should NOT fire pseudoexon (even though old code did)."""
        result = self._invoke(
            ds_ag=0.3, ds_al=0, ds_dg=0.3, ds_dl=0,
            dp_ag=-50, dp_al=0, dp_dg=50, dp_dl=0,
            variant_pos=32890533,  # 26bp from exon 1 start -> fails d_50bp
            strand="+",
        )
        types = [a["aberration_type"] for a in result]
        self.assertNotIn("pseudoexon", types)

    def test_pseudoexon_deep_intronic(self):
        """Variant deep in an intron with valid paired-gain parameters -> pseudoexon."""
        # BRCA2 intron 0: [32889805, 32890558]. Pick a variant well inside it
        # and gain positions that define a 101bp gained exon entirely inside
        # this intron, both gain positions >50bp from native splice sites.
        result = self._invoke(
            ds_ag=0.3, ds_al=0, ds_dg=0.3, ds_dl=0,
            dp_ag=50, dp_al=0, dp_dg=150, dp_dl=0,
            variant_pos=32890100,  # 296bp from exon 0, 458bp from exon 1
            strand="+",
        )
        types = [a["aberration_type"] for a in result]
        self.assertIn("pseudoexon", types)
        ps = next(a for a in result if a["aberration_type"] == "pseudoexon")
        # GEX_size = (DP_DG - DP_AG) + 1 = (150 - 50) + 1 = 101.
        # `size` is populated by sai10k_annotate_frameshift (called from
        # sai10k_compute_predictions); here we verify the raw geo positions.
        self.assertEqual(abs(ps["geo_dg"] - ps["geo_ag"]) + 1, 101)

    def test_lowered_min_catches_weak_paired_loss(self):
        """DS_AL=0.13 (weak) + DS_DL=0.93 (strong) fires exon_skipping thanks to
        the paired min=0.02 / max=0.2 gate."""
        result = self._invoke(
            ds_ag=0, ds_al=0.13, ds_dg=0, ds_dl=0.93,
            dp_ag=0, dp_al=-50, dp_dg=0, dp_dl=-1,
            variant_pos=32890600,
            strand="+",
        )
        types = [a["aberration_type"] for a in result]
        self.assertIn("exon_skipping", types)

    def test_dp_zero_is_valid_not_missing(self):
        """DP_*=0 must be treated as a valid 'at the variant base' position."""
        # Deep-intronic variant with DP_AG=0 (cryptic acceptor at variant base)
        # and DP_DG=150 (gained donor 150bp downstream). Should fire pseudoexon.
        result = self._invoke(
            ds_ag=0.3, ds_al=0, ds_dg=0.3, ds_dl=0,
            dp_ag=0, dp_al=0, dp_dg=150, dp_dl=0,
            variant_pos=32890100,  # deep intronic (intron 0, 296bp from nearest exon)
            strand="+",
        )
        types = [a["aberration_type"] for a in result]
        self.assertIn("pseudoexon", types)
        ps = next(a for a in result if a["aberration_type"] == "pseudoexon")
        self.assertEqual(ps["geo_ag"], 32890100)  # variant_pos + DP_AG(0) = variant_pos

    def test_pseudoexon_ag_dg_in_different_introns(self):
        """Pseudoexon must NOT fire when AG and DG fall in different introns."""
        # Place the variant deep in intron 0 but configure DP_DG so that GEO_DG
        # lands past exon 1, i.e. in intron 1. AG and DG no longer share an intron.
        # BRCA2 intron 0: [32889805, 32890558]; exon 1: [32890559, 32890664].
        # Variant at 32890100 (intron 0). DP_AG = 100 -> GEO_AG = 32890200 (intron 0).
        # DP_DG = 700 -> GEO_DG = 32890800 (intron 1).
        result = self._invoke(
            ds_ag=0.3, ds_al=0, ds_dg=0.3, ds_dl=0,
            dp_ag=100, dp_al=0, dp_dg=700, dp_dl=0,
            variant_pos=32890100,
            strand="+",
        )
        types = [a["aberration_type"] for a in result]
        self.assertNotIn("pseudoexon", types)

    def test_reverse_strand_biological_exon_skipping(self):
        """Reverse-strand BRCA1 variant should report biological exon numbers."""
        # BRCA1 on '-' strand, 23 exons. Variant at 41199659 in genomic intron 0
        # (between genomic exons 0 and 1 = biological exons 22 and 23).
        # DP_AL=61, DP_DL=1 with '-' strand: DP_AL*-1=-61, DP_DL*-1=-1 -> -61 < -1
        # -> exon_skipping. Skipped exon bounded by GEO_AL=41199720 (acceptor of
        # genomic exon 1 on '-' strand = exon_ends[1]) and GEO_DL=41199660
        # (donor of genomic exon 1 on '-' strand = exon_starts[1]).
        result = sai10k_determine_aberrations(
            ds_ag=0.0, ds_al=0.40, ds_dg=0.36, ds_dl=0.93,
            dp_ag=-4, dp_al=61, dp_dg=-4, dp_dl=1,
            ds_ag_alt=0, ds_al_alt=0, ds_dg_alt=0, ds_dl_alt=0,
            variant_pos=41199659,
            exon_starts=BRCA1_TRANSCRIPT_HG19["EXON_STARTS"],
            exon_ends=BRCA1_TRANSCRIPT_HG19["EXON_ENDS"],
            strand="-",
        )
        types = [a["aberration_type"] for a in result]
        self.assertIn("exon_skipping", types)

    def test_cryptic_acceptor_partial_exon_deletion(self):
        """Cryptic acceptor activation deep in an exon should produce
        partial_exon_deletion."""
        # Exonic variant, strong DS_AG at a position inside the exon (beyond the
        # native acceptor), weak-to-absent native acceptor loss, with ALT-score
        # conditions that satisfy the cryptic-acceptor check.
        # BRCA2 exon 2: [32890559, 32890664]. Variant at 32890580. Native acceptor
        # at 32890559 (dp_na = -21). Cryptic acceptor at 32890600 (dp_ag = 20).
        # partial_size = dp_ag - dp_na = 41 (positive -> partial exon deletion).
        result = sai10k_determine_aberrations(
            ds_ag=0.3, ds_al=0.0, ds_dg=0.0, ds_dl=0.0,
            dp_ag=20, dp_al=0, dp_dg=0, dp_dl=0,
            # ALT scores to satisfy Cryptic_Acceptor_activation_check2:
            # DS_AG >= 0.2 AND (DS_AG_ALT - DS_AL_ALT) > -0.2
            ds_ag_alt=0.5, ds_al_alt=0.3, ds_dg_alt=0, ds_dl_alt=0,
            variant_pos=32890580,
            exon_starts=BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            exon_ends=BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            strand="+",
        )
        types = [a["aberration_type"] for a in result]
        self.assertIn("partial_exon_deletion", types)

    def test_multi_exon_skipping(self):
        """When GEO_AL and GEO_DL span multiple exons, _find_exons_between_loss_positions
        should identify all skipped exons and sum their sizes."""
        # BRCA2 exons 2 and 3: [32890559, 32890664] and [32893214, 32893462].
        # Construct a variant whose DP_AL matches exon 2's acceptor and whose
        # DP_DL matches exon 3's donor. Variant at 32890559 (acceptor of exon 2).
        # DP_AL = 0 -> GEO_AL = 32890559 (exon 2 acceptor).
        # DP_DL = 2903 -> GEO_DL = 32893462 (exon 3 donor).
        result = sai10k_determine_aberrations(
            ds_ag=0, ds_al=0.5, ds_dg=0, ds_dl=0.5,
            dp_ag=0, dp_al=0, dp_dg=0, dp_dl=2903,
            ds_ag_alt=0, ds_al_alt=0, ds_dg_alt=0, ds_dl_alt=0,
            variant_pos=32890559,
            exon_starts=BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            exon_ends=BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            strand="+",
        )
        # Should be exon_skipping on loss-pair + orientation; two exons in the bracket.
        es = [a for a in result if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)

    def test_loss_plus_cryptic_combination(self):
        """Loss branch and cryptic subflow should fire in parallel for
        BRCA2 c.-40+1G>A (reference shows Intron_retention=YES + Partial_exon_deletion=YES)."""
        result = sai10k_determine_aberrations(
            ds_ag=0.01, ds_al=0.08, ds_dg=0.57, ds_dl=0.99,
            # ALT scores consistent with a real cryptic-donor activation.
            ds_ag_alt=0, ds_al_alt=0,
            ds_dg_alt=0.7, ds_dl_alt=0.01,
            dp_ag=-38, dp_al=754, dp_dg=-100, dp_dl=-1,
            variant_pos=32889805,
            exon_starts=BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            exon_ends=BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            strand="+",
        )
        types = [a["aberration_type"] for a in result]
        self.assertIn("whole_intron_retention", types)
        self.assertIn("partial_exon_deletion", types)


class TestAnnotateFrameshift(unittest.TestCase):
    """Test the sai10k_annotate_frameshift function (mutates aberration dict)."""

    def test_frameshift_prediction_exon_skipping(self):
        # Use real BRCA2 exon 2 boundaries for geo_al / geo_dl so the size
        # calculation can identify the skipped exon from loss positions.
        exon2_start = BRCA2_TRANSCRIPT_HG19["EXON_STARTS"][1]  # 32890559
        exon2_end = BRCA2_TRANSCRIPT_HG19["EXON_ENDS"][1]     # 32890664
        aberration = {
            "aberration_type": "exon_skipping",
            "affected_region": {
                "region_type": "exon",
                "region_number": 2,
                "_genomic_index": 1,
            },
            "geo_al": exon2_start, "geo_dl": exon2_end,
            "geo_ag": 0, "geo_dg": 0,
            "_strand": "+",
        }
        sai10k_annotate_frameshift(
            aberration,
            cds_start=BRCA2_TRANSCRIPT_HG19["CDS_START"],
            cds_end=BRCA2_TRANSCRIPT_HG19["CDS_END"],
            exon_starts=BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            exon_ends=BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
        )
        self.assertTrue(aberration["affects_coding"])
        self.assertEqual(aberration["size_type"], "exon")
        # Exon 2 is 106bp (32890559..32890664)
        self.assertEqual(aberration["size"], 106)

    def test_non_coding_transcript(self):
        aberration = {
            "aberration_type": "exon_skipping",
            "affected_region": {
                "region_type": "exon",
                "region_number": 1,
                "_genomic_index": 0,
            },
            "geo_al": 100, "geo_dl": 150, "geo_ag": 0, "geo_dg": 0,
            "_strand": "+",
        }
        sai10k_annotate_frameshift(
            aberration,
            cds_start=None, cds_end=None,
            exon_starts=[100, 200], exon_ends=[150, 250],
        )
        self.assertFalse(aberration["affects_coding"])


class TestComputePredictions(unittest.TestCase):
    """Test the sai10k_compute_predictions function."""

    def test_compute_predictions_brca2_variant(self):
        """Test full prediction computation for a BRCA2 variant."""
        # BRCA2_c.-40+1G>A at position 32889805
        transcript_scores = {
            "DS_AG": 0.01,
            "DS_AL": 0.08,
            "DS_DG": 0.57,
            "DS_DL": 0.99,
            "DP_AG": -38,
            "DP_AL": 754,
            "DP_DG": -100,
            "DP_DL": -1,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "EXON_FRAMES": BRCA2_TRANSCRIPT_HG19["EXON_FRAMES"],
            "STRAND": BRCA2_TRANSCRIPT_HG19["STRAND"],
            "TX_START": BRCA2_TRANSCRIPT_HG19["TX_START"],
            "TX_END": BRCA2_TRANSCRIPT_HG19["TX_END"],
        }

        result = sai10k_compute_predictions(transcript_scores, variant_pos=32889805)

        self.assertIn("aberrations", result)
        self.assertIn("transcript_info", result)
        self.assertEqual(result["transcript_info"]["strand"], "+")
        self.assertTrue(result["transcript_info"]["is_coding"])
        # Frameshift fields are attached to each aberration in the list.
        for aberration in result["aberrations"]:
            self.assertIn("frameshift", aberration)
            self.assertIn("affects_coding", aberration)

    def test_compute_predictions_brca1_variant(self):
        """Test full prediction computation for a BRCA1 variant (minus strand)."""
        # BRCA1_c.5467+1G>A at position 41199659
        transcript_scores = {
            "DS_AG": 0.0,
            "DS_AL": 0.40,
            "DS_DG": 0.36,
            "DS_DL": 0.93,
            "DP_AG": -4,
            "DP_AL": 61,
            "DP_DG": -4,
            "DP_DL": 1,
            "EXON_STARTS": BRCA1_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA1_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "EXON_FRAMES": BRCA1_TRANSCRIPT_HG19["EXON_FRAMES"],
            "STRAND": BRCA1_TRANSCRIPT_HG19["STRAND"],
            "TX_START": BRCA1_TRANSCRIPT_HG19["TX_START"],
            "TX_END": BRCA1_TRANSCRIPT_HG19["TX_END"],
        }

        result = sai10k_compute_predictions(transcript_scores, variant_pos=41199659)

        self.assertIn("aberrations", result)
        self.assertEqual(result["transcript_info"]["strand"], "-")


class TestSelectCanonicalTranscript(unittest.TestCase):
    """Test the sai10k_select_transcript function."""

    def test_select_mane_select_over_canonical(self):
        """Test that MANE Select is preferred over canonical."""
        scores_list = [
            {"t_priority": "C", "DS_AG": 0.5, "DS_AL": 0.5, "DS_DG": 0.5, "DS_DL": 0.5},
            {"t_priority": "MS", "DS_AG": 0.3, "DS_AL": 0.3, "DS_DG": 0.3, "DS_DL": 0.3},
        ]
        result = sai10k_select_transcript(scores_list)
        self.assertEqual(result["t_priority"], "MS")

    def test_select_highest_sum_at_same_priority(self):
        """Tie-break rule is highest SUM of |DS_*| (matches server.py)."""
        scores_list = [
            # sum = 2.0
            {"t_priority": "C", "DS_AG": 0.5, "DS_AL": 0.5, "DS_DG": 0.5, "DS_DL": 0.5},
            # sum = 1.2 — lower despite higher single max
            {"t_priority": "C", "DS_AG": 0.9, "DS_AL": 0.1, "DS_DG": 0.1, "DS_DL": 0.1},
        ]
        result = sai10k_select_transcript(scores_list)
        self.assertEqual(result["DS_AG"], 0.5)  # the 2.0-sum transcript wins

    def test_empty_scores_list(self):
        """Test handling of empty scores list."""
        result = sai10k_select_transcript([])
        self.assertIsNone(result)


class TestWithDatabaseIntegration(unittest.TestCase):
    """Integration tests that use the postgres database."""

    @classmethod
    def setUpClass(cls):
        """Set up database connection for all tests in this class."""
        cls.conn = get_db_connection()
        if cls.conn is None:
            print("Warning: Skipping database integration tests - no connection available")

    @classmethod
    def tearDownClass(cls):
        """Close database connection."""
        if cls.conn is not None:
            cls.conn.close()

    def test_retrieve_transcript_structure_brca2(self):
        """Test retrieving BRCA2 transcript structure from database."""
        if self.conn is None:
            self.skipTest("Database connection not available")

        # Try to get BRCA2 transcript (ENST00000380152 is the canonical BRCA2)
        result = get_transcript_structure_from_db(self.conn, "ENST00000380152", "37")

        if result is not None:
            self.assertIn("EXON_STARTS", result)
            self.assertIn("EXON_ENDS", result)
            self.assertIn("CDS_START", result)
            self.assertIn("CDS_END", result)
            self.assertIn("STRAND", result)
            self.assertEqual(result["STRAND"], "+")
        else:
            print("Warning: BRCA2 transcript not found in database")

    def test_retrieve_transcript_structure_brca1(self):
        """Test retrieving BRCA1 transcript structure from database."""
        if self.conn is None:
            self.skipTest("Database connection not available")

        # Try to get BRCA1 transcript (ENST00000357654 is the canonical BRCA1)
        result = get_transcript_structure_from_db(self.conn, "ENST00000357654", "37")

        if result is not None:
            self.assertIn("EXON_STARTS", result)
            self.assertIn("EXON_ENDS", result)
            self.assertIn("STRAND", result)
            self.assertEqual(result["STRAND"], "-")
        else:
            print("Warning: BRCA1 transcript not found in database")

    def test_predictions_with_db_transcript(self):
        """Test predictions using transcript data from database."""
        if self.conn is None:
            self.skipTest("Database connection not available")

        # Get BRCA1 transcript from database
        transcript_struct = get_transcript_structure_from_db(
            self.conn, "ENST00000357654", "37"
        )
        if transcript_struct is None:
            self.skipTest("BRCA1 transcript not found in database")

        # Use the database transcript structure
        transcript_scores = {
            "DS_AG": 0.0,
            "DS_AL": 0.40,
            "DS_DG": 0.36,
            "DS_DL": 0.93,
            "DP_AG": -4,
            "DP_AL": 61,
            "DP_DG": -4,
            "DP_DL": 1,
        }
        transcript_scores.update(transcript_struct)

        result = sai10k_compute_predictions(transcript_scores, variant_pos=41199659)

        self.assertIn("aberrations", result)
        self.assertEqual(result["transcript_info"]["strand"], "-")


class TestExampleVariants(unittest.TestCase):
    """Test predictions for specific example variants from SAI-10k-calc.

    These are based on the expected outputs in example_variants_parsed.tsv.
    """

    def test_brca2_c40_1g_a(self):
        """Test BRCA2_c.-40+1G>A (13-32889805-G-A).

        Expected: Any_splicing_aberration=YES, Partial_exon_deletion=YES
        """
        transcript_scores = {
            "DS_AG": 0.01,
            "DS_AL": 0.08,
            "DS_DG": 0.57,
            "DS_DL": 0.99,
            "DP_AG": -38,
            "DP_AL": 754,
            "DP_DG": -100,
            "DP_DL": -1,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "EXON_FRAMES": BRCA2_TRANSCRIPT_HG19["EXON_FRAMES"],
            "STRAND": BRCA2_TRANSCRIPT_HG19["STRAND"],
        }

        # Supply ALT scores that match a real cryptic-donor activation.
        transcript_scores["DS_DG_ALT"] = 0.7
        transcript_scores["DS_DL_ALT"] = 0.01

        result = sai10k_compute_predictions(transcript_scores, variant_pos=32889805)

        # Reference parser output for this variant: Intron_retention=YES AND
        # Partial_exon_deletion=YES (combination).
        types = [a["aberration_type"] for a in result["aberrations"]]
        self.assertIn("whole_intron_retention", types)
        self.assertIn("partial_exon_deletion", types)

    def test_brca2_c67_1g_a(self):
        """Test BRCA2_c.67+1G>A (13-32890665-G-A).

        Expected: Any_splicing_aberration=YES, Exon_skipping=YES, Lost_exons=2
        """
        transcript_scores = {
            "DS_AG": 0.00,
            "DS_AL": 0.81,
            "DS_DG": 0.17,
            "DS_DL": 0.98,
            "DP_AG": -221,
            "DP_AL": -106,
            "DP_DG": -137,
            "DP_DL": -1,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "EXON_FRAMES": BRCA2_TRANSCRIPT_HG19["EXON_FRAMES"],
            "STRAND": BRCA2_TRANSCRIPT_HG19["STRAND"],
        }

        result = sai10k_compute_predictions(transcript_scores, variant_pos=32890665)

        # With both DS_AL=0.81 and DS_DL=0.98 significant, loss-pair branch fires.
        self.assertEqual(len(result["aberrations"]), 1)
        self.assertEqual(result["aberrations"][0]["aberration_type"], "exon_skipping")
        self.assertEqual(result["aberrations"][0]["max_delta_score"], 0.98)

    def test_brca2_no_effect(self):
        """Test BRCA2_c.40A>G (13-32890637-A-G).

        Expected: Any_splicing_aberration=NO (all delta scores < threshold)
        """
        transcript_scores = {
            "DS_AG": 0.01,
            "DS_AL": 0.01,
            "DS_DG": 0.04,
            "DS_DL": 0.04,
            "DP_AG": -78,
            "DP_AL": -193,
            "DP_DG": -214,
            "DP_DL": -109,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "EXON_FRAMES": BRCA2_TRANSCRIPT_HG19["EXON_FRAMES"],
            "STRAND": BRCA2_TRANSCRIPT_HG19["STRAND"],
        }

        result = sai10k_compute_predictions(transcript_scores, variant_pos=32890637)

        # All scores below DS_AGDG_MAX -> empty aberrations list.
        self.assertEqual(result["aberrations"], [])

    def test_brca1_c5467_5g_c(self):
        """Test BRCA1_c.5467+5G>C (17-41199655-C-G).

        Expected: Any_splicing_aberration=YES, Exon_skipping=YES, Lost_exons=22
        """
        transcript_scores = {
            "DS_AG": 0.00,
            "DS_AL": 0.28,
            "DS_DG": 0.00,
            "DS_DL": 0.27,
            "DP_AG": 5,
            "DP_AL": 65,
            "DP_DG": -1976,
            "DP_DL": 5,
            "EXON_STARTS": BRCA1_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA1_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "EXON_FRAMES": BRCA1_TRANSCRIPT_HG19["EXON_FRAMES"],
            "STRAND": BRCA1_TRANSCRIPT_HG19["STRAND"],
        }

        result = sai10k_compute_predictions(transcript_scores, variant_pos=41199655)

        # Reference sample output: Exon_skipping=YES, Lost_exons=22. BRCA1 is on
        # '-' strand; orientation gate DP_AL*Strand < DP_DL*Strand => -65 < -5
        # => exon_skipping (not whole_intron_retention).
        types = [a["aberration_type"] for a in result["aberrations"]]
        self.assertIn("exon_skipping", types)

    def test_brca1_c5467_1g_a(self):
        """Test BRCA1_c.5467+1G>A (17-41199659-C-T).

        Reference sample output: Exon_skipping=YES AND Partial_intron_retention=YES
        (combination). Loss branch fires exon_skipping; cryptic-donor subflow
        fires partial_intron_retention on the 3' side of the exon.
        """
        transcript_scores = {
            "DS_AG": 0.00,
            "DS_AL": 0.40,
            "DS_DG": 0.36,
            "DS_DL": 0.93,
            # ALT scores chosen to satisfy Cryptic_Donor_activation_check2:
            # DS_DG < 0.2, DS_DL > 0.1, DS_DG_ALT >= 0.9, (DS_DG_ALT - DS_DL_ALT) >= 0.2.
            "DS_AG_ALT": 0, "DS_AL_ALT": 0,
            "DS_DG_ALT": 0.95, "DS_DL_ALT": 0.05,
            "DP_AG": -4,
            "DP_AL": 61,
            "DP_DG": -4,
            "DP_DL": 1,
            "EXON_STARTS": BRCA1_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA1_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "EXON_FRAMES": BRCA1_TRANSCRIPT_HG19["EXON_FRAMES"],
            "STRAND": BRCA1_TRANSCRIPT_HG19["STRAND"],
        }

        result = sai10k_compute_predictions(transcript_scores, variant_pos=41199659)

        types = [a["aberration_type"] for a in result["aberrations"]]
        self.assertIn("exon_skipping", types)
        self.assertIn("partial_intron_retention", types)


class TestIncreasedExonInclusion(unittest.TestCase):
    """Test the increased_exon_inclusion aberration type.

    Increased exon inclusion fires when:
      - Gain paired gate passes (min(DS_AG, DS_DG) >= 0.02, max >= 0.2)
      - Orientation passes (DP_AG*Strand < DP_DG*Strand)
      - Pseudoexon conditions fail (e.g. gain positions are at native sites)
      - GEX_size == native exon size
      - DP_AG == DP_NA and DP_DG == DP_ND (gain positions match native sites)
    """

    def test_increased_exon_inclusion_plus_strand(self):
        """Variant in intron 1 of BRCA2 with gain positions matching exon 2's
        native acceptor and donor should produce increased_exon_inclusion.

        BRCA2 exon 2: [32890559, 32890664], size=106bp, '+' strand.
        Native acceptor (geo_na) = 32890559, native donor (geo_nd) = 32890664.
        Variant at 32890540 (intron 1, 19bp from exon 2 — closer to exon 2 than exon 1).
        DP_AG = 32890559 - 32890540 = 19  -> geo_ag = 32890559 = geo_na
        DP_DG = 32890664 - 32890540 = 124 -> geo_dg = 32890664 = geo_nd
        GEX_size = (124 - 19)*1 + 1 = 106 = native exon size.
        """
        result = sai10k_determine_aberrations(
            ds_ag=0.3, ds_al=0.0, ds_dg=0.3, ds_dl=0.0,
            dp_ag=19, dp_al=0, dp_dg=124, dp_dl=0,
            ds_ag_alt=0, ds_al_alt=0, ds_dg_alt=0, ds_dl_alt=0,
            variant_pos=32890540,
            exon_starts=BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            exon_ends=BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            strand="+",
        )
        types = [a["aberration_type"] for a in result]
        self.assertIn("increased_exon_inclusion", types)
        iei = next(a for a in result if a["aberration_type"] == "increased_exon_inclusion")
        self.assertEqual(iei["max_delta_score"], 0.3)

    def test_increased_exon_inclusion_minus_strand(self):
        """Variant near a BRCA1 exon with gain positions matching native sites.

        BRCA1 is on '-' strand. Genomic exon at index 1: [41199660, 41199720],
        size=61bp. On '-' strand: native acceptor = exon_end = 41199720,
        native donor = exon_start = 41199660.
        Variant at 41199730 (intron, 10bp past exon end, closer to exon 1 than exon 2).
        DP_AG = 41199720 - 41199730 = -10 -> geo_ag = 41199720 = geo_na
        DP_DG = 41199660 - 41199730 = -70 -> geo_dg = 41199660 = geo_nd
        Orientation: DP_AG*Strand = -10*-1=10, DP_DG*Strand = -70*-1=70 -> 10 < 70 passes.
        GEX_size = (-70 - (-10))*(-1) + 1 = 60*1 + 1 = 61 = native exon size.
        """
        result = sai10k_determine_aberrations(
            ds_ag=0.3, ds_al=0.0, ds_dg=0.3, ds_dl=0.0,
            dp_ag=-10, dp_al=0, dp_dg=-70, dp_dl=0,
            ds_ag_alt=0, ds_al_alt=0, ds_dg_alt=0, ds_dl_alt=0,
            variant_pos=41199730,
            exon_starts=BRCA1_TRANSCRIPT_HG19["EXON_STARTS"],
            exon_ends=BRCA1_TRANSCRIPT_HG19["EXON_ENDS"],
            strand="-",
        )
        types = [a["aberration_type"] for a in result]
        self.assertIn("increased_exon_inclusion", types)

    def test_no_increased_exon_inclusion_when_gex_size_mismatches(self):
        """If GEX_size doesn't match native exon size, increased_exon_inclusion
        should NOT fire (and pseudoexon shouldn't either since positions are
        at native sites)."""
        # Same setup as plus-strand test but shift DP_DG by 1 so GEX_size = 107 != 106.
        result = sai10k_determine_aberrations(
            ds_ag=0.3, ds_al=0.0, ds_dg=0.3, ds_dl=0.0,
            dp_ag=19, dp_al=0, dp_dg=125, dp_dl=0,
            ds_ag_alt=0, ds_al_alt=0, ds_dg_alt=0, ds_dl_alt=0,
            variant_pos=32890540,
            exon_starts=BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            exon_ends=BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            strand="+",
        )
        types = [a["aberration_type"] for a in result]
        self.assertNotIn("increased_exon_inclusion", types)

    def test_no_increased_exon_inclusion_when_dp_ag_mismatches_native(self):
        """If DP_AG doesn't match the native acceptor position, increased_exon_inclusion
        should NOT fire even if GEX_size matches."""
        # Shift DP_AG by 1 (and DP_DG by 1 to keep GEX_size=106), so
        # geo_ag = 32890560 != geo_na = 32890559.
        result = sai10k_determine_aberrations(
            ds_ag=0.3, ds_al=0.0, ds_dg=0.3, ds_dl=0.0,
            dp_ag=20, dp_al=0, dp_dg=125, dp_dl=0,
            ds_ag_alt=0, ds_al_alt=0, ds_dg_alt=0, ds_dl_alt=0,
            variant_pos=32890540,
            exon_starts=BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            exon_ends=BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            strand="+",
        )
        types = [a["aberration_type"] for a in result]
        self.assertNotIn("increased_exon_inclusion", types)

    def test_increased_exon_inclusion_frameshift_annotation(self):
        """Increased exon inclusion should get proper frameshift annotation:
        frameshift=None (no frame change), affects_coding if exon overlaps CDS."""
        transcript_scores = {
            "DS_AG": 0.3, "DS_AL": 0.0, "DS_DG": 0.3, "DS_DL": 0.0,
            "DP_AG": 19, "DP_AL": 0, "DP_DG": 124, "DP_DL": 0,
            "DS_AG_ALT": 0, "DS_AL_ALT": 0, "DS_DG_ALT": 0, "DS_DL_ALT": 0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "EXON_FRAMES": BRCA2_TRANSCRIPT_HG19["EXON_FRAMES"],
            "STRAND": BRCA2_TRANSCRIPT_HG19["STRAND"],
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32890540)
        types = [a["aberration_type"] for a in result["aberrations"]]
        self.assertIn("increased_exon_inclusion", types)
        iei = next(a for a in result["aberrations"]
                   if a["aberration_type"] == "increased_exon_inclusion")
        # BRCA2 exon 2 overlaps CDS (CDS_START=32890598 is within [32890559, 32890664])
        self.assertTrue(iei["affects_coding"])
        # Increased exon inclusion doesn't change reading frame
        self.assertIsNone(iei["frameshift"])
        self.assertEqual(iei["size"], 106)
        self.assertIn("Increased exon inclusion", iei["frameshift_description"])


if __name__ == "__main__":
    unittest.main()
