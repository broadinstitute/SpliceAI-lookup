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
    _translate_dna,
    _reverse_complement,
    _find_difference_point,
    _format_aa_change_and_detect,
    _build_protein_windows,
    _apply_variant_to_altered_exon_seqs,
    _consensus_filtered_in_tx_order,
    _compute_consensus_context,
    _translate_spans,
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

    Mirrors the batched query pattern in server.get_transcript_structures
    (`WHERE transcript_id = ANY(%s)`) but unpacks the single-id result —
    server.py can't be imported here because its module-level startup runs
    GENOME_VERSION validation. The query SQL is kept in lockstep with
    production so an integration-test pass is meaningful evidence the
    production code path is also wiring up correctly.

    Args:
        conn: Database connection.
        transcript_id: Transcript ID (e.g., "ENST00000357654" or "NM_000059.4").
        genome_version: "37" or "38".

    Returns:
        dict with transcript structure fields, or None if not found.
    """
    if conn is None:
        return None

    transcript_id_without_version = transcript_id.split(".")[0]
    table_name = f"transcripts_hg{genome_version}"

    cursor = conn.cursor()
    cursor.execute(
        f"""SELECT transcript_id, strand, cds_start, cds_end, exon_starts, exon_ends
           FROM {table_name} WHERE transcript_id = ANY(%s)""",
        ([transcript_id_without_version],),
    )
    rows = cursor.fetchall()
    cursor.close()

    if not rows:
        return None

    _, strand, cds_start, cds_end, exon_starts_str, exon_ends_str = rows[0]

    # genePred uses 0-based half-open coordinates. Convert to 1-based closed.
    exon_starts_1based = [int(s) + 1 for s in exon_starts_str.rstrip(",").split(",") if s]
    exon_ends_1based = [int(s) for s in exon_ends_str.rstrip(",").split(",") if s]

    return {
        "EXON_STARTS": exon_starts_1based,
        "EXON_ENDS": exon_ends_1based,
        "CDS_START": cds_start + 1 if cds_start is not None else None,
        "CDS_END": cds_end if cds_end is not None else None,
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

    def test_intronic_variant_at_intron_midpoint_picks_upstream_exon(self):
        """Equidistant intronic variant: upstream exon wins on a tie.

        BRCA2 intron 5 (between genomic exons 5 and 6) spans 32900288..32900378
        (length 91). The integer midpoint 32900333 is exactly 46 bp from each
        flanking exon. With the `<=` tie-break in sai10k_find_affected_region,
        the variant is associated with the upstream exon's donor side, matching
        `_get_native_splice_sites`'s tie-break and avoiding cross-branch
        disagreement.
        """
        exon_starts = BRCA2_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA2_TRANSCRIPT_HG19["EXON_ENDS"]
        # Sanity-check distances are equal at this position.
        self.assertEqual(32900333 - exon_ends[4], exon_starts[5] - 32900333)
        result = sai10k_find_affected_region(32900333, exon_starts, exon_ends, strand="+")
        self.assertEqual(result["region_type"], "intron")
        self.assertEqual(result["nearest_boundary"], "donor")
        self.assertEqual(result["_genomic_index"], 4)


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

    def test_multi_site_cryptic_donor_gains(self):
        """Two independent cryptic donor-gain sites (PALB2 c.47A>G-like) should
        each yield a separate gain_B_donor aberration with its own alt_score,
        geo_dg and _partial_size."""
        # BRCA2 exon 2: [32890559, 32890664], native donor at 32890664 on '+'.
        # Variant at 32890650 (exonic, dp_nd = 14). Two donor-gain candidates:
        #   site 1 (max, at GEO_DG=32890640, dp_dg=-10): delta 0.20 / alt 0.51
        #           -> partial_size = (14 - (-10)) = 24 (partial exon deletion).
        #   site 2 (GEO_DG=32890630): delta 0.14 / alt 0.94
        #           -> partial_size = (14 - (-20)) = 34 (partial exon deletion).
        # ds_dl=0.3 / ds_dl_alt=0.4 let both pass _cryptic_gain_activation.
        all_non_zero_scores = [
            {"pos": "32890640", "RA": "0.00", "AA": "0.00", "RD": "0.31", "AD": "0.51"},
            {"pos": "32890630", "RA": "0.00", "AA": "0.00", "RD": "0.80", "AD": "0.94"},
        ]
        result = sai10k_determine_aberrations(
            ds_ag=0.0, ds_al=0.0, ds_dg=0.2, ds_dl=0.3,
            dp_ag=0, dp_al=0, dp_dg=-10, dp_dl=0,
            ds_ag_alt=0.0, ds_al_alt=0.0, ds_dg_alt=0.51, ds_dl_alt=0.4,
            variant_pos=32890650,
            exon_starts=BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            exon_ends=BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            strand="+",
            all_non_zero_scores=all_non_zero_scores,
        )
        donors = [a for a in result if a["_branch"] == "gain_B_donor"]
        self.assertEqual(len(donors), 2)
        # Each sibling carries its own alt_score, geo_dg and _partial_size.
        self.assertEqual({a["geo_dg"] for a in donors}, {32890640, 32890630})
        self.assertEqual({a["_partial_size"] for a in donors}, {24, 34})
        self.assertEqual({round(a["alt_score"], 2) for a in donors}, {0.51, 0.94})
        for a in donors:
            self.assertIsNotNone(a.get("alt_score"))

    def test_single_site_cryptic_donor_gain_regression(self):
        """Same setup as the multi-site case but with only the single strongest
        donor site (and also with all_non_zero_scores=None) must yield exactly
        ONE gain_B_donor aberration, matching pre-change behavior."""
        single_site = [
            {"pos": "32890640", "RA": "0.00", "AA": "0.00", "RD": "0.31", "AD": "0.51"},
        ]
        for scores in (single_site, None):
            result = sai10k_determine_aberrations(
                ds_ag=0.0, ds_al=0.0, ds_dg=0.2, ds_dl=0.3,
                dp_ag=0, dp_al=0, dp_dg=-10, dp_dl=0,
                ds_ag_alt=0.0, ds_al_alt=0.0, ds_dg_alt=0.51, ds_dl_alt=0.4,
                variant_pos=32890650,
                exon_starts=BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
                exon_ends=BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
                strand="+",
                all_non_zero_scores=scores,
            )
            donors = [a for a in result if a["_branch"] == "gain_B_donor"]
            self.assertEqual(len(donors), 1)
            self.assertEqual(donors[0]["aberration_type"], "partial_exon_deletion")
            self.assertEqual(donors[0]["geo_dg"], 32890640)
            self.assertEqual(donors[0]["_partial_size"], 24)
            self.assertEqual(round(donors[0]["alt_score"], 2), 0.51)

    def test_sibling_cryptic_donor_gain_boundary_delta(self):
        """Regression: a sibling donor gain whose per-site delta is exactly 0.20
        via two-decimal scores (AD='0.30', RD='0.10') must still be emitted. Naive
        float subtraction gives 0.19999999999999998, which fails the >= 0.2 branch-2
        threshold in _cryptic_gain_activation and would drop the site; the per-site
        delta is rounded to 2 decimals to avoid that."""
        all_non_zero_scores = [
            # strongest donor site (excluded from the sibling loop; emitted above)
            {"pos": "32890640", "RA": "0.00", "AA": "0.00", "RD": "0.31", "AD": "0.51"},
            # sibling with per-site delta exactly 0.20 (0.30 - 0.10) and alt 0.30,
            # so it only passes branch2 (>= 0.2) once the delta is rounded.
            {"pos": "32890630", "RA": "0.00", "AA": "0.00", "RD": "0.10", "AD": "0.30"},
        ]
        result = sai10k_determine_aberrations(
            ds_ag=0.0, ds_al=0.0, ds_dg=0.2, ds_dl=0.3,
            dp_ag=0, dp_al=0, dp_dg=-10, dp_dl=0,
            ds_ag_alt=0.0, ds_al_alt=0.0, ds_dg_alt=0.51, ds_dl_alt=0.4,
            variant_pos=32890650,
            exon_starts=BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            exon_ends=BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            strand="+",
            all_non_zero_scores=all_non_zero_scores,
        )
        donors = [a for a in result if a["_branch"] == "gain_B_donor"]
        self.assertEqual(len(donors), 2)
        boundary = [a for a in donors if a["geo_dg"] == 32890630]
        self.assertEqual(len(boundary), 1)
        self.assertEqual(round(boundary[0]["alt_score"], 2), 0.30)


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
            "STRAND": BRCA2_TRANSCRIPT_HG19["STRAND"],
            "TX_START": BRCA2_TRANSCRIPT_HG19["TX_START"],
            "TX_END": BRCA2_TRANSCRIPT_HG19["TX_END"],
        }

        result = sai10k_compute_predictions(transcript_scores, variant_pos=32889805)

        self.assertIn("aberrations", result)
        self.assertIn("transcript_info", result)
        self.assertEqual(result["transcript_info"]["strand"], "+")
        self.assertTrue(result["transcript_info"]["is_coding"])
        # A broken predictor returning {"aberrations": []} used to pass the old
        # assertions. BRCA2_c.-40+1G>A has DS_DL=0.99 and DS_DG=0.57, well above
        # the paper's DS_ALDL_MIN=0.02 / DS_AGDG_MIN=0.02 thresholds, so we
        # should always get at least one classified aberration.
        self.assertGreater(
            len(result["aberrations"]), 0,
            f"Expected at least one aberration for BRCA2 c.-40+1G>A, got: {result}",
        )
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
            "STRAND": BRCA1_TRANSCRIPT_HG19["STRAND"],
            "TX_START": BRCA1_TRANSCRIPT_HG19["TX_START"],
            "TX_END": BRCA1_TRANSCRIPT_HG19["TX_END"],
        }

        result = sai10k_compute_predictions(transcript_scores, variant_pos=41199659)

        self.assertIn("aberrations", result)
        self.assertEqual(result["transcript_info"]["strand"], "-")
        # Same rationale as BRCA2: DS_DL=0.93 plus DS_AL=0.40 and DS_DG=0.36
        # triples are well above threshold, so we should always classify at
        # least one aberration.
        self.assertGreater(
            len(result["aberrations"]), 0,
            f"Expected at least one aberration for BRCA1 c.5467+1G>A, got: {result}",
        )
        for aberration in result["aberrations"]:
            self.assertIn("frameshift", aberration)
            self.assertIn("affects_coding", aberration)


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

        self.assertIsNotNone(result, "BRCA2 (ENST00000380152) transcript missing from database")
        self.assertIn("EXON_STARTS", result)
        self.assertIn("EXON_ENDS", result)
        self.assertIn("CDS_START", result)
        self.assertIn("CDS_END", result)
        self.assertIn("STRAND", result)
        self.assertEqual(result["STRAND"], "+")

    def test_retrieve_transcript_structure_brca1(self):
        """Test retrieving BRCA1 transcript structure from database."""
        if self.conn is None:
            self.skipTest("Database connection not available")

        # Try to get BRCA1 transcript (ENST00000357654 is the canonical BRCA1)
        result = get_transcript_structure_from_db(self.conn, "ENST00000357654", "37")

        self.assertIsNotNone(result, "BRCA1 (ENST00000357654) transcript missing from database")
        self.assertIn("EXON_STARTS", result)
        self.assertIn("EXON_ENDS", result)
        self.assertIn("STRAND", result)
        self.assertEqual(result["STRAND"], "-")

    def test_predictions_with_db_transcript(self):
        """Test predictions using transcript data from database."""
        if self.conn is None:
            self.skipTest("Database connection not available")

        # Get BRCA1 transcript from database
        transcript_struct = get_transcript_structure_from_db(
            self.conn, "ENST00000357654", "37"
        )
        self.assertIsNotNone(transcript_struct, "BRCA1 (ENST00000357654) transcript missing from database")

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
        # BRCA1 c.5467+1G>A — strong DS_DL=0.93 must produce at least one aberration.
        # Mirrors the in-memory equivalent (TestExampleVariants) so a regression that
        # returns an empty list is not silently accepted.
        self.assertGreater(len(result["aberrations"]), 0)


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
            "STRAND": BRCA2_TRANSCRIPT_HG19["STRAND"],
        }

        result = sai10k_compute_predictions(transcript_scores, variant_pos=32890665)

        # With both DS_AL=0.81 and DS_DL=0.98 significant, loss-pair branch fires.
        self.assertEqual(len(result["aberrations"]), 1)
        self.assertEqual(result["aberrations"][0]["aberration_type"], "exon_skipping")

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
        frameshift=False for coding events (no frame change),
        frameshift=None for non-coding events, affects_coding if exon overlaps CDS."""
        transcript_scores = {
            "DS_AG": 0.3, "DS_AL": 0.0, "DS_DG": 0.3, "DS_DL": 0.0,
            "DP_AG": 19, "DP_AL": 0, "DP_DG": 124, "DP_DL": 0,
            "DS_AG_ALT": 0, "DS_AL_ALT": 0, "DS_DG_ALT": 0, "DS_DL_ALT": 0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
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
        self.assertFalse(iei["frameshift"])
        self.assertEqual(iei["size"], 106)
        self.assertEqual(iei["description"]["label"], "Increased exon inclusion")


class TestExonSkippingFallback(unittest.TestCase):
    """Fallback path in sai10k_annotate_frameshift for exon_skipping when the
    reference exact-match lookup (_find_exons_between_loss_positions) fails.

    The fallback declares single-exon skipping of the variant's affected exon
    whenever at least one of GEO_AL/GEO_DL matches that exon's own native
    splice site. Covers two reachable cases:

    * Exonic variants where SpliceAI's DP_AL/DP_DL peak lands interior to the
      affected exon (e.g. CHEK2 22:28695127 C>T, where DP_DL=0 anchors at the
      native donor but DP_AL points 50bp inside the exon).
    * Intronic splice-site variants where one peak lands at a native site and
      the other on a cryptic position (e.g. RAD51C c.404+2T>C, where DP_DL=-2
      anchors at exon 2's native donor but DP_AL points to a cryptic acceptor
      well inside the upstream exon).

    The four (strand, nearest_boundary) cells of the intronic case are each
    exercised by their own dedicated test below so the strand × side mapping
    in the fallback can't drift.
    """

    def test_fallback_uses_affected_exon_when_geo_al_interior(self):
        """Variant at exon 4's donor (BRCA2 + strand). DP_DL=0 anchors at the
        native donor, but DP_AL lands 50bp inside the exon so the reference
        exact-match fails — the fallback should recover exon 4."""
        # BRCA2 exon 4 (1-based 4, genomic index 3): [32899213, 32899321], 109 bp.
        # Variant at 32899321 (exon donor on + strand = last base of exon 4).
        # DP_DL = 0  -> GEO_DL = 32899321 == exon_ends[3]  (native donor match)
        # DP_AL = -50 -> GEO_AL = 32899271 (50bp into exon 4; doesn't match any boundary)
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.30, "DS_DG": 0.00, "DS_DL": 0.80,
            "DP_AG": 0, "DP_AL": -50, "DP_DG": 0, "DP_DL": 0,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32899321)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        # 109 bp exon 4, 109 % 3 = 1 -> frameshift.
        self.assertEqual(es[0]["size"], 109)
        self.assertTrue(es[0]["frameshift"])
        self.assertEqual(es[0]["description"], {
            "label": "Exon 4 skipping",
            "size_bp": 109,
            "size_is_coding": True,
            "consequence": None,
            "status": "frameshift",
            "introduces_stop_codon": False,
            "extends_past_native_stop": False,
        })

    def test_no_fallback_when_intronic_variant_has_no_native_anchor(self):
        """Intronic variant where neither GEO_AL nor GEO_DL lands on a native
        splice site of any flanking exon — fallback should NOT fire and the
        'could not be mapped' message should be emitted."""
        # BRCA2 intron 1 (between exon 1 end 32889804 and exon 2 start 32890559).
        # Variant 40bp into intron 1; DP_AL/DP_DL chosen so neither GEO_AL nor
        # GEO_DL hits any exon boundary.
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.30, "DS_DG": 0.00, "DS_DL": 0.80,
            "DP_AG": 0, "DP_AL": 100, "DP_DG": 0, "DP_DL": 200,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32889844)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        self.assertEqual(es[0]["description"], {
            "label": "Exon skipping predicted but skipped exon(s)",
            "size_bp": None,
            "size_is_coding": None,
            "consequence": None,
            "status": "could not be mapped",
            "introduces_stop_codon": False,
            "extends_past_native_stop": False,
        })

    def test_fallback_uses_affected_exon_when_intronic_variant_donor_at_native_plus_strand(self):
        """`+` strand intronic variant +2 from the native donor (RAD51C
        c.404+2T>C analogue). Exercises the `idx = intron_idx` branch on
        `+` strand."""
        # BRCA2 intron 4 (between exon 4 and exon 5): genomic_index 3,
        # spans 32899322..32900237. Variant at 32899323 (=exon_ends[3] + 2).
        # DP_DL = -2 -> GEO_DL = 32899321 == exon_ends[3]  (native donor of exon 4)
        # DP_AL = -100 -> GEO_AL = 32899223 (interior of exon 4, no native match)
        # nearest_boundary='donor', intron_idx=3, idx=3 -> exon 4 (1-based).
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.30, "DS_DG": 0.00, "DS_DL": 0.80,
            "DP_AG": 0, "DP_AL": -100, "DP_DG": 0, "DP_DL": -2,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32899323)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        # Exon 4 is fully coding (109 bp); 109 % 3 = 1 -> frameshift.
        self.assertEqual(es[0]["size"], 109)
        self.assertTrue(es[0]["frameshift"])
        self.assertEqual(es[0]["description"], {
            "label": "Exon 4 skipping",
            "size_bp": 109,
            "size_is_coding": True,
            "consequence": None,
            "status": "frameshift",
            "introduces_stop_codon": False,
            "extends_past_native_stop": False,
        })

    def test_fallback_uses_affected_exon_when_intronic_variant_acceptor_at_native_plus_strand(self):
        """`+` strand intronic variant -2 from the native acceptor of the
        downstream exon. Exercises the `idx = intron_idx + 1` branch on
        `+` strand."""
        # BRCA2 intron 3 (between exon 3 and exon 4): genomic_index 2,
        # spans 32893463..32899212. Variant at 32899211 (=exon_starts[3] - 2).
        # DP_AL = +2  -> GEO_AL = 32899213 == exon_starts[3]  (native acceptor of exon 4)
        # DP_DL = +100 -> GEO_DL = 32899311 (interior of exon 4, no native match)
        # nearest_boundary='acceptor', intron_idx=2, idx=3 -> exon 4 (1-based).
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.30, "DS_DG": 0.00, "DS_DL": 0.80,
            "DP_AG": 0, "DP_AL": 2, "DP_DG": 0, "DP_DL": 100,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32899211)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        self.assertEqual(es[0]["size"], 109)
        self.assertTrue(es[0]["frameshift"])
        self.assertEqual(es[0]["description"], {
            "label": "Exon 4 skipping",
            "size_bp": 109,
            "size_is_coding": True,
            "consequence": None,
            "status": "frameshift",
            "introduces_stop_codon": False,
            "extends_past_native_stop": False,
        })

    def test_fallback_uses_affected_exon_when_intronic_variant_acceptor_at_native_minus_strand(self):
        """`-` strand intronic variant near a transcript-acceptor (which sits
        on the genomic-upstream side of the intron). Exercises the
        `idx = intron_idx` branch on `-` strand."""
        # BRCA1 intron between genomic exons 5 and 6: genomic_index 5, spans
        # 41215391..41215890. Variant at 41215392 (=exon_ends[5] + 2; 2bp into
        # the intron from genomic exon 5, which is transcript exon 18).
        # On '-' strand the native ACCEPTOR of genomic exon 5 sits at
        # exon_ends[5] = 41215390.
        # DP_AL = -2  -> GEO_AL = 41215390 == exon_ends[5]  (native acceptor)
        # DP_DL = -50 -> GEO_DL = 41215342 (no native match)
        # nearest_boundary='acceptor' on '-' strand maps to idx = intron_idx,
        # i.e. genomic exon 5 = biological exon 23-5 = 18.
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.30, "DS_DG": 0.00, "DS_DL": 0.80,
            "DP_AG": 0, "DP_AL": -2, "DP_DG": 0, "DP_DL": -50,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA1_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA1_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "-",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=41215392)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        # Genomic exon 5: [41215350, 41215390], 41 bp, fully coding. 41 % 3 = 2 -> frameshift.
        self.assertEqual(es[0]["size"], 41)
        self.assertTrue(es[0]["frameshift"])
        self.assertEqual(es[0]["description"], {
            "label": "Exon 18 skipping",
            "size_bp": 41,
            "size_is_coding": True,
            "consequence": None,
            "status": "frameshift",
            "introduces_stop_codon": False,
            "extends_past_native_stop": False,
        })

    def test_fallback_uses_affected_exon_when_intronic_variant_donor_at_native_minus_strand(self):
        """`-` strand intronic variant near a transcript-donor (which sits on
        the genomic-downstream side of the intron). Exercises the
        `idx = intron_idx + 1` branch on `-` strand."""
        # Same BRCA1 intron (genomic_index 5, 41215391..41215890). Variant at
        # 41215889 (=exon_starts[6] - 2). On '-' strand the native DONOR of
        # genomic exon 6 sits at exon_starts[6] = 41215891.
        # DP_DL = +2  -> GEO_DL = 41215891 == exon_starts[6]  (native donor)
        # DP_AL = +50 -> GEO_AL = 41215939 (interior of exon 6, no native match)
        # nearest_boundary='donor' on '-' strand maps to idx = intron_idx + 1,
        # i.e. genomic exon 6 = biological exon 23-6 = 17.
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.30, "DS_DG": 0.00, "DS_DL": 0.80,
            "DP_AG": 0, "DP_AL": 50, "DP_DG": 0, "DP_DL": 2,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA1_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA1_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "-",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=41215889)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        # Genomic exon 6: [41215891, 41215968], 78 bp, fully coding. 78 % 3 = 0 -> in-frame.
        self.assertEqual(es[0]["size"], 78)
        self.assertFalse(es[0]["frameshift"])
        self.assertEqual(es[0]["description"], {
            "label": "Exon 17 skipping",
            "size_bp": 78,
            "size_is_coding": True,
            "consequence": None,
            "status": "in-frame",
            "introduces_stop_codon": False,
            "extends_past_native_stop": False,
        })


class TestDeltaType(unittest.TestCase):
    """delta_type (name of the Δ column carrying the overall max Δ score) is
    computed per-variant by sai10k_compute_predictions and copied onto every
    aberration."""

    def _base_scores(self):
        return {
            "DS_AG": 0.00, "DS_AL": 0.00, "DS_DG": 0.00, "DS_DL": 0.00,
            "DP_AG": 0, "DP_AL": 0, "DP_DG": 0, "DP_DL": 0,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }

    def test_donor_loss_delta_type(self):
        """DS_DL is the max → delta_type = 'donor loss'. BRCA2 c.67+1G>A."""
        scores = self._base_scores()
        scores.update({"DS_AL": 0.81, "DS_DL": 0.98, "DP_AL": -106, "DP_DL": -1})
        result = sai10k_compute_predictions(scores, variant_pos=32890665)
        self.assertTrue(result["aberrations"])
        self.assertEqual(result["aberrations"][0]["delta_type"], "donor loss")

    def test_acceptor_loss_delta_type(self):
        """DS_AL is the max → delta_type = 'acceptor loss'."""
        scores = self._base_scores()
        # Variant at BRCA2 exon 2 donor (32890664); DP_DL=0, DP_AL matches exon 2 acceptor.
        scores.update({
            "DS_AL": 0.54, "DS_DL": 0.32,
            "DP_AL": -105,  # 32890664 + (-105) = 32890559 = exon_starts[1]
            "DP_DL": 0,
        })
        result = sai10k_compute_predictions(scores, variant_pos=32890664)
        self.assertTrue(result["aberrations"])
        self.assertEqual(result["aberrations"][0]["delta_type"], "acceptor loss")

    def test_acceptor_gain_delta_type(self):
        """DS_AG is the max → delta_type = 'acceptor gain'. Cryptic-acceptor
        partial_exon_deletion fixture."""
        scores = self._base_scores()
        scores.update({
            "DS_AG": 0.30, "DP_AG": 20,
            "DS_AG_ALT": 0.5, "DS_AL_ALT": 0.3,
        })
        result = sai10k_compute_predictions(scores, variant_pos=32890580)
        self.assertTrue(result["aberrations"])
        self.assertEqual(result["aberrations"][0]["delta_type"], "acceptor gain")

    def test_all_aberrations_share_variant_level_delta_type(self):
        """When one variant produces multiple aberrations (e.g. combination
        row), they should all carry the same delta_type since it's variant-level,
        not per-aberration."""
        scores = self._base_scores()
        # BRCA2 c.-40+1G>A setup: DS_DL=0.99 loss + cryptic donor partial combo.
        scores.update({
            "DS_AG": 0.01, "DS_AL": 0.08, "DS_DG": 0.57, "DS_DL": 0.99,
            "DP_AG": -38, "DP_AL": 754, "DP_DG": -100, "DP_DL": -1,
            "DS_DG_ALT": 0.7, "DS_DL_ALT": 0.01,
        })
        result = sai10k_compute_predictions(scores, variant_pos=32889805)
        self.assertGreaterEqual(len(result["aberrations"]), 2)
        delta_types = {a["delta_type"] for a in result["aberrations"]}
        self.assertEqual(delta_types, {"donor loss"})


class TestPotentialPrefix(unittest.TestCase):
    """'Potential ' prefix on exon_skipping / whole_intron_retention when the
    weaker of the two loss Δ scores is in [0.02, 0.2) — flags loss-pair
    asymmetry to the reader (Canson reviewer feedback, 2026-04-17)."""

    def test_potential_exon_skipping_on_weak_paired_loss(self):
        """BRCA2 Figure 8 profile: DS_AL=0.13 weak partner + DS_DL=0.93 strong.
        Expect description.label = 'Potential exon N skipping'."""
        # BRCA2 exon 2 (1-based 2): [32890559, 32890664], 106 bp, frameshift.
        # Variant at 32890600 inside exon 2. DP_AL=-41 matches exon_starts[1],
        # DP_DL=64 matches exon_ends[1], so reference lookup succeeds directly.
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.13, "DS_DG": 0.00, "DS_DL": 0.93,
            "DP_AG": 0, "DP_AL": -41, "DP_DG": 0, "DP_DL": 64,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32890600)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        # 106 bp, 106 % 3 = 1 -> frameshift.
        self.assertEqual(es[0]["description"], {
            "label": "Potential exon 2 skipping",
            "size_bp": 67,
            "size_is_coding": True,
            "consequence": None,
            "status": "start codon lost",
            "introduces_stop_codon": False,
            "extends_past_native_stop": False,
        })

    def test_no_potential_prefix_when_both_loss_scores_strong(self):
        """Both DS_AL and DS_DL ≥ 0.2 → no 'Potential' prefix."""
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.50, "DS_DG": 0.00, "DS_DL": 0.93,
            "DP_AG": 0, "DP_AL": -41, "DP_DG": 0, "DP_DL": 64,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32890600)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        self.assertFalse(es[0]["description"]["label"].startswith("Potential"))
        self.assertEqual(es[0]["description"], {
            "label": "Exon 2 skipping",
            "size_bp": 67,
            "size_is_coding": True,
            "consequence": None,
            "status": "start codon lost",
            "introduces_stop_codon": False,
            "extends_past_native_stop": False,
        })


class TestPartialShiftDescriptions(unittest.TestCase):
    """Reformatted line-2 strings for partial_exon_deletion / partial_intron_retention:
    '<Donor|Acceptor> shift (Nbp partial ...) - <frameshift|in-frame|non-coding>'."""

    def test_acceptor_shift_partial_exon_deletion(self):
        """Cryptic acceptor activation inside an exon → 'Acceptor shift (...)'."""
        # Reuse the cryptic_acceptor fixture at BRCA2 exon 2.
        transcript_scores = {
            "DS_AG": 0.30, "DS_AL": 0.00, "DS_DG": 0.00, "DS_DL": 0.00,
            "DP_AG": 20, "DP_AL": 0, "DP_DG": 0, "DP_DL": 0,
            "DS_AG_ALT": 0.5, "DS_AL_ALT": 0.3, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32890580)
        pd = [a for a in result["aberrations"] if a["aberration_type"] == "partial_exon_deletion"]
        self.assertEqual(len(pd), 1)
        # partial_size = dp_ag - dp_na = 20 - (-21) = 41 bp; 41 % 3 = 2 -> frameshift.
        self.assertEqual(pd[0]["size"], 41)
        self.assertEqual(pd[0]["description"], {
            "label": "Acceptor shift",
            "size_bp": 2,
            "size_is_coding": True,
            "consequence": "partial exon 2 deletion",
            "status": "start codon lost",
            "introduces_stop_codon": False,
            "extends_past_native_stop": False,
        })

    def test_acceptor_shift_partial_intron_retention(self):
        """gain_B_acceptor on a + strand transcript: cryptic acceptor in the
        intron upstream of the affected biological exon, so the retained intron
        number = bio_exon - 1."""
        # BRCA2 (+ strand): variant just before exon 3 acceptor (32893214) with
        # a cryptic acceptor 50bp upstream → partial_intron_retention on intron 2.
        transcript_scores = {
            "DS_AG": 0.80, "DS_AL": 0.30, "DS_DG": 0.00, "DS_DL": 0.00,
            "DS_AG_ALT": 0.95, "DS_AL_ALT": 0.05, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "DP_AG": -50, "DP_AL": 1, "DP_DG": 0, "DP_DL": 0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32893213)
        pir = [a for a in result["aberrations"] if a["aberration_type"] == "partial_intron_retention"]
        self.assertEqual(len(pir), 1)
        self.assertEqual(pir[0]["affected_intron_number"], 2)
        self.assertEqual(pir[0]["description"], {
            "label": "Acceptor shift",
            "size_bp": 51,
            "size_is_coding": True,
            "consequence": "partial intron 2 retention",
            "status": "in-frame",
            "introduces_stop_codon": False,
            "extends_past_native_stop": False,
        })

    def test_donor_shift_partial_intron_retention(self):
        """gain_B_donor on a - strand transcript: cryptic donor in the intron
        downstream (in transcript order) of the affected biological exon, so
        the retained intron number = bio_exon."""
        # BRCA1 c.5467+1G>A fixture (- strand): variant at 41199659 produces
        # a partial_intron_retention on intron 22 via the donor branch.
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.40, "DS_DG": 0.36, "DS_DL": 0.93,
            "DS_AG_ALT": 0, "DS_AL_ALT": 0,
            "DS_DG_ALT": 0.95, "DS_DL_ALT": 0.05,
            "DP_AG": -4, "DP_AL": 61, "DP_DG": -4, "DP_DL": 1,
            "EXON_STARTS": BRCA1_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA1_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": BRCA1_TRANSCRIPT_HG19["STRAND"],
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=41199659)
        pir = [a for a in result["aberrations"] if a["aberration_type"] == "partial_intron_retention"]
        self.assertEqual(len(pir), 1)
        self.assertEqual(pir[0]["affected_intron_number"], 22)
        self.assertEqual(pir[0]["description"], {
            "label": "Donor shift",
            "size_bp": 5,
            "size_is_coding": True,
            "consequence": "partial intron 22 retention",
            "status": "frameshift",
            "introduces_stop_codon": False,
            "extends_past_native_stop": False,
        })


class TestStopCodonLost(unittest.TestCase):
    """Skipping or partially deleting an exon that contains cds_end should
    label the aberration 'stop codon lost' (priority over frameshift/in-frame)."""

    def test_exon_skipping_last_exon_contains_cds_end(self):
        """Skip BRCA2 exon 27 (last exon, 32972299-32974405; contains CDS_END=32972907).
        Loss peaks at exon 27 acceptor (DP_AL=0) and exon 27 donor (DP_DL=2106)."""
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.80, "DS_DG": 0.00, "DS_DL": 0.50,
            "DP_AG": 0, "DP_AL": 0, "DP_DG": 0, "DP_DL": 2106,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32972299)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        self.assertTrue(es[0]["stop_codon_lost"])
        self.assertFalse(es[0]["start_codon_lost"])
        self.assertIsNone(es[0]["frameshift"])  # suppressed when codon is lost
        self.assertEqual(es[0]["description"]["status"], "stop codon lost")
        self.assertTrue(es[0]["description"]["size_is_coding"])

    def test_no_partial_exon_deletion_in_terminal_exon(self):
        """Donor-side cryptic shift inside BRCA2 exon 27 (the last exon) must NOT
        produce partial_exon_deletion: there is no biological next exon, so the
        R reference's Cryptic_Donor_orientation is NA and the partial event is
        suppressed (spliceAI_parser.R:954, 963-968)."""
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.00, "DS_DG": 0.30, "DS_DL": 0.00,
            "DP_AG": 0, "DP_AL": 0, "DP_DG": 300, "DP_DL": 0,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.5, "DS_DL_ALT": 0.3,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32972500)
        pd = [a for a in result["aberrations"] if a["aberration_type"] == "partial_exon_deletion"]
        self.assertEqual(len(pd), 0)


class TestPseudoexonBoundaryMatchesRParser(unittest.TestCase):
    """Pin that the pseudoexon AG/DG distance-from-intron-boundary threshold
    matches spliceAI_parser.R lines 211-214. R uses `intronStart + 50 < AGpos`
    where intronStart = prev_eEnd + 1, which for integer positions means
    `AGpos - prev_eEnd > 51`. This implementation uses > 51 bp distance to
    match R's effective threshold; a gain at exactly 51 bp from a flanking
    boundary must be REJECTED (previously Python accepted, R rejected)."""

    # Shared synthetic transcript fixture for the two boundary tests.
    #   exon 1: [1, 100]        (donor at 100)
    #   exon 2: [1001, 1100]
    #   exon 3: [2001, 2100]
    # Intron 1 = [101, 1000] — plenty of room to place both AG and DG.
    FIXTURE_EXON_STARTS = [1, 1001, 2001]
    FIXTURE_EXON_ENDS = [100, 1100, 2100]

    def _run_with_ag_at(self, ag_dist_from_upstream):
        """Run Subflow A with AG placed `ag_dist_from_upstream` bp past the upstream
        exon's donor, DG 200 bp further into the intron, and a variant halfway
        between them. Returns the list of pseudoexon aberrations emitted."""
        # Place variant at pos 250 (well inside intron 1). Upstream exon ends at 100.
        variant_pos = 250
        geo_ag = 100 + ag_dist_from_upstream
        geo_dg = geo_ag + 200    # gex_size = 201 bp, safely within [25, 500]
        transcript_scores = {
            "DS_AG": 0.30, "DS_AL": 0.00, "DS_DG": 0.40, "DS_DL": 0.00,
            "DP_AG": geo_ag - variant_pos,
            "DP_AL": 0,
            "DP_DG": geo_dg - variant_pos,
            "DP_DL": 0,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": self.FIXTURE_EXON_STARTS,
            "EXON_ENDS": self.FIXTURE_EXON_ENDS,
            "CDS_START": 1,
            "CDS_END": 2100,
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        return [a for a in result["aberrations"] if a["aberration_type"] == "pseudoexon"]

    def test_ag_at_51bp_from_upstream_boundary_rejected(self):
        """AG exactly 51 bp downstream of the upstream exon's donor: must NOT
        fire pseudoexon (matches R's `intronStart + 50 < AGpos` at the boundary)."""
        self.assertEqual(len(self._run_with_ag_at(51)), 0,
                         "AG at exactly 51 bp from upstream exon boundary must be rejected (matches R)")

    def test_ag_at_52bp_from_upstream_boundary_accepted(self):
        """AG at 52 bp (one more) must FIRE pseudoexon."""
        self.assertEqual(len(self._run_with_ag_at(52)), 1,
                         "AG at 52 bp from upstream exon boundary must fire pseudoexon")

    # --- Mirror the + strand boundary pinning on the − strand ---
    #
    # On − strand, the biologically upstream exon is at HIGHER genomic coordinate
    # than the gained exon; "AG 51 bp from the upstream exon's donor" therefore
    # corresponds to AG being 51 bp BELOW the genomically-next exon's start.
    MINUS_STRAND_EXON_STARTS = [1, 1001, 2001]
    MINUS_STRAND_EXON_ENDS = [100, 1100, 2100]

    def _run_minus_strand_with_ag_at(self, ag_dist_from_upstream):
        """Same synthetic transcript but strand='-'. Places AG at
        `ag_dist_from_upstream` bp from the biologically-upstream exon's donor
        (genomically downstream exon's start, i.e. below exon_starts[2] = 2001)
        and DG 200 bp deeper into the intron in transcript order."""
        variant_pos = 1500
        # "biologically upstream donor" on − strand = exon_starts[2] = 2001.
        geo_ag = 2001 - ag_dist_from_upstream   # below genomically-next exon start
        geo_dg = geo_ag - 200                    # deeper in transcript order (lower coord)
        transcript_scores = {
            "DS_AG": 0.30, "DS_AL": 0.00, "DS_DG": 0.40, "DS_DL": 0.00,
            "DP_AG": geo_ag - variant_pos,
            "DP_AL": 0,
            "DP_DG": geo_dg - variant_pos,
            "DP_DL": 0,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": self.MINUS_STRAND_EXON_STARTS,
            "EXON_ENDS": self.MINUS_STRAND_EXON_ENDS,
            "CDS_START": 1,
            "CDS_END": 2100,
            "STRAND": "-",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        return [a for a in result["aberrations"] if a["aberration_type"] == "pseudoexon"]

    def test_minus_strand_ag_at_51bp_rejected(self):
        """− strand: AG exactly 51 bp from the biologically-upstream exon boundary must be rejected."""
        self.assertEqual(len(self._run_minus_strand_with_ag_at(51)), 0,
                         "− strand AG at 51 bp from upstream boundary must be rejected (matches R)")

    def test_minus_strand_ag_at_52bp_accepted(self):
        """− strand: AG at 52 bp must fire pseudoexon."""
        self.assertEqual(len(self._run_minus_strand_with_ag_at(52)), 1,
                         "− strand AG at 52 bp must fire pseudoexon")

    # --- Symmetric DG-from-downstream-boundary boundary tests on + strand ---
    #
    # Pseudoexon activation also requires DG to be > 51 bp from the downstream
    # exon boundary of its host intron. R's line 211 uses the symmetric condition
    # `intronEndAdj − 50 > DGpos`. Test DG at 51 (reject) and 52 (accept).

    def _run_with_dg_at(self, dg_dist_from_downstream):
        """Place DG `dg_dist_from_downstream` bp before the downstream exon's
        acceptor (= 1001 in the synthetic fixture), with AG fixed 200 bp further
        upstream so gex_size stays in [25, 500]."""
        variant_pos = 500
        geo_dg = 1001 - dg_dist_from_downstream
        geo_ag = geo_dg - 200
        transcript_scores = {
            "DS_AG": 0.30, "DS_AL": 0.00, "DS_DG": 0.40, "DS_DL": 0.00,
            "DP_AG": geo_ag - variant_pos,
            "DP_AL": 0,
            "DP_DG": geo_dg - variant_pos,
            "DP_DL": 0,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": self.FIXTURE_EXON_STARTS,
            "EXON_ENDS": self.FIXTURE_EXON_ENDS,
            "CDS_START": 1,
            "CDS_END": 2100,
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        return [a for a in result["aberrations"] if a["aberration_type"] == "pseudoexon"]

    def test_dg_at_51bp_from_downstream_boundary_rejected(self):
        """DG exactly 51 bp before the downstream exon's acceptor must NOT fire pseudoexon."""
        self.assertEqual(len(self._run_with_dg_at(51)), 0,
                         "DG at exactly 51 bp from downstream exon boundary must be rejected (matches R)")

    def test_dg_at_52bp_from_downstream_boundary_accepted(self):
        """DG at 52 bp must fire pseudoexon."""
        self.assertEqual(len(self._run_with_dg_at(52)), 1,
                         "DG at 52 bp from downstream exon boundary must fire pseudoexon")


class TestMinusStrandCodonLost(unittest.TestCase):
    """On − strand genes, cds_start is the GENOMIC lower CDS bound, which is
    the biological STOP codon (not start); similarly cds_end is the biological
    START codon. These tests pin that start_codon_lost / stop_codon_lost are
    computed biologically, not genomically."""

    def test_minus_strand_skip_exon_containing_biological_atg(self):
        """BRCA1 is − strand. Biological ATG lies at genomic CDS_END = 41276113,
        which falls inside the exon at genomic index 21 (biological exon 2):
        [41276034, 41276132]. Skipping it must report 'start codon lost'."""
        exon_starts = BRCA1_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA1_TRANSCRIPT_HG19["EXON_ENDS"]
        idx = 21  # genomic index of the exon containing CDS_END
        # Place variant at the biological acceptor (exon_ends[idx]). Loss peaks
        # at biological acceptor and donor (exon_starts[idx]) of this exon.
        variant_pos = exon_ends[idx]
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.80, "DS_DG": 0.00, "DS_DL": 0.50,
            "DP_AG": 0,
            "DP_AL": 0,                           # GEO_AL = exon_ends[idx]
            "DP_DG": 0,
            "DP_DL": exon_starts[idx] - variant_pos,   # GEO_DL = exon_starts[idx]
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": exon_starts,
            "EXON_ENDS": exon_ends,
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "-",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        self.assertTrue(es[0]["start_codon_lost"],
                        "skipping the exon containing the biological ATG must set start_codon_lost=True on − strand")
        self.assertFalse(es[0]["stop_codon_lost"])
        self.assertEqual(es[0]["description"]["status"], "start codon lost")

    def test_minus_strand_skip_exon_containing_biological_stop(self):
        """Biological stop codon lies at genomic CDS_START = 41197695, inside
        the first exon (genomic index 0 = biological last exon on − strand):
        [41196312, 41197819]. Skipping it must report 'stop codon lost'."""
        exon_starts = BRCA1_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA1_TRANSCRIPT_HG19["EXON_ENDS"]
        idx = 0
        # On − strand: biological acceptor = exon_ends[idx], biological donor = exon_starts[idx].
        # For exon skipping on − strand, orientation check is DP_AL*strand < DP_DL*strand,
        # i.e. DP_AL > DP_DL (strand = −1). Place the variant at the biological acceptor.
        variant_pos = exon_ends[idx]
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.80, "DS_DG": 0.00, "DS_DL": 0.50,
            "DP_AG": 0,
            "DP_AL": 0,                                # GEO_AL = exon_ends[idx]
            "DP_DG": 0,
            "DP_DL": exon_starts[idx] - variant_pos,   # GEO_DL = exon_starts[idx]
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": exon_starts,
            "EXON_ENDS": exon_ends,
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "-",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        self.assertTrue(es[0]["stop_codon_lost"],
                        "skipping the exon containing the biological stop codon must set stop_codon_lost=True on − strand")
        self.assertFalse(es[0]["start_codon_lost"])
        self.assertEqual(es[0]["description"]["status"], "stop codon lost")


class TestCodonInMiddleExonOfMultiExonSkip(unittest.TestCase):
    """When a multi-exon skip brackets ≥3 exons and the codon lies in a
    non-boundary exon, the iteration over `genomic_indices` must still detect
    codon loss. Pins that the `any(contains_* for i in ...)` loop reaches
    interior exons, not just the two boundaries."""

    def test_start_codon_in_middle_of_three_exon_skip(self):
        # Synthetic transcript with 5 exons. Skip exons 2-4 (genomic idx 1-3).
        # ATG sits inside exon 3 (middle of the skipped bracket).
        exon_starts = [1, 201, 401, 601, 801]
        exon_ends = [100, 300, 500, 700, 900]
        cds_start = 450           # inside exon 3 → ATG = [450, 452]
        cds_end = 880             # inside exon 5
        # Loss peaks: acceptor of exon 2 (= 201), donor of exon 4 (= 700).
        variant_pos = 201
        transcript_scores = {
            "DS_AG": 0.0, "DS_AL": 0.8, "DS_DG": 0.0, "DS_DL": 0.8,
            "DP_AG": 0, "DP_AL": 0, "DP_DG": 0, "DP_DL": 700 - 201,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": exon_starts, "EXON_ENDS": exon_ends,
            "CDS_START": cds_start, "CDS_END": cds_end, "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        # The bracket is exons 2-4; ATG is in exon 3 (middle). The loop must
        # reach the middle exon, not just the boundary ones.
        self.assertTrue(es[0]["start_codon_lost"])
        self.assertFalse(es[0]["stop_codon_lost"])
        self.assertEqual(es[0]["description"]["status"], "start codon lost")


class TestCodonAtExonEndpoint(unittest.TestCase):
    """Pin that contains_start_codon / contains_stop_codon use inclusive
    bounds: a codon position exactly AT an exon endpoint must count as
    contained. Guards against accidental strict-inequality regression."""

    def test_start_codon_at_exon_start_triggers_codon_lost(self):
        # cds_start exactly at exon_starts[1] — ATG spans [201, 203].
        exon_starts = [1, 201, 401]
        exon_ends = [100, 300, 500]
        transcript_scores = {
            "DS_AG": 0.0, "DS_AL": 0.8, "DS_DG": 0.0, "DS_DL": 0.8,
            "DP_AG": 0, "DP_AL": 0, "DP_DG": 0, "DP_DL": 300 - 201,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": exon_starts, "EXON_ENDS": exon_ends,
            "CDS_START": 201,     # = exon_starts[1] (boundary)
            "CDS_END": 480,
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=201)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        self.assertTrue(es[0]["start_codon_lost"],
                        "codon at exact exon boundary must count as contained")


class TestPartialExonDeletionMinusStrandCodonLost(unittest.TestCase):
    """Mirror of TestPartialShiftDescriptions + TestStopCodonLost but on −
    strand. The cryptic segment bounds + codon-position swap on − strand is
    a recently-added feature; these tests pin it."""

    def test_minus_strand_acceptor_shift_loses_start_codon(self):
        """BRCA1 − strand. Biological ATG at genomic CDS_END = 41276113,
        inside exon 21 [41276034, 41276132]. Variant at that exon's biological
        acceptor (= exon_ends[21] = 41276132). Cryptic acceptor shift INTO
        the exon: GEO_AG < 41276132 (lower genomic coord on − strand).
        Set GEO_AG at 41276100 (32 bp into the exon). The deletion segment
        on − strand acceptor shift = [GEO_AG + 1, native_acceptor] =
        [41276101, 41276132], which contains CDS_END = 41276113."""
        variant_pos = 41276132
        transcript_scores = {
            "DS_AG": 0.30, "DS_AL": 0.00, "DS_DG": 0.00, "DS_DL": 0.00,
            "DP_AG": 41276100 - variant_pos,   # = -32
            "DP_AL": 0, "DP_DG": 0, "DP_DL": 0,
            "DS_AG_ALT": 0.5, "DS_AL_ALT": 0.3, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA1_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA1_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "-",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        ped = [a for a in result["aberrations"] if a["aberration_type"] == "partial_exon_deletion"]
        self.assertEqual(len(ped), 1)
        self.assertTrue(ped[0]["start_codon_lost"],
                        "− strand partial_exon_deletion crossing CDS_END must report start codon lost")
        self.assertFalse(ped[0]["stop_codon_lost"])
        self.assertEqual(ped[0]["description"]["status"], "start codon lost")


class TestPartialExonDeletionMinusStrandMoreCases(unittest.TestCase):
    """Mirror the + strand partial_exon_deletion codon-lost tests across the
    remaining three − strand combinations:
      acceptor shift × stop codon,
      donor    shift × start codon,
      donor    shift × stop  codon.
    Only the acceptor × start case is pinned by TestPartialExonDeletionMinusStrandCodonLost."""

    def test_minus_strand_donor_shift_in_terminal_exon_suppressed(self):
        """BRCA1 − strand. exon_starts[0] = 41196312 is the biologically LAST
        exon (− strand reverses biological order). A donor shift in the last
        biological exon has no downstream exon to align with, so the R
        reference's Cryptic_Donor_orientation is NA (next_eEnd is NA) and the
        partial event is suppressed (spliceAI_parser.R:954, 963-968)."""
        variant_pos = 41196312
        exon_starts = BRCA1_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA1_TRANSCRIPT_HG19["EXON_ENDS"]
        geo_dg = 41197700
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.00, "DS_DG": 0.30, "DS_DL": 0.00,
            "DP_AG": 0, "DP_AL": 0,
            "DP_DG": geo_dg - variant_pos,
            "DP_DL": 0,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.5, "DS_DL_ALT": 0.3,
            "EXON_STARTS": exon_starts, "EXON_ENDS": exon_ends,
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "-",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        ped = [a for a in result["aberrations"] if a["aberration_type"] == "partial_exon_deletion"]
        self.assertEqual(len(ped), 0)


class TestUtrOnlyExonSkipping(unittest.TestCase):
    """Skipping an exon that lies entirely in the 5'UTR or 3'UTR should not
    affect the CDS at all: cds_size == 0, affects_coding == False,
    start_codon_lost / stop_codon_lost == None (non-coding branch), and label
    'non-coding' regardless of strand."""

    def test_plus_strand_skip_5utr_exon(self):
        # Synthetic transcript: exons [1,100], [201,300], [401,500].
        # CDS starts in exon 2 (cds_start=220). Exon 1 is entirely 5'UTR.
        # Skip exon 1 by pointing GEO_AL at its acceptor and GEO_DL at its donor.
        variant_pos = 1
        transcript_scores = {
            "DS_AG": 0.0, "DS_AL": 0.8, "DS_DG": 0.0, "DS_DL": 0.8,
            "DP_AG": 0, "DP_AL": 0, "DP_DG": 0, "DP_DL": 99,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": [1, 201, 401], "EXON_ENDS": [100, 300, 500],
            "CDS_START": 220, "CDS_END": 480, "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        self.assertEqual(es[0]["cds_size"], 0)
        self.assertFalse(es[0]["affects_coding"])
        self.assertIsNone(es[0]["start_codon_lost"])
        self.assertIsNone(es[0]["stop_codon_lost"])
        self.assertIsNone(es[0]["frameshift"])
        self.assertEqual(es[0]["description"]["status"], "non-coding change")


class TestPseudoexonNeverLosesCodon(unittest.TestCase):
    """pseudoexon adds bases; start_codon_lost and stop_codon_lost must be
    False when the aberration overlaps the CDS, and None otherwise."""

    def test_pseudoexon_coding_overlap_flags_are_false(self):
        # Use the TestPseudoexonBoundaryMatchesRParser synthetic transcript that
        # does fire a pseudoexon. CDS_START=1 and CDS_END=2100 so the
        # pseudoexon inside intron 1 lies entirely within the CDS range.
        variant_pos = 250
        transcript_scores = {
            "DS_AG": 0.30, "DS_AL": 0.00, "DS_DG": 0.40, "DS_DL": 0.00,
            "DP_AG": (100 + 52) - variant_pos,
            "DP_AL": 0,
            "DP_DG": (100 + 52 + 200) - variant_pos,
            "DP_DL": 0,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": [1, 1001, 2001],
            "EXON_ENDS": [100, 1100, 2100],
            "CDS_START": 1, "CDS_END": 2100, "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        pseudo = [a for a in result["aberrations"] if a["aberration_type"] == "pseudoexon"]
        self.assertEqual(len(pseudo), 1)
        self.assertFalse(pseudo[0]["start_codon_lost"])
        self.assertFalse(pseudo[0]["stop_codon_lost"])


class TestIncreasedExonInclusionNeverLosesCodon(unittest.TestCase):
    """increased_exon_inclusion includes a native exon verbatim; it cannot
    delete any codon."""

    def test_iei_coding_overlap_flags_are_false(self):
        # Synthetic 3-exon transcript on + strand. Exon 2 is entirely in CDS.
        # Variant at exon 2 acceptor, AG at native acceptor (DP_AG = DP_NA),
        # DG at native donor (DP_DG = DP_ND) so IEI fires.
        # exon 2 = [201, 300]; variant at 201; DP_AG = 0; DP_DG = 99.
        variant_pos = 201
        transcript_scores = {
            "DS_AG": 0.30, "DS_AL": 0.00, "DS_DG": 0.40, "DS_DL": 0.00,
            "DP_AG": 0, "DP_AL": 0, "DP_DG": 99, "DP_DL": 0,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": [1, 201, 401],
            "EXON_ENDS": [100, 300, 500],
            "CDS_START": 50, "CDS_END": 480, "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        iei = [a for a in result["aberrations"] if a["aberration_type"] == "increased_exon_inclusion"]
        self.assertEqual(len(iei), 1)
        self.assertFalse(iei[0]["start_codon_lost"])
        self.assertFalse(iei[0]["stop_codon_lost"])


class TestEntireCdsInOneExon(unittest.TestCase):
    """When a single skipped exon contains both cds_start and cds_end (e.g. a
    single-coding-exon gene), priority rules: start_codon_lost=True wins,
    stop_codon_lost=False, frameshift=None, label='start codon lost'."""

    def test_skipped_exon_contains_both_start_and_stop(self):
        # Synthetic transcript: exons [1,100], [201,400], [501,600].
        # Entire CDS [250, 350] lies inside exon 2.
        # Skip exon 2 by setting GEO_AL at exon 2 acceptor (201) and GEO_DL at
        # exon 2 donor (400). Variant at exon 2 acceptor.
        variant_pos = 201
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.80, "DS_DG": 0.00, "DS_DL": 0.50,
            "DP_AG": 0, "DP_AL": 0, "DP_DG": 0, "DP_DL": 199,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": [1, 201, 501],
            "EXON_ENDS": [100, 400, 600],
            "CDS_START": 250,
            "CDS_END": 350,
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        # Contract: at most one codon-lost flag is True (priority: start > stop).
        self.assertTrue(es[0]["start_codon_lost"])
        self.assertFalse(es[0]["stop_codon_lost"])
        self.assertIsNone(es[0]["frameshift"])
        self.assertEqual(es[0]["description"]["status"], "start codon lost")
        self.assertNotEqual(es[0]["description"]["status"], "stop codon lost")


class TestPartialIntronRetentionCodonInvariant(unittest.TestCase):
    """Pin that partial_intron_retention never reports codon loss (the retained
    bases are intronic; the native splice boundary is still removed, so neither
    cds_start nor cds_end can sit inside the retained segment by construction)."""

    def test_partial_intron_retention_codon_lost_never_true(self):
        """Reuse the BRCA1 c.5467+1G>A fixture which DOES produce a
        partial_intron_retention on the 3' side of the exon (verified by
        test_brca1_c5467_1g_a)."""
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.40, "DS_DG": 0.36, "DS_DL": 0.93,
            "DS_AG_ALT": 0, "DS_AL_ALT": 0,
            "DS_DG_ALT": 0.95, "DS_DL_ALT": 0.05,
            "DP_AG": -4, "DP_AL": 61, "DP_DG": -4, "DP_DL": 1,
            "EXON_STARTS": BRCA1_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA1_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": BRCA1_TRANSCRIPT_HG19["STRAND"],
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=41199659)
        pir = [a for a in result["aberrations"] if a["aberration_type"] == "partial_intron_retention"]
        # Guard against a vacuous test: the fixture MUST produce at least one
        # partial_intron_retention for the invariant assertion below to be
        # meaningful.
        self.assertGreaterEqual(len(pir), 1,
                                "fixture must produce a partial_intron_retention to exercise the invariant")
        for ab in pir:
            self.assertIn(ab["start_codon_lost"], (False, None),
                          "partial_intron_retention must never set start_codon_lost=True")
            self.assertIn(ab["stop_codon_lost"], (False, None),
                          "partial_intron_retention must never set stop_codon_lost=True")


class TestCodonLostFieldContract(unittest.TestCase):
    """Pinning tests for the output-dict field contract.

    Every aberration dict must carry cds_size, start_codon_lost, stop_codon_lost
    (possibly None). Types that structurally cannot lose a codon
    (whole_intron_retention, pseudoexon, increased_exon_inclusion,
    partial_intron_retention) must report start_codon_lost == False and
    stop_codon_lost == False when the transcript is coding and the aberration
    overlaps the CDS."""

    def test_whole_intron_retention_never_loses_codons(self):
        """Whole-intron retention adds bases; it must never report codon loss."""
        # BRCA2 intron 2 (between exon 2 end 32890664 and exon 3 start 32893214).
        # Variant at exon 2 donor (32890664) with DS_DL>0.2 AND DS_AL>0.2 and
        # orientation DP_AL > DP_DL -> whole_intron_retention.
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.50, "DS_DG": 0.00, "DS_DL": 0.80,
            "DP_AG": 0, "DP_AL": 2550, "DP_DG": 0, "DP_DL": 0,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32890664)
        wir = [a for a in result["aberrations"] if a["aberration_type"] == "whole_intron_retention"]
        self.assertEqual(len(wir), 1)
        self.assertIn("start_codon_lost", wir[0])
        self.assertIn("stop_codon_lost", wir[0])
        self.assertIn("cds_size", wir[0])
        # BRCA2 intron 2 lies within the CDS (cds_size > 0), so the coding
        # branch runs and set_coding_fields writes explicit False for both
        # codon-lost flags. Tighter than the None-allowed variant so we catch
        # accidental regressions to `None` on coding-overlap WIR.
        self.assertFalse(wir[0]["start_codon_lost"])
        self.assertFalse(wir[0]["stop_codon_lost"])
        self.assertGreater(wir[0]["cds_size"], 0, "fixture must have CDS overlap to exercise the coding branch")

    def test_all_aberrations_carry_uniform_fields(self):
        """Sanity: regardless of type, every emitted aberration dict has the new fields."""
        # Re-use the BRCA2 exon 2 c.-40+1G>A scenario (exon skipping + partial
        # exon deletion combination).
        transcript_scores = {
            "DS_AG": 0.30, "DS_AL": 0.77, "DS_DG": 0.00, "DS_DL": 0.17,
            "DP_AG": 20, "DP_AL": -41, "DP_DG": 0, "DP_DL": 64,
            "DS_AG_ALT": 0.5, "DS_AL_ALT": 0.3, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32890600)
        for ab in result["aberrations"]:
            self.assertIn("cds_size", ab, f"missing cds_size on {ab['aberration_type']}")
            self.assertIn("start_codon_lost", ab, f"missing start_codon_lost on {ab['aberration_type']}")
            self.assertIn("stop_codon_lost", ab, f"missing stop_codon_lost on {ab['aberration_type']}")


class TestExonSkippingFallbackReverseStrand(unittest.TestCase):
    """Reverse-strand coverage for the fallback path in sai10k_annotate_frameshift.

    On '-' strand, the affected exon's biological donor is at exon_starts[idx]
    (low genomic coord) and biological acceptor is at exon_ends[idx]. A
    CHEK2-style variant at the donor with an interior DP_AL peak exercises
    the is_reverse=True branch (sai10k_predictions.py:678)."""

    def test_fallback_on_minus_strand_donor_variant(self):
        """BRCA1 is on '-' strand. Pick an internal exon and place a variant
        exactly at its biological donor (exon_starts on '-' strand)."""
        exon_starts = BRCA1_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA1_TRANSCRIPT_HG19["EXON_ENDS"]
        # Genomic index 4 (biological exon 23 - 4 = 19; BRCA1 has 23 exons):
        # exon_starts[4]=41209069, exon_ends[4]=41209152, size=84 bp (84 % 3 = 0 -> in-frame).
        # On '-' strand, biological donor = exon_starts[4] = 41209069.
        variant_pos = 41209069
        bio_exon = len(exon_starts) - 4  # = 19
        exon_size = exon_ends[4] - exon_starts[4] + 1  # = 84
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.30, "DS_DG": 0.00, "DS_DL": 0.80,
            # DP_DL = 0 -> GEO_DL = variant_pos = exon_starts[4] (native donor on '-' strand)
            # DP_AL = 50 -> GEO_AL = 41209119, 50 bp inside exon (not at exon_ends[4] acceptor)
            "DP_AG": 0, "DP_AL": 50, "DP_DG": 0, "DP_DL": 0,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": exon_starts,
            "EXON_ENDS": exon_ends,
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "-",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        es = [a for a in result["aberrations"] if a["aberration_type"] == "exon_skipping"]
        self.assertEqual(len(es), 1)
        # Orientation on '-' strand for exon_skipping: dp_al*-1 < dp_dl*-1
        # -50 < 0, so exon_skipping (not whole_intron_retention). ✓
        self.assertEqual(es[0]["size"], exon_size)
        self.assertEqual(es[0]["description"], {
            "label": f"Exon {bio_exon} skipping",
            "size_bp": exon_size,
            "size_is_coding": True,
            "consequence": None,
            "status": "in-frame",
            "introduces_stop_codon": False,
            "extends_past_native_stop": False,
        })


class TestNativeSplicesitesTiebreak(unittest.TestCase):
    """Regression: the off-by-one fix in _get_native_splice_sites changed
    tiebreak behavior at exact equidistance from two flanking exons.
    Pin the '<=' upstream-favoring convention."""

    def test_variant_equidistant_picks_upstream(self):
        """BRCA2 intron 1 lies between exon_ends[0]=32889804 and
        exon_starts[1]=32890559 (754 bp gap). Midpoint: variant equidistant
        from both flanking exons. With the '<=' tiebreak, assoc_idx = 0
        (upstream)."""
        # Pick a gap where both distances fit inside 50 bp of the boundary
        # so the loss branch fires. Use an artificially small intron for the
        # test so the variant is within d_50bp of both flanking exons.
        exon_starts = [100, 221, 400]  # intron 1 spans (140, 220) = 80 bp
        exon_ends =   [140, 300, 500]
        # Midpoint of intron 1: 180 (equidistant from exon_ends[0]=140 and exon_starts[1]=221:
        # dist_to_upstream = 180 - 140 = 40; dist_to_downstream = 221 - 180 = 41.
        # Wait — these aren't quite equal. Use 180 and verify:
        # Actually for exactly-equal: variant_pos - exon_ends[0] == exon_starts[1] - variant_pos
        #   variant_pos = (exon_ends[0] + exon_starts[1]) / 2 = (140 + 221) / 2 = 180.5.
        # Integer midpoint isn't achievable with these coords. Use exon_starts[1]=221.
        # Adjust so midpoint is an integer: exon_ends[0]=140, exon_starts[1]=220 → midpoint 180.
        exon_ends = [140, 300, 500]
        exon_starts = [100, 220, 400]
        # variant_pos=180: dist_to_upstream = 180-140=40, dist_to_downstream = 220-180=40. Tie.

        # Build a scenario where the loss branch fires and assoc_idx determines
        # which exon's splice sites become "native" for the cryptic subflow.
        # We don't need to assert on cryptic output; we can directly invoke
        # the _get_native_splice_sites helper.
        from sai10k_predictions import _get_native_splice_sites, sai10k_find_affected_region
        region = sai10k_find_affected_region(180, exon_starts, exon_ends, strand="+")
        native = _get_native_splice_sites(region, 180, exon_starts, exon_ends, strand="+")
        self.assertEqual(native["assoc_idx"], 0)  # upstream exon due to `<=` tiebreak
        # Native donor of upstream exon on + strand = exon_ends[0] = 140
        self.assertEqual(native["geo_nd"], 140)


class TestPotentialWholeIntronRetention(unittest.TestCase):
    """Complements TestPotentialPrefix by exercising the 'Potential whole
    intron retention' label (sai10k_annotate_frameshift:726)."""

    def test_potential_whole_intron_retention_on_weak_paired_loss(self):
        """Variant in intron 1 with asymmetric loss pair and orientation
        favoring whole_intron_retention over exon_skipping."""
        exon_starts = BRCA2_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA2_TRANSCRIPT_HG19["EXON_ENDS"]
        # Intron 1 spans [32889805, 32890558] (754 bp). Place variant within
        # d_50bp of exon 2's acceptor. DP_AL matches exon 2 acceptor
        # (32890559), DP_DL matches exon 1 donor (32889804) — orientation
        # gives whole_intron_retention (dp_al*1 > dp_dl*1).
        variant_pos = 32890549  # 10 bp upstream of exon 2 acceptor
        transcript_scores = {
            "DS_AG": 0.0, "DS_AL": 0.93, "DS_DG": 0.0, "DS_DL": 0.13,  # weak donor partner
            "DP_AG": 0, "DP_AL": 10, "DP_DG": 0, "DP_DL": -745,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": exon_starts,
            "EXON_ENDS": exon_ends,
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        wir = [a for a in result["aberrations"] if a["aberration_type"] == "whole_intron_retention"]
        self.assertEqual(len(wir), 1)
        # Variant lies in intron 1 of BRCA2 (between exon 1 and exon 2 on +
        # strand), so the label embeds the biological intron number.
        self.assertTrue(wir[0]["description"]["label"].startswith("Potential intron 1 retention"))

    def test_suppressed_when_lost_sites_dont_match_native(self):
        """BRCA2 c.-40+5G>T at SpliceAI max_distance=500: the 500bp window
        doesn't reach exon 2's native acceptor (754bp downstream), so DP_AL
        peaks at a cryptic position 269bp downstream with weak DS_AL=0.03.
        Orientation still routes to whole_intron_retention, but geo_dl and
        geo_al do not jointly map to the native donor and native acceptor of
        a single intron — the prediction must be suppressed rather than
        emitted with a misleading 273bp size."""
        exon_starts = BRCA2_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA2_TRANSCRIPT_HG19["EXON_ENDS"]
        # c.-40+5 lies in intron 1, 5 bp downstream of exon 1's donor
        # (32889804). With max_distance=500, DP_AL peaks at a cryptic
        # acceptor 269 bp downstream of the variant, NOT at the native
        # acceptor of exon 2 (which is 754 bp downstream).
        variant_pos = 32889809
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.03, "DS_DG": 0.61, "DS_DL": 0.89,
            "DP_AG": 0, "DP_AL": 269, "DP_DG": -104, "DP_DL": -5,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": exon_starts,
            "EXON_ENDS": exon_ends,
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        types = [a["aberration_type"] for a in result["aberrations"]]
        self.assertNotIn("whole_intron_retention", types)

    def test_minus_strand_whole_intron_retention_native_match(self):
        """'-' strand coverage: BRCA1 transcript intron 17 (genomic intron
        between exon_ends[5]=41215390 and exon_starts[6]=41215891, 500 bp).
        Variant 1 bp into the intron from the donor; DP_DL/DP_AL point to the
        native donor (exon_starts[6]) and native acceptor (exon_ends[5]) —
        helper must accept and emit whole_intron_retention."""
        exon_starts = BRCA1_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA1_TRANSCRIPT_HG19["EXON_ENDS"]
        # On '-' strand the transcript donor of intron 17 sits at the LOW
        # genomic coord of the upstream-in-transcript exon (exon_starts[6]),
        # and the acceptor at the HIGH coord of the downstream exon
        # (exon_ends[5]). variant_pos = 41215890 is donor+1 in transcript order.
        variant_pos = 41215890
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.78, "DS_DG": 0.00, "DS_DL": 0.89,
            # DP_DL=1 -> geo_dl=41215891=exon_starts[6] (native donor)
            # DP_AL=-500 -> geo_al=41215390=exon_ends[5] (native acceptor)
            "DP_AG": 0, "DP_AL": -500, "DP_DG": 0, "DP_DL": 1,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": exon_starts,
            "EXON_ENDS": exon_ends,
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "-",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        wir = [a for a in result["aberrations"] if a["aberration_type"] == "whole_intron_retention"]
        self.assertEqual(len(wir), 1)
        self.assertEqual(wir[0]["size"], 500)

    def test_minus_strand_whole_intron_retention_suppressed(self):
        """'-' strand suppression: same BRCA1 fixture with DP_DL still pointing
        to the native donor but DP_AL pointing into the intron interior (not a
        native acceptor). Orientation still routes to WIR but the helper must
        reject — prediction suppressed."""
        exon_starts = BRCA1_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA1_TRANSCRIPT_HG19["EXON_ENDS"]
        variant_pos = 41215890
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.03, "DS_DG": 0.00, "DS_DL": 0.89,
            # DP_DL=1 -> native donor; DP_AL=-300 -> 41215590, intron interior.
            "DP_AG": 0, "DP_AL": -300, "DP_DG": 0, "DP_DL": 1,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": exon_starts,
            "EXON_ENDS": exon_ends,
            "CDS_START": BRCA1_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA1_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "-",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=variant_pos)
        types = [a["aberration_type"] for a in result["aberrations"]]
        self.assertNotIn("whole_intron_retention", types)


class TestDonorShiftPartialDescription(unittest.TestCase):
    """The donor-shift branch (_branch='gain_B_donor') uses seg_a = geo_nd
    (native donor) and shift_type='Donor'. Previously only the acceptor-shift
    branch had a description test."""

    def test_donor_shift_partial_exon_deletion(self):
        """Cryptic donor inside the exon upstream of the native donor should
        render as 'Donor shift (Nbp partial exon deletion) - ...'."""
        # BRCA2 exon 2: [32890559, 32890664], native donor at 32890664 on '+' strand.
        # Variant at 32890650 (14 bp upstream of native donor, still inside the exon).
        # DP_DG = -10 -> GEO_DG = 32890640 (cryptic donor 24 bp into the exon).
        # partial_size = (dp_nd - dp_dg) * strand = (14 - (-10)) * 1 = 24, positive -> partial exon deletion.
        # Need cryptic_donor_check: DS_DG >= 0.2 AND donor_diff > -0.2.
        transcript_scores = {
            "DS_AG": 0.0, "DS_AL": 0.0, "DS_DG": 0.35, "DS_DL": 0.0,
            "DP_AG": 0, "DP_AL": 0, "DP_DG": -10, "DP_DL": 0,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0,
            "DS_DG_ALT": 0.5, "DS_DL_ALT": 0.3,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32890650)
        pd = [a for a in result["aberrations"] if a["aberration_type"] == "partial_exon_deletion"]
        self.assertEqual(len(pd), 1)
        self.assertEqual(pd[0]["size"], 24)
        # 24 % 3 == 0 -> in-frame
        self.assertEqual(pd[0]["description"], {
            "label": "Donor shift",
            "size_bp": 24,
            "size_is_coding": True,
            "consequence": "partial exon 2 deletion",
            "status": "in-frame",
            "introduces_stop_codon": False,
            "extends_past_native_stop": False,
        })


class TestDeltaTypeTieBreak(unittest.TestCase):
    """When two or more Δ scores tie at the max, Python's max() returns the
    first item in the iterable. The order in sai10k_compute_predictions is:
    acceptor loss, donor loss, acceptor gain, donor gain. Pin this convention."""

    def test_al_dl_tie_picks_acceptor_loss(self):
        """DS_AL == DS_DL == 0.50 → delta_type = 'acceptor loss' (first)."""
        transcript_scores = {
            "DS_AG": 0.0, "DS_AL": 0.50, "DS_DG": 0.0, "DS_DL": 0.50,
            "DP_AG": 0, "DP_AL": -41, "DP_DG": 0, "DP_DL": 64,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32890600)
        self.assertTrue(result["aberrations"])
        self.assertEqual(result["aberrations"][0]["delta_type"], "acceptor loss")


# ===========================================================================
# Premature-stop detection tests (PDF feedback comment #5).
# ===========================================================================

# Path to a local hg38 FASTA. Used by integration tests below; tests that need
# this file are skipped when it isn't present so the suite still runs in
# environments without a reference genome.
HG38_FASTA_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "hg38.fa"
)

# RAD51C transcript ENST00000583539.5 in hg38 (extracted from
# gencode.v49.GRCh38.txt). + strand, 8 exons. The motivating example for PDF
# feedback comment #5 — variant chr17:58695191 T>C is c.404+2T>C, expected to
# produce stop_codon_introduced=True with aa_change="TQL[WQNKVFSF*]".
RAD51C_TRANSCRIPT_HG38 = {
    "transcript_id": "ENST00000583539.5",
    "gene_name": "RAD51C",
    "STRAND": "+",
    "TX_START": 58692573,
    "TX_END": 58732766,
    "CDS_START": 58692644,
    "CDS_END": 58732559,
    "EXON_STARTS": [58692573, 58694931, 58696693, 58703196, 58709859,
                    58720746, 58724040, 58732484],
    "EXON_ENDS":   [58692788, 58695189, 58696859, 58703329, 58709990,
                    58720812, 58724100, 58732766],
}


class TestTranslateDna(unittest.TestCase):
    def test_basic_translation(self):
        self.assertEqual(_translate_dna("ATGAAATAA"), "MK*")

    def test_lower_case_input(self):
        self.assertEqual(_translate_dna("atgaaataa"), "MK*")

    def test_partial_codon_at_end(self):
        # Trailing non-multiple-of-3 bases are dropped (no padding).
        self.assertEqual(_translate_dna("ATGAA"), "M")

    def test_unknown_codon_yields_X(self):
        self.assertEqual(_translate_dna("ATGNNN"), "MX")


class TestReverseComplement(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(_reverse_complement("ATCG"), "CGAT")

    def test_n_passthrough(self):
        self.assertEqual(_reverse_complement("ATNCG"), "CGNAT")


class TestFindDifferencePoint(unittest.TestCase):
    def test_identical(self):
        self.assertEqual(_find_difference_point("ABC", "ABC", "forward"), "identical")

    def test_first_position_differs(self):
        self.assertEqual(_find_difference_point("XBC", "ABC", "forward"), 1)

    def test_third_position_differs(self):
        self.assertEqual(_find_difference_point("ABXY", "ABCD", "forward"), 3)

    def test_prefix(self):
        # Shorter sequence is a prefix of the longer one: returns len(shorter)+1.
        self.assertEqual(_find_difference_point("ABC", "ABCDE", "forward"), 4)

    def test_reverse_walk(self):
        # Walk from the end. seq_a "ABCDX" vs seq_b "ABCDY": last char differs.
        # Reversed: "XDCBA" vs "YDCBA", differs at position 1 in reversed coords,
        # which back-translates to position len-1+1 = len = 5 in forward coords.
        self.assertEqual(_find_difference_point("ABCDX", "ABCDY", "reverse"), 5)


class TestFormatAaChangeAndDetect(unittest.TestCase):
    """Unit tests for the AA-comparison / divergence-window logic.

    These hand-craft altered/consensus AA strings to exercise the algorithm
    independent of FASTA / transcript inputs.
    """

    def test_identical_returns_no_change(self):
        # Both proteins end in '*' (the natural terminal stop); after the
        # trailing-* strip the strings are identical.
        result = _format_aa_change_and_detect("ABCDEFG*", "ABCDEFG*")
        self.assertEqual(result, (False, None, False))

    def test_premature_stop_introduced(self):
        # Forward divergence at AA 4. Altered hits a premature stop a few AAs
        # in. Both end with the natural '*' which is stripped first.
        altered = "ABC*XYZ*"
        consensus = "ABCDEFGHIJK*"
        stop_introduced, aa_change, extends_past = _format_aa_change_and_detect(altered, consensus)
        self.assertTrue(stop_introduced)
        self.assertFalse(extends_past)
        # The aa_change dict should carry the prefix "ABC" and a changed_aa
        # segment ending in "*".
        self.assertEqual(aa_change["prefix_aa"], "ABC")
        self.assertTrue(aa_change["changed_aa"].endswith("*"))

    def test_native_stop_shifted_earlier_not_misclassified(self):
        # Critical regression: an in-frame deletion that shifts the native stop
        # earlier in transcription order leaves a trailing '*' in the altered
        # protein. Without the trailing-* strip, .find('*') would misclassify
        # the natural stop as introduced.
        altered = "ABCDEFG*"           # native stop just shifted earlier in position
        consensus = "ABCDEFGHIJ*"      # consensus with the natural stop at the end
        stop_introduced, _, _ = _format_aa_change_and_detect(altered, consensus)
        self.assertFalse(stop_introduced)

    def test_frameshift_introduces_stop(self):
        # Frameshift case: altered translation diverges at AA 4 and hits a PTC
        # a few residues later. The consensus carries its natural terminal '*'
        # (always stripped); the altered '*' must be detected as introduced.
        altered = "ABCXY*"               # frameshift garbage, then PTC
        consensus = "ABCDEFGHIJK*"
        stop_introduced, aa_change, extends_past = _format_aa_change_and_detect(
            altered, consensus, frameshift=True,
        )
        self.assertTrue(stop_introduced)
        self.assertFalse(extends_past)
        self.assertEqual(aa_change["prefix_aa"], "ABC")
        self.assertTrue(aa_change["changed_aa"].endswith("*"))

    def test_frameshift_extends_past_native_stop(self):
        # Frameshift case: altered translation diverges and never hits a PTC
        # within the translated CDS-clamped spans. Mirrors reference's
        # "protein sequence extends beyond native stop site".
        altered = "ABCXYZUVW"            # no '*' anywhere
        consensus = "ABCDEFGHIJK*"
        stop_introduced, aa_change, extends_past = _format_aa_change_and_detect(
            altered, consensus, frameshift=True,
        )
        self.assertFalse(stop_introduced)
        self.assertIsNone(aa_change)
        self.assertTrue(extends_past)

    def test_frameshift_natural_stop_in_altered_is_real_ptc(self):
        # Sanity check: for frameshift, the trailing-`*` strip must NOT apply
        # to altered_aa. If altered_aa happens to end in '*' (because the PTC
        # lands at the very end of the translated spans), it must be detected
        # as introduced. Without the strip-disable, the PTC would be silenced.
        altered = "ABCXYZ*"
        consensus = "ABCDEFGHIJK*"
        stop_introduced, _, _ = _format_aa_change_and_detect(
            altered, consensus, frameshift=True,
        )
        self.assertTrue(stop_introduced)

    def test_in_frame_default_unchanged(self):
        # Defaulting frameshift=False keeps existing in-frame behavior intact:
        # both natural stops are stripped, divergence-window logic runs.
        altered = "ABCDEFG*"
        consensus = "ABCDEFGHIJ*"
        stop_introduced, _, _ = _format_aa_change_and_detect(altered, consensus)
        self.assertFalse(stop_introduced)


class TestApplyVariantToAlteredExonSeqs(unittest.TestCase):
    """The variant-application helper supports indels and indel-vs-splice
    boundary detection. These tests use a stub fasta object — no real genome
    file required."""

    class _StubFasta:
        def __init__(self, sequences):
            # sequences: dict of chrom -> sequence string (1-based, supports slicing).
            self._seqs = sequences

        def keys(self):
            return list(self._seqs.keys())

        def __getitem__(self, chrom):
            return self._StubChrom(self._seqs[chrom])

        class _StubChrom:
            def __init__(self, seq):
                self._seq = seq

            def __getitem__(self, slc):
                return self._seq[slc]

    def _stub(self, **kwargs):
        return self._StubFasta(kwargs)

    def test_snv_substitution(self):
        # Spans = [(11, 20)], exon_seqs = ["ACGTACGTAC"], variant at pos 13 T>G.
        spans = [(11, 20)]
        exon_seqs = ["ACGTACGTAC"]  # bases at genomic 11..20
        # Stub fasta returns the same sequence so reference matches.
        # FASTA contains a region [10..20] -> "_ACGTACGTAC" (1-based: pos 11 = "A").
        # We supply 0-based slice indexing: pos[10:11] = "A".
        # Build a "chr1" sequence such that 1-based pos 13 = "G". Our exon string starts at pos 11.
        # exon[0]=A(pos11), exon[1]=C(pos12), exon[2]=G(pos13), exon[3]=T(pos14),...
        # ref at pos 13 should be "G".
        # The stub exposes the chromosome as a string indexed 0-based, so we
        # need bases [0..20] with pos11=A, pos12=C, pos13=G, etc.
        chrom_seq = "X" * 10 + "ACGTACGTAC"  # pos 11..20
        fasta = self._stub(chr1=chrom_seq)
        # 1-based variant: pos 13 ref=G alt=T. cds_start/cds_end far away.
        result = _apply_variant_to_altered_exon_seqs(
            spans, list(exon_seqs), 13, "G", "T",
            fasta, "1", cds_start=1, cds_end=10000,
        )
        self.assertEqual(result, "0_AC" + "T" + "TACGTAC")  # pos 13 G -> T at offset 2

    def test_reference_mismatch(self):
        chrom_seq = "X" * 10 + "ACGTACGTAC"
        fasta = self._stub(chr1=chrom_seq)
        # Claim pos 13 = "A", but actual genome has "G".
        result = _apply_variant_to_altered_exon_seqs(
            [(11, 20)], ["ACGTACGTAC"], 13, "A", "T",
            fasta, "1", cds_start=1, cds_end=10000,
        )
        self.assertEqual(result, "reference mismatch")

    def test_variant_outside_kept_spans(self):
        # variant footprint at pos 30 (outside [(11, 20)]).
        chrom_seq = "X" * 30 + "G"
        fasta = self._stub(chr1=chrom_seq)
        result = _apply_variant_to_altered_exon_seqs(
            [(11, 20)], ["ACGTACGTAC"], 31, "G", "T",
            fasta, "1", cds_start=1, cds_end=10000,
        )
        self.assertEqual(result, "cannot determine")

    def test_indel_straddles_splice_boundary(self):
        # 3-base ref starting at pos 19, footprint = [19, 20, 21]. Span ends at 20.
        chrom_seq = "X" * 10 + "ACGTACGTAC" + "TGA"  # pos 21 = "T"
        fasta = self._stub(chr1=chrom_seq)
        result = _apply_variant_to_altered_exon_seqs(
            [(11, 20)], ["ACGTACGTAC"], 19, "ACT", "G",
            fasta, "1", cds_start=1, cds_end=10000,
        )
        self.assertEqual(result, "variant straddles splice boundary")

    def test_anchored_insertion_at_trailing_edge_treated_as_straddle(self):
        # VCF-anchored insertion at the last base of a kept span: ref='C', alt='CTAA' at
        # pos=20 in span (11, 20). The anchor 'C' is the last exonic base; by VCF
        # convention the inserted "TAA" sits past the donor in the intron and must NOT
        # be spliced into the exon sequence (which would falsely look like a 3-AA
        # insertion that could carry a premature stop). Expected: straddle sentinel.
        chrom_seq = "X" * 10 + "ACGTACGTAC"  # pos 20 = "C"
        fasta = self._stub(chr1=chrom_seq)
        result = _apply_variant_to_altered_exon_seqs(
            [(11, 20)], ["ACGTACGTAC"], 20, "C", "CTAA",
            fasta, "1", cds_start=1, cds_end=10000,
        )
        self.assertEqual(result, "variant straddles splice boundary")

    def test_indel_substitution_inside_exon(self):
        # 3-base deletion (in-frame): ref="GTA" alt="" at pos 13 (positions 13,14,15).
        chrom_seq = "X" * 10 + "ACGTACGTAC"  # pos 13="G", 14="T", 15="A"
        fasta = self._stub(chr1=chrom_seq)
        result = _apply_variant_to_altered_exon_seqs(
            [(11, 20)], ["ACGTACGTAC"], 13, "GTA", "",
            fasta, "1", cds_start=1, cds_end=10000,
        )
        # pos 13 is at offset 2 in the exon. exon = "AC|GTA|CGTAC" — after
        # deleting "GTA": exon[:2] + "" + exon[5:] = "AC" + "CGTAC" = "ACCGTAC".
        self.assertEqual(result, "0_ACCGTAC")

    def test_variant_at_start_codon(self):
        # cds_start=13. A SNV at pos 14 (cds_start+1) lands on the start codon.
        chrom_seq = "X" * 10 + "ACGTACGTAC"
        fasta = self._stub(chr1=chrom_seq)
        result = _apply_variant_to_altered_exon_seqs(
            [(11, 20)], ["ACGTACGTAC"], 14, "T", "C",
            fasta, "1", cds_start=13, cds_end=10000,
        )
        self.assertEqual(result, "impacts native start or stop site")

    def test_variant_at_high_cds_codon(self):
        # Locks in the high-window (cds_end - 2 .. cds_end) branch. On a
        # '-' strand transcript this 3 bp window holds the biological ATG;
        # on '+' strand it holds the stop codon. Either way the function
        # must return the same sentinel — proves detection is symmetric
        # and not strand-dependent.
        chrom_seq = "X" * 10 + "ACGTACGTAC"  # pos 17 = "G"
        fasta = self._stub(chr1=chrom_seq)
        result = _apply_variant_to_altered_exon_seqs(
            [(11, 20)], ["ACGTACGTAC"], 17, "G", "A",
            fasta, "1", cds_start=1, cds_end=18,
        )
        self.assertEqual(result, "impacts native start or stop site")


class TestComputeConsensusContextSliceIdentity(unittest.TestCase):
    """Lock down the invariant the Tier 1 hoist depends on:

        consensus_aa_full[(cum_cds_bases[i] + 2) // 3:]
            == _translate_spans(consensus[i:], cum_cds_bases[i] % 3, ...)

    for every row_anchor i and for all three first_eframe residues. If anyone
    changes the eFrame leading-base trim in _translate_spans without updating
    the slice formula in _compute_consensus_context, this test catches it.
    """

    class _StubChrom:
        def __init__(self, seq): self._seq = seq
        def __getitem__(self, slc): return self._seq[slc]

    class _StubFasta:
        def __init__(self, seq): self._seq = seq
        def keys(self): return ["chr1"]
        def __getitem__(self, _): return TestComputeConsensusContextSliceIdentity._StubChrom(self._seq)

    def _check(self, exon_starts, exon_ends, cds_start, cds_end, strand, seq):
        scores = {
            "EXON_STARTS": exon_starts,
            "EXON_ENDS": exon_ends,
            "CDS_START": cds_start,
            "CDS_END": cds_end,
            "STRAND": strand,
        }
        fasta = self._StubFasta(seq)
        ctx = _compute_consensus_context(scores, fasta, "1")
        self.assertIsNotNone(ctx)
        consensus = ctx["consensus"]
        cum = ctx["cum_cds_bases"]
        full = ctx["consensus_aa_full"]
        for i in range(len(consensus)):
            spans = [(ex["cds_lo"], ex["cds_hi"]) for ex in consensus[i:]]
            expected = _translate_spans(spans, cum[i] % 3, fasta, "1", strand)
            sliced = full[(cum[i] + 2) // 3:]
            self.assertEqual(sliced, expected,
                f"slice mismatch at row_anchor={i}, eframe={cum[i] % 3}")

    def test_plus_strand_all_eframe_residues(self):
        # 4 CDS-only exons with sizes 10, 11, 12, 13 -> cumulative 0,10,21,33,46.
        # Residues mod 3 at the boundaries: 0, 1, 0, 0, 1 -- covers eframe 0 and 1.
        # Pad with one more exon of size 14 to reach a residue-2 anchor.
        # cum -> 0, 10, 21, 33, 46, 60 -> mod 3: 0, 1, 0, 0, 1, 0. Need a size that
        # makes a residue-2: e.g. sizes 10, 11, 13 give cum 0,10,21,34 -> mod 3: 0,1,0,1.
        # Use sizes 11, 11, 11 -> cum 0,11,22,33 -> mod 3: 0,2,1,0. All three residues hit.
        # Place exons at 101..111, 201..211, 301..311. CDS = whole.
        seq = "X" * 100 + "ACGTACGTACG" + "X" * 89 + "TTTAAACCCGG" + "X" * 89 + "GGGCCCAAATT"
        self._check([101, 201, 301], [111, 211, 311], 101, 311, "+", seq)

    def test_minus_strand_all_eframe_residues(self):
        # Same exon layout on - strand; tx order reverses, so row_anchor 0 = highest-coord exon.
        seq = "X" * 100 + "ACGTACGTACG" + "X" * 89 + "TTTAAACCCGG" + "X" * 89 + "GGGCCCAAATT"
        self._check([101, 201, 301], [111, 211, 311], 101, 311, "-", seq)

    def test_cds_clamped_first_and_last_exon(self):
        # First and last exons partly UTR; verify clamping doesn't break the identity.
        # Exon sizes (after clamp): 5, 11, 8 -> cum 0, 5, 16, 24 -> mod 3: 0, 2, 1, 0.
        seq = "X" * 100 + "ACGTACGTACG" + "X" * 89 + "TTTAAACCCGG" + "X" * 89 + "GGGCCCAAATT"
        # CDS_START=107 (clip exon1's first 6 bases), CDS_END=308 (clip exon3's last 3 bases).
        self._check([101, 201, 301], [111, 211, 311], 107, 308, "+", seq)

    def test_context_built_when_tiny_interior_exon_does_not_become_anchor(self):
        # CDS exon sizes [4, 1, 30] → cum [0, 4, 5, 35] → mod 3 [0, 1, 2, 2].
        # An interior 1-bp exon at cum%3==1 is fine for the *consensus* (the
        # full translation is computed at eFrame=0, no trim). It is only the
        # per-aberration altered translation that has trouble when the
        # row_anchor lands on that 1-bp exon — and that's guarded inside
        # _detect_premature_stop, not here.
        scores = {
            "EXON_STARTS": [101, 201, 301],
            "EXON_ENDS":   [104, 201, 330],
            "CDS_START": 101,
            "CDS_END":   330,
            "STRAND": "+",
        }
        seq = "X" * 100 + "ACGT" + "X" * 96 + "A" + "X" * 99 + ("ACGT" * 8)[:30]
        ctx = _compute_consensus_context(scores, self._StubFasta(seq), "1")
        self.assertIsNotNone(ctx)
        self.assertEqual(ctx["cum_cds_bases"], [0, 4, 5, 35])

    def test_whole_intron_retention_with_tiny_downstream_exon_detects_stop(self):
        # Regression for the over-broad 1-bp guard. whole_intron_retention
        # anchors at flanking_idx=0 (the consensus exon BEFORE the retained
        # intron), so first_eframe=0 and the trim hazard is not in play —
        # detection should produce M[*]. Constructed so the retained intron
        # contains TAA in-frame after ATG.
        seq = list("C" * 330)
        # exon 1 (101..104) = "ATGT"
        for pos, base in zip(range(101, 105), "ATGT"):
            seq[pos - 1] = base
        # intron contributes "AAC..." after the variant (pos 105 C→A makes
        # the codon TAA across the splice junction in the altered transcript).
        seq[104] = "C"   # genome pos 105 ref C
        seq[105] = "A"   # genome pos 106
        # tiny exon at pos 201, 1bp 'G'
        seq[200] = "G"
        # exon 3 (301..330) repeating "CAG" -> Q codons, no extra stop
        for i, pos in enumerate(range(301, 331)):
            seq[pos - 1] = "CAG"[i % 3]
        scores = {
            "DS_AG": 0.0, "DS_AL": 0.8, "DS_DG": 0.0, "DS_DL": 0.8,
            "DP_AG": 0, "DP_AL": 96, "DP_DG": 0, "DP_DL": -1,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0, "DS_DG_ALT": 0.0, "DS_DL_ALT": 0.0,
            "EXON_STARTS": [101, 201, 301],
            "EXON_ENDS":   [104, 201, 330],
            "CDS_START": 101, "CDS_END": 330, "STRAND": "+",
        }
        result = sai10k_compute_predictions(
            scores, variant_pos=105, chrom="1", ref="C", alt="A",
            fasta=self._StubFasta("".join(seq)),
        )
        wirs = [a for a in result["aberrations"]
                if a["aberration_type"] == "whole_intron_retention"]
        self.assertEqual(len(wirs), 1)
        ab = wirs[0]
        self.assertEqual(ab["stop_codon_introduced"], True,
                         f"aa_change={ab.get('aa_change')}")
        self.assertEqual(ab["aa_change"], {"prefix_aa": "M", "changed_aa": "*"})


class TestConsensusFilteredInTxOrder(unittest.TestCase):
    def test_plus_strand_drops_utr_only(self):
        # 3-exon transcript on + strand. CDS = [101, 200]. Exon 1 fully UTR (1-100),
        # exon 2 partially CDS (50-150), exon 3 fully CDS (180-200).
        scores = {
            "EXON_STARTS": [1, 50, 180],
            "EXON_ENDS":   [100, 150, 200],
            "CDS_START": 101,
            "CDS_END": 200,
            "STRAND": "+",
        }
        out = _consensus_filtered_in_tx_order(scores)
        self.assertEqual(len(out), 2)  # UTR-only exon 1 dropped
        self.assertEqual(out[0]["cds_lo"], 101)  # exon 2 clamped to cds_start
        self.assertEqual(out[0]["cds_hi"], 150)
        self.assertEqual(out[0]["eNum"], 2)
        self.assertEqual(out[1]["eNum"], 3)

    def test_minus_strand_tx_order_reversed(self):
        # Same transcript but on - strand; tx-order index 0 is the highest-coord exon.
        scores = {
            "EXON_STARTS": [1, 50, 180],
            "EXON_ENDS":   [100, 150, 200],
            "CDS_START": 1,
            "CDS_END": 200,
            "STRAND": "-",
        }
        out = _consensus_filtered_in_tx_order(scores)
        self.assertEqual(len(out), 3)
        self.assertEqual(out[0]["orig_start"], 180)  # genomically-highest = first in tx order
        self.assertEqual(out[0]["eNum"], 1)
        self.assertEqual(out[2]["orig_start"], 1)
        self.assertEqual(out[2]["eNum"], 3)


class TestPrematureStopBackwardCompat(unittest.TestCase):
    """When fasta is None (or chrom/ref/alt aren't passed), detection is
    silently skipped: stop_codon_introduced and aa_change are None on every
    aberration, and description.introduces_stop_codon stays False."""

    def test_no_fasta_yields_none(self):
        # BRCA2 c.67+1G>A — known to produce an exon_skipping aberration.
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.81, "DS_DG": 0.17, "DS_DL": 0.98,
            "DP_AG": -221, "DP_AL": -106, "DP_DG": -137, "DP_DL": -1,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": BRCA2_TRANSCRIPT_HG19["STRAND"],
        }
        result = sai10k_compute_predictions(transcript_scores, variant_pos=32890665)
        self.assertTrue(result["aberrations"])
        for ab in result["aberrations"]:
            self.assertIsNone(ab["stop_codon_introduced"])
            self.assertIsNone(ab["aa_change"])
            # The description PTC flag must NOT have been flipped on.
            self.assertFalse(ab["description"]["introduces_stop_codon"])

    def test_no_fasta_skips_frameshift_aberration_too(self):
        # When fasta=None, detection short-circuits regardless of aberration
        # frame status (frameshift aberrations are now handled when fasta is
        # available — see _format_aa_change_and_detect).
        # Construct a frameshift exon_skipping: BRCA2 exon 5 (32900238..32900287,
        # size 50, 50 % 3 == 2). Variant at the acceptor splice site.
        transcript_scores = {
            "DS_AG": 0.00, "DS_AL": 0.95, "DS_DG": 0.00, "DS_DL": 0.95,
            "DP_AG": 0, "DP_AL": -1, "DP_DG": 0, "DP_DL": 49,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }
        result = sai10k_compute_predictions(
            transcript_scores, variant_pos=32900238,
            chrom="13", ref="G", alt="A", fasta=None,
        )
        # At least one frameshift aberration must be present so this test
        # actually exercises the no-FASTA short-circuit on a frameshift case.
        self.assertTrue(any(ab.get("frameshift") is True for ab in result["aberrations"]),
                        "expected at least one frameshift aberration in this fixture")
        for ab in result["aberrations"]:
            self.assertIsNone(ab["stop_codon_introduced"])
            self.assertIsNone(ab["aa_change"])


class TestPrematureStopDescriptionMapping(unittest.TestCase):
    """The mapping from _detect_premature_stop's outputs onto the
    aberration['description'] booleans (sai10k_predictions.py around L1419-1427)
    has three live branches:
        stop_introduced=True                       -> introduces_stop_codon=True
        frameshift=True and extends_past=True      -> extends_past_native_stop=True
        otherwise                                  -> both stay False
    The integration test (TestPrematureStopRad51c) covers only the in-frame PTC
    branch via a real-genome variant. These mock-based tests exercise the
    frameshift-PTC and extends-past branches without needing exotic real
    variants. They also reuse the BRCA2 frameshift exon_skipping fixture from
    TestPrematureStopBackwardCompat to keep the aberration shape realistic."""

    def _frameshift_es_scores(self):
        # Same fixture as TestPrematureStopBackwardCompat.test_no_fasta_skips_
        # frameshift_aberration_too: BRCA2 exon 5 (size 50, 50 % 3 == 2) -> a
        # frameshift exon_skipping aberration. We need a frameshift aberration
        # going into the detection block so the mapping branches are reached.
        return {
            "DS_AG": 0.00, "DS_AL": 0.95, "DS_DG": 0.00, "DS_DL": 0.95,
            "DP_AG": 0, "DP_AL": -1, "DP_DG": 0, "DP_DL": 49,
            "EXON_STARTS": BRCA2_TRANSCRIPT_HG19["EXON_STARTS"],
            "EXON_ENDS": BRCA2_TRANSCRIPT_HG19["EXON_ENDS"],
            "CDS_START": BRCA2_TRANSCRIPT_HG19["CDS_START"],
            "CDS_END": BRCA2_TRANSCRIPT_HG19["CDS_END"],
            "STRAND": "+",
        }

    def test_frameshift_ptc_sets_introduces_stop_codon(self):
        # _detect_premature_stop returns stop_introduced=True for a frameshift
        # aberration -> description['introduces_stop_codon'] must flip True.
        # Pass a non-None fasta sentinel so the detection block isn't skipped.
        from unittest.mock import patch
        import sai10k_predictions as sai
        fake_aa_change = {"prefix_aa": "ABC", "changed_aa": "XY*"}
        fake_window = {
            "prefix_hidden_aa": 10, "visible_aa": "ABC", "changed_aa": None,
            "suffix_hidden_aa": 100, "total_aa": 200,
        }
        with patch.object(sai, "_detect_premature_stop",
                          return_value=(True, fake_aa_change, fake_window, fake_window, False)):
            result = sai10k_compute_predictions(
                self._frameshift_es_scores(), variant_pos=32900238,
                chrom="13", ref="G", alt="A", fasta=object(),
            )
        frameshift_aberrations = [
            ab for ab in result["aberrations"]
            if ab.get("frameshift") is True and ab.get("affects_coding") is True
        ]
        self.assertTrue(frameshift_aberrations,
                        "expected at least one coding frameshift aberration in fixture")
        for ab in frameshift_aberrations:
            self.assertEqual(ab["stop_codon_introduced"], True)
            self.assertEqual(ab["aa_change"], fake_aa_change)
            self.assertTrue(ab["description"]["introduces_stop_codon"])
            self.assertFalse(ab["description"]["extends_past_native_stop"])

    def test_extends_past_sets_extends_past_native_stop(self):
        # _detect_premature_stop returns stop_introduced=False, extends_past=True
        # for a frameshift aberration -> description['extends_past_native_stop']
        # must flip True (and introduces_stop_codon must stay False).
        from unittest.mock import patch
        import sai10k_predictions as sai
        with patch.object(sai, "_detect_premature_stop",
                          return_value=(False, None, None, None, True)):
            result = sai10k_compute_predictions(
                self._frameshift_es_scores(), variant_pos=32900238,
                chrom="13", ref="G", alt="A", fasta=object(),
            )
        frameshift_aberrations = [
            ab for ab in result["aberrations"]
            if ab.get("frameshift") is True and ab.get("affects_coding") is True
        ]
        self.assertTrue(frameshift_aberrations,
                        "expected at least one coding frameshift aberration in fixture")
        for ab in frameshift_aberrations:
            self.assertEqual(ab["stop_codon_introduced"], False)
            self.assertIsNone(ab["aa_change"])
            self.assertFalse(ab["description"]["introduces_stop_codon"])
            self.assertTrue(ab["description"]["extends_past_native_stop"])


@unittest.skipUnless(os.path.exists(HG38_FASTA_PATH),
                     f"hg38 FASTA not found at {HG38_FASTA_PATH}")
class TestPrematureStopRad51c(unittest.TestCase):
    """Integration test for the motivating PDF example: NM_058216.3(RAD51C):c.404+2T>C
    (chr17:58695191 T>C, hg38). + strand transcript. Expected:
        stop_codon_introduced = True
        aa_change             = "TQL[WQNKVFSF*]"
    Also serves as the X4 chr-prefix regression test (chrom arrives as "17"
    while the FASTA's contigs are "chr*")."""

    @classmethod
    def setUpClass(cls):
        import pyfastx
        cls.fasta = pyfastx.Fasta(HG38_FASTA_PATH)

    def _build_scores(self):
        # SpliceAI scores chosen to trigger the partial_intron_retention
        # cryptic-donor branch:
        #   - DS_DG = 0.9 (>= MIN_DELTA_SCORE; ds_ag < ds_dg, geo_dg != native donor).
        #   - DP_DG = 25 -> cryptic donor at chr17:58695216 (27 bp past native donor at 58695189).
        #   - DP_DL = -2 -> native donor at chr17:58695189; partial_size = -27.
        #   - DS_DG_ALT >= ALT_SCORE_STRONG; donor_diff > ALT_DIFF_LOW.
        return {
            "DS_AG": 0.0, "DS_AL": 0.0, "DS_DG": 0.9, "DS_DL": 0.0,
            "DP_AG": 0, "DP_AL": 0, "DP_DG": 25, "DP_DL": -2,
            "DS_AG_ALT": 0.0, "DS_AL_ALT": 0.0,
            "DS_DG_ALT": 0.95, "DS_DL_ALT": 0.05,
            "EXON_STARTS": RAD51C_TRANSCRIPT_HG38["EXON_STARTS"],
            "EXON_ENDS": RAD51C_TRANSCRIPT_HG38["EXON_ENDS"],
            "CDS_START": RAD51C_TRANSCRIPT_HG38["CDS_START"],
            "CDS_END": RAD51C_TRANSCRIPT_HG38["CDS_END"],
            "STRAND": "+",
        }

    def test_rad51c_c404_2t_c_introduces_stop(self):
        result = sai10k_compute_predictions(
            self._build_scores(), variant_pos=58695191,
            chrom="17", ref="T", alt="C", fasta=self.fasta,
        )
        partials = [a for a in result["aberrations"]
                    if a["aberration_type"] == "partial_intron_retention"]
        self.assertEqual(len(partials), 1, f"expected 1 partial_intron_retention, "
                                            f"got {[a['aberration_type'] for a in result['aberrations']]}")
        ab = partials[0]
        self.assertEqual(ab["size"], 27)
        self.assertTrue(ab["affects_coding"])
        self.assertFalse(ab["frameshift"])
        self.assertEqual(ab["stop_codon_introduced"], True,
                         f"aa_change={ab['aa_change']}, description={ab['description']}")
        self.assertEqual(ab["aa_change"], {"prefix_aa": "TQL", "changed_aa": "WQNKVFSF*"})
        # In-frame PTC: introduces_stop_codon=True with frameshift=False so the
        # frontend renders the " but introduces a stop codon" join.
        self.assertTrue(ab["description"]["introduces_stop_codon"])
        self.assertFalse(ab["frameshift"])
        # Windowed protein view (±15 aa flanking the change, with truncation
        # markers and total residue counts). Bracket aligns to the first
        # divergent residue; PTC terminates the predicted line. Pinned values:
        #   prefix_offset_aa = 119 (residues upstream of row_anchor)
        #   total WT = 346 AA, total predicted = 142 AA (terminates at PTC).
        self.assertEqual(
            ab["wt_protein_window"],
            {
                "prefix_hidden_aa": 119,
                "visible_aa": "TTEICGAPGVGKTQLCMQLAVDVQIPECFG",
                "changed_aa": None,
                "trailing_visible_aa": "",
                "suffix_hidden_aa": 197,
                "total_aa": 346,
            },
        )
        self.assertEqual(
            ab["altered_protein_window"],
            {
                "prefix_hidden_aa": 119,
                "visible_aa": "TTEICGAPGVGKTQL",
                "changed_aa": "WQNKVFSF*",
                "trailing_visible_aa": "",
                "suffix_hidden_aa": 0,
                "total_aa": 142,
            },
        )

    def test_reference_mismatch_skipped(self):
        # Same variant but claim ref='A' (genome has 'T'). Detection must be
        # skipped: stop_codon_introduced=None and description unchanged.
        result = sai10k_compute_predictions(
            self._build_scores(), variant_pos=58695191,
            chrom="17", ref="A", alt="C", fasta=self.fasta,
        )
        partials = [a for a in result["aberrations"]
                    if a["aberration_type"] == "partial_intron_retention"]
        self.assertEqual(len(partials), 1)
        ab = partials[0]
        self.assertIsNone(ab["stop_codon_introduced"])
        self.assertIsNone(ab["aa_change"])
        self.assertIsNone(ab["wt_protein_window"])
        self.assertIsNone(ab["altered_protein_window"])
        self.assertFalse(ab["description"]["introduces_stop_codon"])


class TestBuildProteinWindowsInFrame(unittest.TestCase):
    """Direct tests of _build_protein_windows in in-frame (no-PTC) mode.

    Use a small flank (3) so the window boundaries fit in short fixture
    sequences. prefix_offset_aa=10 simulates ten WT residues upstream of
    `consensus_aa`, exercising prefix_hidden_aa accounting.
    """

    def test_substitution_brackets_both_sides(self):
        # WT: A B C D E F G   altered: A B C X Y G
        # forward_diff at idx 3 (D vs X); matching tail "G".
        wt, alt = _build_protein_windows(
            altered_aa="ABCXYG", consensus_aa="ABCDEFG",
            frameshift=False, prefix_offset_aa=10, flank=3,
            stop_introduced=False,
        )
        self.assertEqual(wt, {
            "prefix_hidden_aa": 10,
            "visible_aa": "ABC",
            "changed_aa": "DEF",
            "trailing_visible_aa": "G",
            "suffix_hidden_aa": 0,
            "total_aa": 17,
        })
        self.assertEqual(alt, {
            "prefix_hidden_aa": 10,
            "visible_aa": "ABC",
            "changed_aa": "XY",
            "trailing_visible_aa": "G",
            "suffix_hidden_aa": 0,
            "total_aa": 16,
        })

    def test_pure_deletion_marks_altered_changed_empty(self):
        # WT: A B C D E F G   altered: A B C F G  (deletes "DE")
        wt, alt = _build_protein_windows(
            altered_aa="ABCFG", consensus_aa="ABCDEFG",
            frameshift=False, prefix_offset_aa=0, flank=3,
            stop_introduced=False,
        )
        # WT has the deleted residues bracketed.
        self.assertEqual(wt["changed_aa"], "DE")
        # Altered has '' (sentinel — frontend renders a `|` marker).
        self.assertEqual(alt["changed_aa"], "")
        self.assertEqual(wt["total_aa"], 7)
        self.assertEqual(alt["total_aa"], 5)
        # Visible flanks identical on both sides of the change.
        self.assertEqual(wt["visible_aa"], alt["visible_aa"])
        self.assertEqual(wt["trailing_visible_aa"], alt["trailing_visible_aa"])
        self.assertEqual(wt["trailing_visible_aa"], "FG")

    def test_pure_insertion_marks_wt_changed_empty(self):
        # WT: A B C F G   altered: A B C D E F G  (inserts "DE")
        wt, alt = _build_protein_windows(
            altered_aa="ABCDEFG", consensus_aa="ABCFG",
            frameshift=False, prefix_offset_aa=0, flank=3,
            stop_introduced=False,
        )
        self.assertEqual(wt["changed_aa"], "")        # `|` site marker
        self.assertEqual(alt["changed_aa"], "DE")
        self.assertEqual(wt["total_aa"], 5)
        self.assertEqual(alt["total_aa"], 7)

    def test_identical_returns_none(self):
        wt, alt = _build_protein_windows(
            altered_aa="ABCDEFG", consensus_aa="ABCDEFG",
            frameshift=False, prefix_offset_aa=0,
            stop_introduced=False,
        )
        self.assertIsNone(wt)
        self.assertIsNone(alt)

    def test_tandem_repeat_insertion_left_anchored(self):
        # WT: A B C D E F   altered: A B C D D E F  (insert 'D' adjacent to existing 'D').
        # The forward and reverse diffs cross because the inserted D matches
        # the existing D. Without re-anchoring, both `changed_aa` collapse to ''
        # and the inserted residue silently lands in trailing_visible_aa.
        # The fix shifts the change to a left-aligned canonical form: WT side
        # gets '' (rendered as `|` deletion-site marker), altered side gets the
        # inserted 'D' bracketed.
        wt, alt = _build_protein_windows(
            altered_aa="ACDDEF", consensus_aa="ACDEF",
            frameshift=False, prefix_offset_aa=0, flank=3,
            stop_introduced=False,
        )
        self.assertEqual(wt["changed_aa"], "")
        self.assertEqual(alt["changed_aa"], "D")
        # Visible flanks are now identical on both sides — the residues that
        # used to leak into trailing_visible_aa as 'DEF' on the altered side
        # are now correctly rendered as 'EF'.
        self.assertEqual(wt["trailing_visible_aa"], "EF")
        self.assertEqual(alt["trailing_visible_aa"], "EF")
        self.assertEqual(wt["total_aa"], 5)
        self.assertEqual(alt["total_aa"], 6)

    def test_tandem_repeat_deletion_left_anchored(self):
        # WT: A C E E E F G   altered: A C E E F G (delete one of three Es).
        # Symmetric ambiguity: the deleted E matches an adjacent E. Re-anchor
        # so WT shows the deleted E bracketed and altered gets the '' marker.
        wt, alt = _build_protein_windows(
            altered_aa="ACEEFG", consensus_aa="ACEEEFG",
            frameshift=False, prefix_offset_aa=0, flank=3,
            stop_introduced=False,
        )
        self.assertEqual(wt["changed_aa"], "E")
        self.assertEqual(alt["changed_aa"], "")
        self.assertEqual(wt["total_aa"], 7)
        self.assertEqual(alt["total_aa"], 6)

    def test_natural_stop_strip_applied(self):
        # Both end with '*' (natural stop). After stripping, sequences differ
        # only in the middle. In in-frame mode altered's '*' is also stripped
        # (frameshift=False).
        wt, alt = _build_protein_windows(
            altered_aa="ABCXYG*", consensus_aa="ABCDEFG*",
            frameshift=False, prefix_offset_aa=0, flank=3,
            stop_introduced=False,
        )
        self.assertEqual(wt["changed_aa"], "DEF")
        self.assertEqual(alt["changed_aa"], "XY")
        # totals exclude the stripped natural stop
        self.assertEqual(wt["total_aa"], 7)
        self.assertEqual(alt["total_aa"], 6)


if __name__ == "__main__":
    unittest.main()
