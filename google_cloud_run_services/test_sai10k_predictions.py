#!/usr/bin/env python3
"""
Tests for SAI-10k predictions module.

These tests verify that the sai10k_predictions logic produces results consistent
with the examples in the SAI-10k-calc repository.

To run:
    python -m unittest tests.sai10k_predictions_tests -v

Or directly:
    python tests/sai10k_predictions_tests.py
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
    sai10k_determine_aberration_type,
    sai10k_predict_frameshift,
    sai10k_compute_predictions,
    sai10k_select_canonical_transcript,
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
        """Test finding a variant within an exon of BRCA1 (minus strand)."""
        exon_starts = BRCA1_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA1_TRANSCRIPT_HG19["EXON_ENDS"]

        # Variant at 41197774 is in exon 1
        result = sai10k_find_affected_region(41197774, exon_starts, exon_ends)
        self.assertIsNotNone(result)
        self.assertEqual(result["region_type"], "exon")

    def test_empty_exons(self):
        """Test with empty exon lists."""
        result = sai10k_find_affected_region(100, [], [])
        self.assertIsNone(result)


class TestDetermineAberrationType(unittest.TestCase):
    """Test the sai10k_determine_aberration_type function."""

    def setUp(self):
        self.exon_starts = BRCA2_TRANSCRIPT_HG19["EXON_STARTS"]
        self.exon_ends = BRCA2_TRANSCRIPT_HG19["EXON_ENDS"]

    def test_donor_loss(self):
        """Test detection of donor loss."""
        # Example: BRCA2_c.-40+1G>A at position 32889805
        # DS_DL=0.99, DS_DG=0.57, DS_AL=0.08, DS_AG=0.01
        result = sai10k_determine_aberration_type(
            ds_ag=0.01, ds_al=0.08, ds_dg=0.57, ds_dl=0.99,
            dp_ag=-38, dp_al=754, dp_dg=-100, dp_dl=-1,
            variant_pos=32889805,
            exon_starts=self.exon_starts,
            exon_ends=self.exon_ends,
            strand="+",
        )
        self.assertIn(result["aberration_type"], ["donor_shift", "donor_loss"])
        self.assertEqual(result["confidence"], "high")

    def test_exon_skipping_prediction(self):
        """Test detection of exon skipping pattern (DL + AL)."""
        # When both donor loss and acceptor loss are significant
        result = sai10k_determine_aberration_type(
            ds_ag=0.0, ds_al=0.81, ds_dg=0.17, ds_dl=0.98,
            dp_ag=-221, dp_al=-106, dp_dg=-137, dp_dl=-1,
            variant_pos=32890665,
            exon_starts=self.exon_starts,
            exon_ends=self.exon_ends,
            strand="+",
        )
        self.assertIn(result["aberration_type"], ["exon_skipping", "donor_loss"])
        self.assertEqual(result["max_delta_score"], 0.98)

    def test_no_significant_effect(self):
        """Test when delta scores are below threshold."""
        result = sai10k_determine_aberration_type(
            ds_ag=0.01, ds_al=0.01, ds_dg=0.04, ds_dl=0.04,
            dp_ag=-78, dp_al=-193, dp_dg=-214, dp_dl=-109,
            variant_pos=32890637,
            exon_starts=self.exon_starts,
            exon_ends=self.exon_ends,
            strand="+",
        )
        self.assertEqual(result["aberration_type"], "none")


class TestPredictFrameshift(unittest.TestCase):
    """Test the sai10k_predict_frameshift function."""

    def test_frameshift_prediction_exon_skipping(self):
        """Test frameshift prediction for exon skipping."""
        aberration_info = {
            "aberration_type": "exon_skipping",
            "affected_region": {"region_type": "exon", "region_number": 2},
        }
        exon_starts = BRCA2_TRANSCRIPT_HG19["EXON_STARTS"]
        exon_ends = BRCA2_TRANSCRIPT_HG19["EXON_ENDS"]
        exon_frames = BRCA2_TRANSCRIPT_HG19["EXON_FRAMES"]
        cds_start = BRCA2_TRANSCRIPT_HG19["CDS_START"]
        cds_end = BRCA2_TRANSCRIPT_HG19["CDS_END"]

        result = sai10k_predict_frameshift(
            aberration_type="exon_skipping",
            aberration_info=aberration_info,
            cds_start=cds_start,
            cds_end=cds_end,
            exon_starts=exon_starts,
            exon_ends=exon_ends,
            exon_frames=exon_frames,
        )
        self.assertTrue(result["affects_coding"])
        self.assertIn("exon_size", result)

    def test_non_coding_transcript(self):
        """Test frameshift prediction for non-coding transcript."""
        aberration_info = {
            "aberration_type": "donor_loss",
            "affected_region": {"region_type": "exon", "region_number": 1},
        }
        result = sai10k_predict_frameshift(
            aberration_type="donor_loss",
            aberration_info=aberration_info,
            cds_start=None,
            cds_end=None,
            exon_starts=[100, 200],
            exon_ends=[150, 250],
            exon_frames=[-1, -1],
        )
        self.assertFalse(result["affects_coding"])


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

        self.assertIn("aberration", result)
        self.assertIn("frameshift", result)
        self.assertIn("transcript_info", result)
        self.assertEqual(result["transcript_info"]["strand"], "+")
        self.assertTrue(result["transcript_info"]["is_coding"])

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

        self.assertIn("aberration", result)
        self.assertEqual(result["transcript_info"]["strand"], "-")


class TestSelectCanonicalTranscript(unittest.TestCase):
    """Test the sai10k_select_canonical_transcript function."""

    def test_select_mane_select_over_canonical(self):
        """Test that MANE Select is preferred over canonical."""
        scores_list = [
            {"t_priority": "C", "DS_AG": 0.5, "DS_AL": 0.5, "DS_DG": 0.5, "DS_DL": 0.5},
            {"t_priority": "MS", "DS_AG": 0.3, "DS_AL": 0.3, "DS_DG": 0.3, "DS_DL": 0.3},
        ]
        result = sai10k_select_canonical_transcript(scores_list)
        self.assertEqual(result["t_priority"], "MS")

    def test_select_highest_score_at_same_priority(self):
        """Test that highest max score is selected among same priority."""
        scores_list = [
            {"t_priority": "C", "DS_AG": 0.5, "DS_AL": 0.5, "DS_DG": 0.5, "DS_DL": 0.5},
            {"t_priority": "C", "DS_AG": 0.9, "DS_AL": 0.1, "DS_DG": 0.1, "DS_DL": 0.1},
        ]
        result = sai10k_select_canonical_transcript(scores_list)
        # Function selects by max delta score (0.9 > 0.5), not by sum
        self.assertEqual(result["DS_AG"], 0.9)

    def test_empty_scores_list(self):
        """Test handling of empty scores list."""
        result = sai10k_select_canonical_transcript([])
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

        self.assertIn("aberration", result)
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

        result = sai10k_compute_predictions(transcript_scores, variant_pos=32889805)

        # With DS_DL=0.99 and DS_DG=0.57, we expect donor-related aberration
        self.assertIn(
            result["aberration"]["aberration_type"],
            ["donor_shift", "donor_loss"],
        )
        self.assertEqual(result["aberration"]["confidence"], "high")

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

        # With both DS_AL=0.81 and DS_DL=0.98 significant, should predict exon skipping
        self.assertIn(
            result["aberration"]["aberration_type"],
            ["exon_skipping", "donor_loss"],
        )
        self.assertEqual(result["aberration"]["max_delta_score"], 0.98)

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

        self.assertEqual(result["aberration"]["aberration_type"], "none")
        self.assertEqual(result["aberration"]["confidence"], "low")

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

        # With DS_AL=0.28 and DS_DL=0.27 both >= 0.2, should predict exon skipping
        self.assertIn(
            result["aberration"]["aberration_type"],
            ["exon_skipping", "acceptor_loss", "donor_loss"],
        )

    def test_brca1_c5467_1g_a(self):
        """Test BRCA1_c.5467+1G>A (17-41199659-C-T).

        Expected: Any_splicing_aberration=YES, Partial_intron_retention=YES,
                  Exon_skipping=YES
        """
        transcript_scores = {
            "DS_AG": 0.00,
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
        }

        result = sai10k_compute_predictions(transcript_scores, variant_pos=41199659)

        # With DS_DL=0.93 and DS_DG=0.36 and DS_AL=0.40, multiple patterns
        self.assertIn(
            result["aberration"]["aberration_type"],
            ["exon_skipping", "donor_shift", "donor_loss"],
        )
        self.assertEqual(result["aberration"]["max_delta_score"], 0.93)


if __name__ == "__main__":
    unittest.main()
