"""Tests for server.py's model-context-window check.

Both SpliceAI and Pangolin slice a fixed flank of reference sequence around the variant straight
out of pyfastx with no bounds check of their own. A window that runs off the start of a contig
gives pyfastx a negative start index, which it answers by segfaulting -- taking the gunicorn
worker with it and turning the request into a 503. check_model_context_window rejects those
before either model runs. The boundaries asserted here were measured against production:
chr1-11500-C-T scored at distance=6499 and 503'd at distance=6500 on all four services.

server.py can't be imported here (see the note in check_ref_allele_tests), so the function under
test is lifted out of the source with ast and exec'd on its own, against small FASTA files built
in setUpClass.

Run with:  python3 -m unittest check_model_context_window_tests -v
"""

import os
import shutil
import tempfile
import unittest

from check_ref_allele_tests import load_functions_from_server_py

FUNCTIONS_UNDER_TEST = ("_get_fasta", "resolve_fasta_sequence_name", "check_model_context_window")
CONSTANTS_UNDER_TEST = ("MODEL_FLANK_SIZE",)

# Real contig lengths matter here, so these are built to size rather than written out literally.
# chrM is the length of the real hg38 mitochondrion, which is shorter than the models' own
# window: no distance setting makes m.4092 scoreable, and that is the case that was crashing.
CHRM_LENGTH = 16569
CHR1_LENGTH = 40000
MT_LENGTH = 16571  # hg19's mitochondrion, under its hg19 name


def write_fasta(path, contigs):
    """Write a FASTA of `contigs` ({name: length}), filled with a repeating ACGT."""
    with open(path, "w") as f:
        for name, length in contigs.items():
            f.write(f">{name}\n")
            sequence = ("ACGT" * (length // 4 + 1))[:length]
            for i in range(0, length, 60):
                f.write(sequence[i:i + 60] + "\n")


class CheckModelContextWindowTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp_dir = tempfile.mkdtemp(prefix="check_model_context_window_tests_")
        paths = {}
        for genome_version, contigs in (
                ("37", {"1": CHR1_LENGTH, "MT": MT_LENGTH}),
                ("38", {"chr1": CHR1_LENGTH, "chrM": CHRM_LENGTH})):
            path = os.path.join(cls.tmp_dir, f"hg{genome_version}.fa")
            write_fasta(path, contigs)
            paths[genome_version] = path
        cls.namespace = load_functions_from_server_py(
            paths, function_names=FUNCTIONS_UNDER_TEST, constant_names=CONSTANTS_UNDER_TEST)
        cls.flank = cls.namespace["MODEL_FLANK_SIZE"]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp_dir, ignore_errors=True)

    def check(self, chrom, pos, ref_length, distance, genome_version="38"):
        return self.namespace["check_model_context_window"](chrom, pos, ref_length, distance, genome_version)

    # --- the start of a contig: the case that segfaults ---

    def test_boundary_at_the_start_of_a_contig(self):
        # the slice start is pos-MODEL_FLANK_SIZE-distance-1, so the last safe position is
        # exactly MODEL_FLANK_SIZE+distance+1. Measured in production as distance=6499 -> 200
        # and distance=6500 -> 503 for chr1:11500.
        self.assertIsNone(self.check("chr1", 11500, 1, 6498))
        self.assertIsNone(self.check("chr1", 11500, 1, 6499))
        self.assertIsNotNone(self.check("chr1", 11500, 1, 6500))

    def test_boundary_holds_at_distance_zero(self):
        last_bad = self.flank
        self.assertIsNotNone(self.check("chr1", last_bad, 1, 0))
        self.assertIsNone(self.check("chr1", last_bad + 1, 1, 0))

    def test_mitochondrial_variant_at_the_default_distance(self):
        # MT-ND1 spans m.3306-4262, so this is a real variant a user can type. It segfaulted
        # every hg38 service at the default distance of 500.
        message = self.check("M", 4092, 1, 500)
        self.assertIsNotNone(message)
        self.assertIn("M-4092", message)
        self.assertIn("start of M", message)

    def test_no_workable_distance_is_reported_as_such(self):
        # m.4092 is inside the fixed flank, so no distance setting rescues it
        message = self.check("M", 4092, 1, 500)
        self.assertIn("No distance setting is small enough", message)
        self.assertNotIn("Retry with distance=", message)

    def test_a_workable_distance_is_suggested_and_is_the_real_boundary(self):
        message = self.check("chr1", 11500, 1, 10000)
        self.assertIn("Retry with distance=6499 or less.", message)
        # the suggestion has to actually pass the check it came from
        self.assertIsNone(self.check("chr1", 11500, 1, 6499))

    def test_message_names_the_distance_that_was_asked_for(self):
        message = self.check("chr1", 11500, 1, 10000)
        self.assertIn("15,000bp of sequence on each side", message)
        self.assertIn("distance setting of 10,000bp", message)

    # --- the end of a contig: no crash, just a clearer message than either tool gives ---

    def test_boundary_at_the_end_of_a_contig(self):
        # the read runs to pos+ref_length+MODEL_FLANK_SIZE+distance-1
        last_good = CHRM_LENGTH - self.flank - 50
        self.assertIsNone(self.check("chrM", last_good, 1, 50))
        self.assertIsNotNone(self.check("chrM", last_good + 1, 1, 50))

    def test_a_longer_ref_allele_needs_one_more_base(self):
        # Pangolin's process_variant reads to pos+len(ref)+4999+distance, so a 4-base deletion
        # runs 3 bases further than a SNV at the same position
        last_good_for_a_snv = CHRM_LENGTH - self.flank - 50
        self.assertIsNone(self.check("chrM", last_good_for_a_snv, 1, 50))
        self.assertIsNotNone(self.check("chrM", last_good_for_a_snv, 4, 50))
        self.assertIsNone(self.check("chrM", last_good_for_a_snv - 3, 4, 50))

    def test_end_message_suggests_a_workable_distance(self):
        message = self.check("chrM", 11000, 1, 1000)
        self.assertIsNotNone(message)
        self.assertIn("end of chrM", message)
        max_distance = CHRM_LENGTH - 11000 - self.flank
        self.assertIn(f"Retry with distance={max_distance} or less.", message)
        self.assertIsNone(self.check("chrM", 11000, 1, max_distance))

    def test_end_reports_no_workable_distance_when_the_fixed_flank_alone_overruns(self):
        # only 4,569 bases follow m.12000, fewer than the fixed flank, so no distance helps
        message = self.check("chrM", 12000, 1, 50)
        self.assertIn("end of chrM", message)
        self.assertIn("No distance setting is small enough", message)

    # --- a window that fits is left alone ---

    def test_a_window_well_inside_a_contig_is_accepted(self):
        self.assertIsNone(self.check("chrM", 8600, 1, 50))
        self.assertIsNone(self.check("chr1", 20000, 1, 500))

    # --- chromosome naming, matching check_ref_allele's behaviour ---

    def test_hg19_names_are_accepted_for_the_end_check(self):
        self.assertIsNone(self.check("MT", 8600, 1, 50, genome_version="37"))
        self.assertIsNotNone(self.check("MT", 16000, 1, 50, genome_version="37"))

    def test_mitochondrial_alias_resolves_for_the_end_check(self):
        # hg38's FASTA calls it "chrM"; a user typing the hg19 spelling still gets the length
        self.assertIsNotNone(self.check("MT", 12000, 1, 50))
        self.assertIsNotNone(self.check("chrMT", 12000, 1, 50))

    # --- failing open, and the one place it must not ---

    def test_unknown_contig_fails_open_at_the_end_but_not_at_the_start(self):
        # nothing is known about how long chr99 is, so the end check can't speak to it
        self.assertIsNone(self.check("chr99", 30000, 1, 500))
        # ...but the start check is pure arithmetic and still has to reject
        self.assertIsNotNone(self.check("chr99", 100, 1, 500))

    def test_missing_fasta_still_rejects_the_crashing_case(self):
        # a FASTA that failed to open must not turn the guard off and hand the segfault back
        namespace = load_functions_from_server_py(
            {"38": os.path.join(self.tmp_dir, "does_not_exist.fa")},
            function_names=FUNCTIONS_UNDER_TEST, constant_names=CONSTANTS_UNDER_TEST)
        self.assertIsNotNone(namespace["check_model_context_window"]("chrM", 4092, 1, 500, "38"))
        self.assertIsNone(namespace["check_model_context_window"]("chr1", 30000, 1, 500, "38"))


if __name__ == "__main__":
    unittest.main()
