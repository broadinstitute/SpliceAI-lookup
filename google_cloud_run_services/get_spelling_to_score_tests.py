"""Tests for server.py's choice of which spelling of a variant to score.

A variant can be written with extra unchanged bases around the change: chr1-55057513-TG-TA and
chr1-55057511-GCTG-GCTA are both chr1-55057514-G-A. SpliceAI used to score those spellings
differently (https://github.com/broadinstitute/SpliceAI-lookup/issues/137), so server.py trims the
shared bases off before the cache lookup and before scoring. trim_shared_bases does the trimming, and
get_spelling_to_score checks the full REF against the reference first, so a wrong unchanged base is
still reported rather than trimmed away unread.

server.py can't be imported here (see the note in check_ref_allele_tests), so the functions under
test are lifted out of the source with ast and exec'd on their own, against small FASTA files built
in setUpClass.

Run with:  python3 -m unittest get_spelling_to_score_tests -v
"""

import os
import shutil
import tempfile
import unittest

from check_ref_allele_tests import load_functions_from_server_py

FUNCTIONS_UNDER_TEST = (
    "_get_fasta", "genome_display_name", "resolve_fasta_sequence_name", "check_ref_allele",
    "trim_shared_bases", "get_spelling_to_score",
)

# Bases 1-4 are GCTG, the unchanged bases and REF base of the issue #137 report, so the reporter's
# spellings can be replayed here at positions 1-4 instead of 55057511-55057514.
HG19_FASTA = ">1\nGCTGGTGAGC\n"
HG38_FASTA = ">chr1\nGCTGGTGAGC\n"


class TrimSharedBasesTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.namespace = load_functions_from_server_py({}, function_names=("trim_shared_bases",))

    def trim(self, pos, ref, alt):
        return self.namespace["trim_shared_bases"](pos, ref, alt)

    def test_snv_is_unchanged(self):
        self.assertEqual(self.trim(55057514, "G", "A"), (55057514, "G", "A"))

    def test_unchanged_bases_before_the_change_are_trimmed(self):
        self.assertEqual(self.trim(55057513, "TG", "TA"), (55057514, "G", "A"))
        self.assertEqual(self.trim(55057511, "GCTG", "GCTA"), (55057514, "G", "A"))

    def test_unchanged_bases_after_the_change_are_trimmed(self):
        self.assertEqual(self.trim(94555, "CGA", "TGA"), (94555, "C", "T"))

    def test_unchanged_bases_on_both_sides_are_trimmed(self):
        self.assertEqual(self.trim(99, "ACGT", "ACAT"), (101, "G", "A"))

    def test_genuine_mnv_is_unchanged(self):
        # both bases change, so there is nothing to trim
        self.assertEqual(self.trim(100, "TG", "CA"), (100, "TG", "CA"))

    def test_indels_with_one_anchor_base_are_unchanged(self):
        # both spellings are in expected_scores.json, whose cache keys must not move
        self.assertEqual(self.trim(1042466, "GGGC", "G"), (1042466, "GGGC", "G"))
        self.assertEqual(self.trim(1042601, "A", "AGAGAG"), (1042601, "A", "AGAGAG"))

    def test_indels_with_extra_shared_bases_keep_one_anchor_base(self):
        self.assertEqual(self.trim(55057512, "CTG", "CT"), (55057513, "TG", "T"))
        self.assertEqual(self.trim(100, "CA", "CAT"), (101, "A", "AT"))

    def test_shared_bases_at_the_end_are_trimmed_before_those_at_the_start(self):
        # CAA>CA could become CA>C or, one base later, AA>A. The SpliceAI fork's own trimming, which
        # this has to agree with, drops the end first.
        self.assertEqual(self.trim(100, "CAA", "CA"), (100, "CA", "C"))

    def test_deletion_insertion_trims_to_the_bases_it_changes(self):
        self.assertEqual(self.trim(100, "CAT", "CGGT"), (101, "A", "GG"))


class GetSpellingToScoreTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp_dir = tempfile.mkdtemp(prefix="get_spelling_to_score_tests_")
        paths = {}
        for genome_version, contents in (("37", HG19_FASTA), ("38", HG38_FASTA)):
            path = os.path.join(cls.tmp_dir, f"hg{genome_version}.fa")
            with open(path, "w") as f:
                f.write(contents)
            paths[genome_version] = path
        cls.namespace = load_functions_from_server_py(paths, function_names=FUNCTIONS_UNDER_TEST)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp_dir, ignore_errors=True)

    def spelling_to_score(self, chrom, pos, ref, alt, genome_version="38"):
        return self.namespace["get_spelling_to_score"](chrom, pos, ref, alt, genome_version)

    def test_padded_spellings_with_a_matching_ref_are_trimmed(self):
        self.assertEqual(self.spelling_to_score("chr1", 3, "TG", "TA"), (4, "G", "A"))
        self.assertEqual(self.spelling_to_score("1", 1, "GCTG", "GCTA"), (4, "G", "A"))
        self.assertEqual(self.spelling_to_score("chr1", 2, "CTG", "CT"), (3, "TG", "T"))

    def test_wrong_unchanged_base_is_left_for_the_ref_check_to_report(self):
        # the reference has T at 3, so trimming AG>AA down to G>A would hide the wrong A
        self.assertEqual(self.spelling_to_score("chr1", 3, "AG", "AA"), (3, "AG", "AA"))

    def test_spelling_with_nothing_to_trim_is_returned_as_given(self):
        # the reference has G at 4; the REF check in get_spliceai_scores reports this one
        self.assertEqual(self.spelling_to_score("chr1", 4, "T", "A"), (4, "T", "A"))
        self.assertEqual(self.spelling_to_score("chr1", 4, "G", "A"), (4, "G", "A"))

    def test_hg19_chromosome_names(self):
        self.assertEqual(self.spelling_to_score("1", 3, "TG", "TA", "37"), (4, "G", "A"))
        self.assertEqual(self.spelling_to_score("chr1", 3, "TG", "TA", "37"), (4, "G", "A"))
        self.assertEqual(self.spelling_to_score("1", 3, "AG", "AA", "37"), (3, "AG", "AA"))

    def test_missing_fasta_fails_open_and_trims(self):
        # like check_ref_allele, a REF that can't be checked doesn't block scoring
        namespace = load_functions_from_server_py(
            {"38": os.path.join(self.tmp_dir, "does_not_exist.fa")}, function_names=FUNCTIONS_UNDER_TEST)
        self.assertEqual(namespace["get_spelling_to_score"]("chr1", 3, "AG", "AA", "38"), (4, "G", "A"))


if __name__ == "__main__":
    unittest.main()
