"""Tests for server.py's mitochondrial annotation/FASTA name alignment.

hg19 named the mitochondrion three different ways across the files the services read: the FASTA
calls it "MT", the SpliceAI annotation "chrM", and the Pangolin database "M". SpliceAI passes one
chromosome name to both its annotation lookup and its FASTA read, through a normalise_chrom() that
only adds or strips a "chr" prefix, so on hg19 no single value satisfied both and every
mitochondrial variant came back as "no scores". align_annotator_mito_chrom renames the annotation's
mito rows in memory to the FASTA's spelling so that one value does.

server.py can't be imported here (see the note in check_ref_allele_tests), so the function under
test is lifted out of the source with ast and exec'd on its own, against a stand-in annotator
holding a numpy CHROM array shaped like the one Annotator builds with
pandas.read_csv(dtype={'CHROM': object}).

Run with:  python3 -m unittest align_annotator_mito_chrom_tests -v
"""

import unittest

import numpy as np

from check_ref_allele_tests import load_functions_from_server_py

FUNCTIONS_UNDER_TEST = ("align_annotator_mito_chrom",)


class FakeAnnotator:
    """Stands in for spliceai.utils.Annotator, which needs a real FASTA and annotation file."""

    def __init__(self, chroms, dtype=object):
        self.chroms = np.array(chroms, dtype=dtype)


def normalise_chrom(source, target):
    """Verbatim copy of spliceai.utils.normalise_chrom (pinned fork 7f36ca84).

    Copied rather than imported because the spliceai package only exists inside the containers.
    The tests below use it to assert the end-to-end property that actually matters: that ONE
    chromosome name now resolves to both the annotation's spelling and the FASTA's.
    """
    def has_prefix(x):
        return x.startswith('chr')

    if has_prefix(source) and not has_prefix(target):
        return source.strip('chr')
    elif not has_prefix(source) and has_prefix(target):
        return 'chr' + source
    return source


class AlignAnnotatorMitoChromTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.namespace = load_functions_from_server_py({}, function_names=FUNCTIONS_UNDER_TEST)

    def align(self, annotator, mito_fasta_name):
        return self.namespace["align_annotator_mito_chrom"](annotator, mito_fasta_name)

    # --- hg19: the build that was broken ---

    def test_hg19_annotation_is_renamed_to_the_fasta_spelling(self):
        annotator = FakeAnnotator(["chr1", "chr1", "chrM", "chrM", "chr2"])
        self.align(annotator, "MT")
        self.assertEqual(list(annotator.chroms), ["chr1", "chr1", "chrMT", "chrMT", "chr2"])

    def test_hg19_one_name_now_satisfies_both_lookups(self):
        # The property the fix exists for. Before it, "MT" normalised to "chrMT" against the
        # annotation (which held "chrM") and found nothing, while "M" would have normalised to
        # "M" against the FASTA (which holds "MT") and failed the read instead.
        annotator = FakeAnnotator(["chr1", "chrM"])
        self.align(annotator, "MT")

        chrom = "MT"  # what MITO_CHROM_NAME remaps the user's chromosome to on hg19
        annotation_target = str(annotator.chroms[0])          # get_name_and_strand normalises to this
        fasta_target = "1"                                    # hg19 FASTA keys are unprefixed
        self.assertEqual(normalise_chrom(chrom, annotation_target), "chrMT")
        self.assertIn("chrMT", [str(c) for c in annotator.chroms])
        self.assertEqual(normalise_chrom(chrom, fasta_target), "MT")

    def test_the_pre_fix_name_would_have_missed(self):
        # Guards the reason the bug was invisible: "chrM" is what the file ships, and "MT"
        # normalised against a prefixed annotation never equals it.
        self.assertNotEqual(normalise_chrom("MT", "chr1"), "chrM")

    # --- hg38: already consistent, must not be disturbed ---

    def test_hg38_is_a_no_op(self):
        annotator = FakeAnnotator(["chr1", "chrM", "chr2"])
        before = list(annotator.chroms)
        self.align(annotator, "M")
        self.assertEqual(list(annotator.chroms), before)

    def test_non_mitochondrial_rows_are_never_touched(self):
        annotator = FakeAnnotator(["chr1", "chrM", "chrX", "chrY", "chr22"])
        self.align(annotator, "MT")
        self.assertEqual([str(c) for c in annotator.chroms if not str(c).upper().endswith("MT")],
                         ["chr1", "chrX", "chrY", "chr22"])

    # --- naming conventions the annotation could plausibly use ---

    def test_unprefixed_annotation_keeps_its_convention(self):
        # An annotation whose first row is bare must not acquire a "chr" on its mito rows,
        # since get_name_and_strand normalises the incoming name against that first row.
        annotator = FakeAnnotator(["1", "M", "2"])
        self.align(annotator, "MT")
        self.assertEqual(list(annotator.chroms), ["1", "MT", "2"])

    def test_every_mito_spelling_is_recognised(self):
        annotator = FakeAnnotator(["chr1", "chrM", "chrMT", "M", "MT"])
        self.align(annotator, "MT")
        self.assertEqual([str(c) for c in annotator.chroms][1:], ["chrMT"] * 4)

    # --- the trap that would silently restore the bug ---

    def test_a_fixed_width_array_does_not_truncate_the_longer_name(self):
        # "chrM" -> "chrMT" grows by one character. Assigning into a '<U4' numpy array would
        # silently write back "chrM" and leave the bug in place with the fix apparently applied.
        annotator = FakeAnnotator(["chr1", "chrM"], dtype="<U4")
        self.align(annotator, "MT")
        self.assertEqual(str(annotator.chroms[1]), "chrMT")

    # --- degenerate inputs ---

    def test_no_mito_sequence_in_the_fasta_is_left_alone(self):
        annotator = FakeAnnotator(["chr1", "chrM"])
        self.align(annotator, None)
        self.assertEqual(list(annotator.chroms), ["chr1", "chrM"])

    def test_annotation_without_mito_rows_is_left_alone(self):
        annotator = FakeAnnotator(["chr1", "chr2"])
        self.align(annotator, "MT")
        self.assertEqual(list(annotator.chroms), ["chr1", "chr2"])

    def test_empty_annotation_does_not_raise(self):
        annotator = FakeAnnotator([])
        self.align(annotator, "MT")
        self.assertEqual(list(annotator.chroms), [])

    def test_running_twice_is_idempotent(self):
        annotator = FakeAnnotator(["chr1", "chrM"])
        self.align(annotator, "MT")
        self.align(annotator, "MT")
        self.assertEqual(list(annotator.chroms), ["chr1", "chrMT"])


if __name__ == "__main__":
    unittest.main()
