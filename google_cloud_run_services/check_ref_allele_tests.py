"""Tests for server.py's REF-allele check.

server.py can't be imported here: at import time it reads GENOME_VERSION/TOOL from the
environment and pulls in the spliceai or pangolin packages, which only exist inside the
containers. So the three functions under test are lifted out of the source with ast and
exec'd on their own, against small FASTA files built in setUpClass. That runs the real
source text rather than a copy of it, so the tests still fail if server.py changes.

Run with:  python3 -m unittest check_ref_allele_tests -v
"""

import ast
import os
import shutil
import tempfile
import textwrap
import threading
import unittest

SERVER_PY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "server.py")
FUNCTIONS_UNDER_TEST = ("_get_fasta", "genome_display_name", "check_ref_allele")

# hg19's FASTA names its sequences without a "chr" prefix and calls the mitochondrion "MT";
# hg38's uses "chr1".."chrM". The check has to cope with both, so each build gets its own file.
HG19_FASTA = ">1\nACGTACGTAC\nGTACGTACGT\n>2\nTTTTAAAACC\n>MT\nGATCACAGGT\n"
HG38_FASTA = ">chr1\nACGTACGTAC\nGTACGTACGT\n>chr2\nTTTTAAAACC\n>chrM\nGATCACAGGT\n"
# Soft-masked (lowercase) reference bases must still compare equal to an uppercase REF.
HG38_SOFTMASKED_CONTIG = ">chr3\nacgtacgtac\n"
# Both real FASTAs hard-mask the chrY pseudoautosomal regions to N, and GENCODE transcripts start
# inside them, so an N reference base has to fail open rather than be reported as a mismatch.
HG38_N_MASKED_CONTIG = ">chr4\nNNNNNNNNNN\n>chr5\nACGNNACGTA\n"


def load_functions_from_server_py(fasta_path_by_genome_version):
    """Exec the functions under test in a namespace of their own.

    Args:
        fasta_path_by_genome_version (dict): what server.py's FASTA_PATH should resolve to,
            keyed by genome version ("37"/"38")

    Return:
        dict: the exec'd namespace, holding the functions plus their module-level state
    """
    with open(SERVER_PY) as f:
        source = f.read()
    tree = ast.parse(source, filename=SERVER_PY)
    sources = {
        node.name: ast.get_source_segment(source, node)
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in FUNCTIONS_UNDER_TEST
    }
    missing = set(FUNCTIONS_UNDER_TEST) - set(sources)
    if missing:
        raise AssertionError(f"server.py no longer defines: {sorted(missing)}")

    namespace = {
        "FASTA_PATH": fasta_path_by_genome_version,
        "FASTA": {},
        "_FASTA_LOCK": threading.Lock(),
        "threading": threading,
    }
    for name in FUNCTIONS_UNDER_TEST:
        exec(sources[name], namespace)
    return namespace


class CheckRefAlleleTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp_dir = tempfile.mkdtemp(prefix="check_ref_allele_tests_")
        paths = {}
        for genome_version, contents in (("37", HG19_FASTA), ("38", HG38_FASTA + HG38_SOFTMASKED_CONTIG + HG38_N_MASKED_CONTIG)):
            path = os.path.join(cls.tmp_dir, f"hg{genome_version}.fa")
            with open(path, "w") as f:
                f.write(contents)
            paths[genome_version] = path
        cls.namespace = load_functions_from_server_py(paths)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp_dir, ignore_errors=True)

    def check(self, chrom, pos, ref, genome_version):
        return self.namespace["check_ref_allele"](chrom, pos, ref, genome_version)

    # --- the REF matches, so scoring goes ahead ---

    def test_matching_ref_returns_no_error(self):
        self.assertIsNone(self.check("chr1", 1, "A", "38"))
        self.assertIsNone(self.check("chr1", 5, "A", "38"))

    def test_matching_multi_base_ref_returns_no_error(self):
        # a deletion's REF spans several bases; all of them have to match.
        # chr1 is ACGT repeated, so 1-4 is "ACGT" and 3-6 is "GTAC".
        self.assertIsNone(self.check("chr1", 1, "ACGT", "38"))
        self.assertIsNone(self.check("chr1", 3, "GTAC", "38"))
        # spanning the FASTA's line break at base 10 must still read straight through
        self.assertIsNone(self.check("chr1", 9, "ACGT", "38"))

    def test_lowercase_reference_bases_still_match(self):
        # soft-masked regions are stored lowercase, but the user types uppercase
        self.assertIsNone(self.check("chr3", 1, "ACGT", "38"))

    def test_lowercase_user_ref_still_matches(self):
        self.assertIsNone(self.check("chr1", 1, "acgt", "38"))

    # --- the REF is wrong, so the user gets told which base is actually there ---

    def test_wrong_single_base_ref_is_reported(self):
        message = self.check("chr1", 1, "T", "38")
        self.assertIsNotNone(message)
        self.assertIn("unexpected reference allele", message)
        self.assertIn("is A, not T", message)
        self.assertIn("hg38", message)

    def test_wrong_multi_base_ref_is_reported(self):
        message = self.check("chr1", 1, "ACGA", "38")
        self.assertIsNotNone(message)
        self.assertIn("is ACGT, not ACGA", message)

    def test_message_names_the_variant_and_position(self):
        message = self.check("chr1", 5, "G", "38")
        self.assertIn("chr1-5-G", message)
        self.assertIn("chr1:5", message)

    def test_hg19_message_says_hg19_not_hg37(self):
        message = self.check("1", 1, "T", "37")
        self.assertIn("hg19", message)
        self.assertNotIn("hg37", message)

    # --- chromosome naming differs between the two builds ---

    def test_hg19_bare_chromosome_name(self):
        self.assertIsNone(self.check("1", 1, "A", "37"))
        self.assertIn("is A, not T", self.check("1", 1, "T", "37"))

    def test_hg19_accepts_a_chr_prefix_the_fasta_does_not_have(self):
        # the front end and the API both allow "chr1" on hg19, where the FASTA says "1"
        self.assertIsNone(self.check("chr1", 1, "A", "37"))
        self.assertIn("is A, not T", self.check("chr1", 1, "T", "37"))

    def test_hg38_accepts_a_bare_name_the_fasta_prefixes(self):
        self.assertIsNone(self.check("1", 1, "A", "38"))
        self.assertIn("is A, not T", self.check("1", 1, "T", "38"))

    def test_mitochondrion_under_each_builds_name(self):
        self.assertIsNone(self.check("MT", 1, "G", "37"))
        self.assertIsNone(self.check("chrM", 1, "G", "38"))

    def test_mitochondrion_m_mt_alias_resolves_on_both_builds(self):
        # The bare name differs between the builds ("MT" in hg19, "M" in hg38), and callers pass
        # either. get_pangolin_scores has no mito remap of its own, so check_ref_allele has to
        # resolve the alias itself or it would silently skip validation for the Pangolin service.
        self.assertIsNone(self.check("M", 1, "G", "37"))     # hg19 FASTA calls it "MT"
        self.assertIsNone(self.check("chrM", 1, "G", "37"))
        self.assertIsNone(self.check("MT", 1, "G", "38"))    # hg38 FASTA calls it "chrM"
        self.assertIsNone(self.check("chrMT", 1, "G", "38"))

    def test_mitochondrion_wrong_ref_is_reported_under_every_alias(self):
        for chrom, genome_version in (("M", "37"), ("MT", "37"), ("chrM", "37"),
                                      ("M", "38"), ("MT", "38"), ("chrM", "38")):
            with self.subTest(chrom=chrom, genome_version=genome_version):
                message = self.check(chrom, 1, "A", genome_version)
                self.assertIsNotNone(message, f"{chrom} on hg{genome_version} skipped validation")
                self.assertIn("is G, not A", message)

    # --- fails open: anything it can't check must not block scoring ---

    def test_unknown_contig_returns_no_error(self):
        self.assertIsNone(self.check("chr99", 1, "A", "38"))

    def test_position_past_the_end_of_the_contig_returns_no_error(self):
        # chr2 is 10bp; a REF running off the end can't be compared, so leave it to the model
        self.assertIsNone(self.check("chr2", 10, "CGG", "38"))
        self.assertIsNone(self.check("chr2", 500, "A", "38"))

    def test_wildly_out_of_range_position_returns_no_error(self):
        # VARIANT_RE accepts positions up to 999,999,999 from a public unauthenticated endpoint, and
        # pyfastx segfaults rather than raising on many out-of-range fetches, so the bound check has
        # to happen before .fetch(). Without it this test kills the interpreter instead of failing.
        for pos in (10 ** 6, 4 * 10 ** 8, 9 * 10 ** 8, 999999999):
            with self.subTest(pos=pos):
                self.assertIsNone(self.check("chr1", pos, "A", "38"))
                self.assertIsNone(self.check("1", pos, "A", "37"))

    def test_n_masked_reference_returns_no_error(self):
        # reporting "the reference allele is N" would be both wrong and blocking
        self.assertIsNone(self.check("chr4", 1, "A", "38"))
        self.assertIsNone(self.check("chr4", 1, "ACGT", "38"))

    def test_ref_spanning_into_an_n_run_returns_no_error(self):
        # chr5 is ACGNNACGTA: a REF overlapping the N run can't be compared either
        self.assertIsNone(self.check("chr5", 3, "GN", "38"))
        self.assertIsNone(self.check("chr5", 1, "ACGTT", "38"))

    def test_ref_beside_an_n_run_is_still_checked(self):
        # the first three bases of chr5 are real, so a mismatch there must still be reported
        self.assertIsNone(self.check("chr5", 1, "ACG", "38"))
        self.assertIn("is ACG, not TTT", self.check("chr5", 1, "TTT", "38"))

    def test_zero_and_negative_positions_return_no_error(self):
        self.assertIsNone(self.check("chr1", 0, "A", "38"))
        self.assertIsNone(self.check("chr1", -5, "A", "38"))

    def test_missing_fasta_returns_no_error(self):
        namespace = load_functions_from_server_py({"38": os.path.join(self.tmp_dir, "does_not_exist.fa")})
        self.assertIsNone(namespace["check_ref_allele"]("chr1", 1, "T", "38"))

    # --- the display name used in messages ---

    def test_genome_display_name(self):
        self.assertEqual(self.namespace["genome_display_name"]("37"), "hg19")
        self.assertEqual(self.namespace["genome_display_name"]("38"), "hg38")


if __name__ == "__main__":
    unittest.main()
