"""Tests for server.py's check that a variant's REF allele fits the window being scored.

SpliceAI reports one score per position from `distance` bases before the variant to `distance` bases
after it, and it skips a record whose changed bases reach past that. The result is an empty score list,
which get_spliceai_scores would otherwise report as the variant not overlapping GENCODE's genes: the
wrong cause, and no help to someone whose deletion simply needs a larger distance. Alleles of the same
length are scored position by position, so they fit at any length, and the check has to let them
through the way the model's own span_fits_in_output_window does.

server.py can't be imported here (see the note in check_ref_allele_tests), so the functions under test
are lifted out of the source with ast and exec'd on their own.

Run with:  python3 -m unittest check_ref_allele_fits_window_tests -v
"""

import unittest

from check_ref_allele_tests import load_functions_from_server_py

FUNCTIONS_UNDER_TEST = ("trim_shared_bases", "check_ref_allele_fits_window")
CONSTANTS_UNDER_TEST = ("MAX_DISTANCE_LIMIT", "DEFAULT_DISTANCE")


class CheckRefAlleleFitsWindowTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.namespace = load_functions_from_server_py(
            {}, function_names=FUNCTIONS_UNDER_TEST, constant_names=CONSTANTS_UNDER_TEST)
        cls.max_distance = cls.namespace["MAX_DISTANCE_LIMIT"]
        cls.default_distance = cls.namespace["DEFAULT_DISTANCE"]

    def check(self, ref, alt, distance):
        return self.namespace["check_ref_allele_fits_window"](ref, alt, distance)

    # --- the changed bases fit, so scoring goes ahead ---

    def test_single_base_variants_fit(self):
        self.assertIsNone(self.check("G", "A", self.default_distance))
        self.assertIsNone(self.check("G", "GAGA", self.default_distance))
        self.assertIsNone(self.check("GAGA", "G", self.default_distance))

    def test_deletion_reaching_the_last_scored_position_fits(self):
        # the window runs `distance` bases past the variant, so a REF of distance + 1 ends on its last base
        self.assertIsNone(self.check("A" * (self.default_distance + 1), "A", self.default_distance))

    def test_alleles_of_the_same_length_fit_at_any_length(self):
        # each position has its own score, so nothing has to be collapsed and no length runs past the window
        self.assertIsNone(self.check("A" * 600, "C" * 600, self.default_distance))
        self.assertIsNone(self.check("A" * 600, "C" * 600, 50))

    def test_shared_bases_are_trimmed_before_measuring(self):
        # "A"*600 + "GT" > "A"*600 + "GC" changes one base, however it is written
        self.assertIsNone(self.check("A" * 600 + "GT", "A" * 600 + "GC", self.default_distance))

    def test_a_smaller_distance_fits_a_smaller_deletion(self):
        self.assertIsNone(self.check("A" * 51, "A", 50))
        self.assertIsNotNone(self.check("A" * 52, "A", 50))

    # --- the changed bases run past the window, so the user is told why and what to retry with ---

    def test_deletion_one_base_past_the_window_is_reported(self):
        message = self.check("A" * (self.default_distance + 2), "A", self.default_distance)
        self.assertIsNotNone(message)
        self.assertIn("502 bases long", message)
        self.assertIn("500 bases on either side", message)

    def test_message_suggests_a_distance_that_would_fit(self):
        message = self.check("A" * 600, "A", self.default_distance)
        self.assertIn("Retry with distance=599 or more", message)
        # the suggestion has to be usable in a URL, so no thousands separator
        self.assertNotIn("distance=5,99", message)

    def test_deletion_insertion_past_the_window_is_reported(self):
        self.assertIsNotNone(self.check("A" * 600, "GC", self.default_distance))

    def test_the_largest_distance_the_api_accepts_is_still_suggested(self):
        # the API rejects only distances above the limit, so a REF needing exactly it can be retried
        self.assertIn(f"Retry with distance={self.max_distance} or more",
                      self.check("A" * (self.max_distance + 1), "A", self.default_distance))

    def test_a_ref_past_what_any_distance_can_reach_is_reported_as_unscoreable(self):
        message = self.check("A" * (self.max_distance + 2), "A", self.default_distance)
        self.assertNotIn("Retry with", message)
        self.assertIn("10,000 bases on either side", message)


if __name__ == "__main__":
    unittest.main()
