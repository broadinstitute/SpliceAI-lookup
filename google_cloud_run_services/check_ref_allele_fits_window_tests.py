"""Tests for server.py's checks that a variant's REF allele fits the window being scored.

The two tools set different limits, so each has its own check here.

SpliceAI reports one score per position from `distance` bases before the variant to `distance` bases
after it, and it skips a record whose changed bases reach past that. The result is an empty score list,
which get_spliceai_scores would otherwise report as the variant not overlapping GENCODE's genes: the
wrong cause, and no help to someone whose deletion simply needs a larger distance. Alleles of the same
length are scored position by position, so they fit at any length, and the check has to let them
through the way the model's own span_fits_in_output_window does.

Pangolin's window instead runs `distance` bases past the END of the REF allele, so the bases a variant
changes always fit and the only limit is on how long the REF allele itself may be.

server.py can't be imported here (see the note in check_ref_allele_tests), so the functions under test
are lifted out of the source with ast and exec'd on their own.

Run with:  python3 -m unittest check_ref_allele_fits_window_tests -v
"""

import unittest

from check_ref_allele_tests import load_functions_from_server_py

FUNCTIONS_UNDER_TEST = ("trim_shared_bases", "check_ref_allele_fits_window",
                        "check_ref_allele_fits_pangolin_window")
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


class CheckRefAlleleFitsPangolinWindowTest(unittest.TestCase):
    """Pangolin's limit is on the REF allele's own length, not on where the changed bases land.

    Its process_variant skips a REF longer than twice the distance ("Deletion too large"), which reaches
    the caller as an empty result and gets reported as the model simply having no scores for the variant.
    """

    @classmethod
    def setUpClass(cls):
        cls.namespace = load_functions_from_server_py(
            {}, function_names=FUNCTIONS_UNDER_TEST, constant_names=CONSTANTS_UNDER_TEST)
        cls.max_distance = cls.namespace["MAX_DISTANCE_LIMIT"]
        cls.default_distance = cls.namespace["DEFAULT_DISTANCE"]

    def check(self, ref, distance):
        return self.namespace["check_ref_allele_fits_pangolin_window"](ref, distance)

    # --- the REF is short enough, so scoring goes ahead ---

    def test_ordinary_variants_fit(self):
        self.assertIsNone(self.check("G", self.default_distance))
        self.assertIsNone(self.check("GAGA", self.default_distance))

    def test_a_ref_of_exactly_twice_the_distance_fits(self):
        self.assertIsNone(self.check("A" * 2 * self.default_distance, self.default_distance))
        self.assertIsNotNone(self.check("A" * (2 * self.default_distance + 1), self.default_distance))

    def test_the_changed_bases_never_have_to_fit(self):
        # unlike SpliceAI, the window runs past the end of the REF allele, so a deletion-insertion
        # whose changed bases sit at the far end of a long REF is still scored
        self.assertIsNone(self.check("A" * 600, 500))

    # --- the REF is too long, so the user is told why and what to retry with ---

    def test_a_ref_past_the_limit_is_reported(self):
        message = self.check("A" * 1001, self.default_distance)
        self.assertIn("1,001 bases long", message)
        self.assertIn("500 bases on either side", message)

    def test_the_suggested_distance_is_the_smallest_one_that_works(self):
        for length in (1001, 1002):
            with self.subTest(length=length):
                message = self.check("A" * length, self.default_distance)
                suggested = int(message.split("distance=")[1].split(" ")[0])
                self.assertIsNone(self.check("A" * length, suggested))
                self.assertIsNotNone(self.check("A" * length, suggested - 1))

    def test_the_suggestion_is_usable_in_a_url(self):
        # a thousands separator would make the suggested distance invalid as a query parameter
        self.assertIn("Retry with distance=5001 or more", self.check("A" * 10001, self.default_distance))

    def test_the_largest_distance_the_api_accepts_is_still_suggested(self):
        self.assertIn(f"Retry with distance={self.max_distance} or more",
                      self.check("A" * (2 * self.max_distance), self.default_distance))

    def test_a_ref_past_what_any_distance_can_reach_is_reported_as_unscoreable(self):
        message = self.check("A" * (2 * self.max_distance + 1), self.default_distance)
        self.assertNotIn("Retry with", message)
        self.assertIn("10,000 bases on either side", message)


if __name__ == "__main__":
    unittest.main()
