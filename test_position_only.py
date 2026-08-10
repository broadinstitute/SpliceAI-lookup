"""
Position-only (REF-only) tests for SpliceAI Lookup.

A bare "chrom-pos" query (no ref/alt) asks for SpliceAI/Pangolin scores computed
on the reference sequence alone, instead of the usual REF-vs-ALT delta scores.
This module covers the two parts of that path that nothing else tests:

  1. server.py's pure helpers -- parse_position / parse_variant precedence and
     get_splicing_scores_cache_key's is_position_only handling. These import
     google_cloud_run_services/server.py, which reads TOOL / GENOME_VERSION at
     import time and imports the matching model package (spliceai or pangolin).
     They skip when neither package is installed, but fail if one is installed
     and too old to import -- a silent skip there would hide a real regression.

  2. The deployed spliceai/pangolin endpoints' REF-only response shape.

To run:
    python3 -m unittest test_position_only -v

The live API tests use the same endpoint configuration as test_api_consistency.py
and default to the production Cloud Run revisions. Until a non-"dev" tag has
deployed REF-only support to production, run them against the dev revisions:

    SPLICEAI_API_ENV=dev python3 -m unittest test_position_only -v
"""

import importlib
import os
import sys
import unittest
from unittest import mock

import requests

from test_api_consistency import API_URLS


# A locus with strong reference acceptor/donor predictions on both tools
# (TRAPPC9, hg38) -- the position half of the chr8-140300616-T-G variant used
# throughout test_api_consistency.py and index.html's examples.
TEST_POSITION = "8-140300616"
TEST_VARIANT = "8-140300616-T-G"


def _import_server():
    """Import google_cloud_run_services/server.py for unit-testing its pure helpers.

    server.py raises at import unless GENOME_VERSION and TOOL are set, and it
    imports the model package named by TOOL. The helpers under test are
    tool-independent, so this tries each tool and keeps whichever import works.

    Only a genuinely absent model package counts as "can't run these tests". A
    package that is installed but too old to provide what server.py imports, and
    any other import-time failure (syntax error, failed env-var check, a defect in
    server.py itself), is raised instead -- otherwise a real regression is reported
    as a skip and nobody notices the tests stopped running.

    Returns:
        The imported server module, or None if no model package is installed.

    Raises:
        ImportError: if a model package is installed but out of date.
        Exception: whatever else server.py raises at import time.
    """
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "google_cloud_run_services"))
    # server.py reads DB_HOST (google_cloud_run_services/server.py) and defaults it
    # to the production Cloud SQL socket. Point it somewhere inert so an unexpected
    # connection attempt during a unit test fails instead of reaching production.
    os.environ.setdefault("DB_HOST", "unused-for-unit-tests")
    os.environ.setdefault("GENOME_VERSION", "38")
    out_of_date = []
    for tool in ("pangolin", "spliceai"):
        os.environ["TOOL"] = tool
        try:
            return importlib.import_module("server")
        except ModuleNotFoundError as e:
            if (e.name or "").split(".")[0] != tool:
                raise
            print(f"The {tool} package is not installed, so server.py cannot be imported with TOOL={tool}.")
            sys.modules.pop("server", None)
        except ImportError as e:
            # Installed, but doesn't provide what server.py imports. Only worth
            # reporting if no other tool works, so hold it until the loop ends.
            out_of_date.append(f"TOOL={tool}: {e}")
            sys.modules.pop("server", None)

    if out_of_date:
        raise ImportError(
            "server.py could not be imported because an installed model package is out of date "
            f"({'; '.join(out_of_date)}). Update it, or uninstall it to skip these tests.")
    return None


try:
    server = _import_server()
    server_import_error = None
except Exception as e:
    # Held rather than raised here: the live API tests below don't need server.py,
    # and a module-scope raise would take them down too.
    server, server_import_error = None, e

_SKIP_REASON = (
    "server.py could not be imported -- install the spliceai or pangolin package "
    "(see google_cloud_run_services/docker/*/requirements.txt) to run these tests"
)


class ServerHelperTestCase(unittest.TestCase):
    """Base for tests that call server.py helpers directly.

    Skips when no model package is installed, but fails loudly when server.py
    exists and could not be imported -- that's a regression, not a missing
    optional dependency.
    """

    @classmethod
    def setUpClass(cls):
        if server_import_error is not None:
            raise server_import_error
        if server is None:
            raise unittest.SkipTest(_SKIP_REASON)


def query_position_api(tool, hg, position, **extra_params):
    """Query a deployed endpoint with a bare chrom-pos position.

    Args:
        tool: "spliceai" or "pangolin".
        hg: "37" or "38".
        position: chrom-pos string, e.g. "8-140300616".
        **extra_params: additional query-string params to send.

    Returns:
        Parsed JSON response dict.
    """
    return requests.get(
        f"{API_URLS[(tool, hg)]}/{tool}/",
        params={"hg": hg, "distance": "500", "variant": position, "raw": position, **extra_params},
        timeout=180,
    ).json()


class TestPositionParsing(ServerHelperTestCase):
    """parse_position and its precedence relative to parse_variant."""

    def test_accepts_the_separators_the_frontend_emits(self):
        """Every separator index.html's POSITION_ONLY_RE accepts must parse server-side."""
        for position in ("chr8-140300616", "8-140300616", "chr8:140300616", "8 140300616"):
            with self.subTest(position=position):
                self.assertEqual(server.parse_position(position), ("8", 140300616))

    def test_strips_chr_prefix_on_non_numeric_chromosomes(self):
        self.assertEqual(server.parse_position("chrX-1000"), ("X", 1000))

    def test_rejects_thousands_separators(self):
        """index.html strips commas before calling the API; the server does not accept them.

        If that frontend normalization is ever dropped, "chr8:140,300,616" would
        reach the API and 400 instead of returning scores.
        """
        with self.assertRaises(ValueError):
            server.parse_position("chr8-140,300,616")

    def test_rejects_a_bare_chromosome(self):
        with self.assertRaises(ValueError):
            server.parse_position("8")

    def test_rejects_a_full_variant(self):
        with self.assertRaises(ValueError):
            server.parse_position(TEST_VARIANT)

    def test_a_full_variant_parses_as_a_variant(self):
        """The endpoint tries parse_variant first, so a well-formed variant is never REF-only."""
        self.assertEqual(server.parse_variant(TEST_VARIANT), ("8", 140300616, "T", "G"))


class TestPositionOnlyCacheKey(ServerHelperTestCase):
    """get_splicing_scores_cache_key's is_position_only handling.

    The flag exists so that bumping SAI10K_VERSION does not needlessly discard
    cached REF-only responses (SAI-10k-calc requires an ALT allele and never
    runs for them). Getting the key wrong silently invalidates cache entries --
    no error, just every request recomputing on the model.
    """

    def test_delta_score_spliceai_key_includes_the_sai10k_version(self):
        self.assertTrue(
            server.get_splicing_scores_cache_key(
                "spliceai", TEST_VARIANT, "38", "500", "0", "basic", False,
            ).endswith(f"__sai10k-{server.SAI10K_VERSION}")
        )

    def test_position_only_spliceai_key_omits_the_sai10k_version(self):
        self.assertNotIn(
            "sai10k",
            server.get_splicing_scores_cache_key("spliceai", TEST_POSITION, "38", "500", "0", "basic", True),
        )

    def test_pangolin_key_omits_the_sai10k_version_in_both_modes(self):
        for is_position_only in (False, True):
            with self.subTest(is_position_only=is_position_only):
                self.assertNotIn(
                    "sai10k",
                    server.get_splicing_scores_cache_key(
                        "pangolin", TEST_POSITION, "38", "500", "0", "basic", is_position_only,
                    ),
                )

    def test_bumping_the_sai10k_version_spares_position_only_entries(self):
        """A SAI10K_VERSION bump must invalidate delta-score entries but not REF-only ones."""
        delta_key = server.get_splicing_scores_cache_key("spliceai", TEST_VARIANT, "38", "500", "0", "basic", False)
        position_only_key = server.get_splicing_scores_cache_key("spliceai", TEST_POSITION, "38", "500", "0", "basic", True)

        with mock.patch.object(server, "SAI10K_VERSION", "v-bumped"):
            self.assertNotEqual(
                delta_key,
                server.get_splicing_scores_cache_key("spliceai", TEST_VARIANT, "38", "500", "0", "basic", False),
                "delta-score cache entries should be invalidated by a SAI10K_VERSION bump",
            )
            self.assertEqual(
                position_only_key,
                server.get_splicing_scores_cache_key("spliceai", TEST_POSITION, "38", "500", "0", "basic", True),
                "REF-only cache entries never contain SAI-10k output and should survive the bump",
            )

    def test_default_argument_produces_the_documented_key_format(self):
        """Pin the exact key format so it can only ever change deliberately.

        Every component is part of the cache identity, so changing any of them orphans all
        existing entries at once. That is intended when CACHE_VERSION or SAI10K_VERSION is
        bumped, and a bug in every other case.
        """
        self.assertEqual(
            server.get_splicing_scores_cache_key("spliceai", TEST_VARIANT, "38", "500", "0", "basic"),
            f"spliceai__{TEST_VARIANT}__hg38__d500__m0__basic__{server.CACHE_VERSION}__sai10k-{server.SAI10K_VERSION}",
        )

    def test_cache_version_appears_in_every_key(self):
        """Bumping CACHE_VERSION must invalidate pangolin and REF-only entries too.

        Those keys carry no SAI-10k suffix, so SAI10K_VERSION alone cannot invalidate them --
        which is the whole reason CACHE_VERSION exists as a separate component.
        """
        for tool_name in ("spliceai", "pangolin"):
            for is_position_only in (False, True):
                self.assertIn(
                    f"__{server.CACHE_VERSION}",
                    server.get_splicing_scores_cache_key(
                        tool_name, TEST_VARIANT, "38", "500", "0", "basic", is_position_only),
                    f"CACHE_VERSION missing from the {tool_name} key (is_position_only={is_position_only})",
                )

    def test_the_flag_changes_the_key_for_an_otherwise_identical_query(self):
        """Hold every other argument fixed, so only is_position_only can explain the difference.

        Comparing a position key against a full-variant key proves nothing: the
        variant string is interpolated into the key, so those two differ whatever
        the flag does.
        """
        self.assertNotEqual(
            server.get_splicing_scores_cache_key("spliceai", TEST_POSITION, "38", "500", "0", "basic", True),
            server.get_splicing_scores_cache_key("spliceai", TEST_POSITION, "38", "500", "0", "basic", False),
        )


class TestPositionOnlyCacheReadWrite(ServerHelperTestCase):
    """The cache read and write paths must derive the same key for the same request.

    If only one of them forwards is_position_only, REF-only responses are written
    under one key and looked up under another, so they never hit the cache.
    """

    def _keys_used_by_read_and_write(self, is_position_only):
        """Return the cache keys that the read and write helpers pass to run_sql."""
        keys = []

        def fake_run_sql(conn, sql, params):
            keys.append(params[0])
            return []

        with mock.patch.object(server, "run_sql", fake_run_sql):
            server.get_splicing_scores_from_cache(
                None, "spliceai", TEST_POSITION, "38", "500", "0", "basic", is_position_only)
            server.add_splicing_scores_to_cache(
                None, "spliceai", TEST_POSITION, "38", "500", "0", "basic", {"scores": []}, is_position_only)

        # Both helpers swallow exceptions, so a stub that never ran shows up as a
        # short list rather than as a failure inside the call.
        self.assertEqual(len(keys), 2, f"expected one run_sql call from each helper, got {keys}")
        return keys

    def test_read_and_write_use_the_position_only_key(self):
        read_key, write_key = self._keys_used_by_read_and_write(True)

        # Comparing the two keys only to each other would pass even if BOTH helpers
        # stopped forwarding the flag, since they would then agree on the same wrong
        # (delta-mode) key. Pin them to the position-only key instead.
        self.assertEqual(read_key, server.get_splicing_scores_cache_key(
            "spliceai", TEST_POSITION, "38", "500", "0", "basic", True))
        self.assertEqual(write_key, read_key)
        self.assertNotIn("sai10k", read_key)

    def test_read_and_write_use_the_delta_score_key(self):
        read_key, write_key = self._keys_used_by_read_and_write(False)

        self.assertEqual(read_key, server.get_splicing_scores_cache_key(
            "spliceai", TEST_POSITION, "38", "500", "0", "basic", False))
        self.assertEqual(write_key, read_key)
        self.assertIn("sai10k", read_key)


class TestPositionOnlyAPI(unittest.TestCase):
    """REF-only responses from the deployed spliceai/pangolin endpoints."""

    def _check_position_only_response(self, tool, score_keys, ref_score_keys):
        """Query the tool with a bare position and assert the REF-only response shape.

        Args:
            tool: "spliceai" or "pangolin".
            score_keys: keys every returned score dict must contain.
            ref_score_keys: the subset of score_keys holding actual REF probabilities,
                at least one of which must be non-trivial at this locus.
        """
        response = query_position_api(tool, "38", TEST_POSITION)

        self.assertNotIn("error", response, f"{tool} API error: {response.get('error')}")
        self.assertTrue(response.get("isPositionOnly"),
                        f"{tool} response should be flagged isPositionOnly, got: {response.get('isPositionOnly')}")
        self.assertEqual(response.get("chrom"), "8")
        self.assertEqual(response.get("pos"), 140300616)

        # An ALT allele is what makes a delta score meaningful, and there isn't
        # one here -- the response must not claim otherwise.
        for absent_key in ("ref", "alt", "mask", "variant_consequence"):
            self.assertNotIn(absent_key, response, f"{tool} REF-only response should not include '{absent_key}'")

        self.assertTrue(response.get("scores"), f"{tool} returned no scores for {TEST_POSITION}")
        for key in score_keys:
            self.assertIn(key, response["scores"][0],
                          f"{tool} REF-only score is missing '{key}': {response['scores'][0]}")

        # This locus was picked because both tools predict a strong reference splice
        # site here (~0.8-0.9), so an all-zero response means the REF-only path is
        # returning structurally valid but empty predictions. The floor is a sanity
        # bound rather than a snapshot, so model updates don't churn it.
        self.assertGreater(
            max(float(response["scores"][0][key]) for key in ref_score_keys), 0.1,
            f"{tool} returned no meaningful REF score at {TEST_POSITION}: {response['scores'][0]}")

        # The IGV REF track reads these; index.html renders an empty track without
        # them. assertTrue rather than assertIsNotNone, so [] and "" fail too.
        for key in ("allNonZeroScores", "allNonZeroScoresStrand", "allNonZeroScoresTranscriptId"):
            self.assertTrue(response.get(key), f"{tool} REF-only response has an empty '{key}'")

    def test_spliceai_returns_reference_acceptor_and_donor_scores(self):
        self._check_position_only_response(
            "spliceai", ("RA_MAX", "RA_MAX_POS", "RD_MAX", "RD_MAX_POS"), ("RA_MAX", "RD_MAX"))

    def test_pangolin_returns_a_reference_splice_site_score(self):
        self._check_position_only_response("pangolin", ("S_REF", "DP_S"), ("S_REF",))

    def test_a_client_supplied_mask_is_not_echoed_back(self):
        """mask has no meaning without an ALT, so echoing it would imply it was applied."""
        for mask in ("0", "1"):
            with self.subTest(mask=mask):
                response = query_position_api("spliceai", "38", TEST_POSITION, mask=mask)
                # An error response also has no "mask" key, so check the request
                # actually took the REF-only path before reading anything into that.
                self.assertNotIn("error", response, f"API error: {response.get('error')}")
                self.assertTrue(response.get("isPositionOnly"))
                self.assertNotIn("mask", response)

    def test_unparseable_input_returns_an_error(self):
        self.assertIn("error", query_position_api("spliceai", "38", "not-a-position"))

    def test_a_full_variant_still_returns_delta_scores(self):
        """Trying parse_variant first must keep well-formed variants on the delta-score path."""
        response = query_position_api("spliceai", "38", TEST_VARIANT, mask="0")

        self.assertNotIn("error", response, f"API error: {response.get('error')}")
        self.assertNotIn("isPositionOnly", response)
        self.assertIn("DS_AG", response["scores"][0])


if __name__ == "__main__":
    unittest.main()
