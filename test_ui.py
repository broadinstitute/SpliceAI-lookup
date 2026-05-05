"""
Playwright UI integration tests for SpliceAI Lookup.

Tests the web UI end-to-end against a configurable base URL.
Requires: pip3 install playwright && python3 -m playwright install chromium

To run:
    python3 -m unittest test_ui -v

To run against a different URL:
    SPLICEAI_LOOKUP_URL=https://spliceailookup.broadinstitute.org python3 -m unittest test_ui -v
"""

import os
import unittest

from playwright.sync_api import sync_playwright, expect


BASE_URL = os.environ.get(
    "SPLICEAI_LOOKUP_URL",
    "https://broadinstitute.github.io/SpliceAI-lookup-dev/index.html",
)

# Set SPLICEAI_API_ENV=dev to route the frontend's spliceai/pangolin API calls
# to the 'dev'-tagged Cloud Run revisions deployed by `build_and_deploy.py --dev`.
# The frontend (index.html) checks the URL hash for apiEnv=<env> at load time;
# this helper appends it to every page.goto() URL.
API_ENV = os.environ.get("SPLICEAI_API_ENV", "")


def _build_url(hash_str=""):
    """Return BASE_URL with apiEnv=<env> merged into the hash when API_ENV is set."""
    if not API_ENV:
        return BASE_URL + hash_str
    if hash_str.startswith("#"):
        return f"{BASE_URL}{hash_str}&apiEnv={API_ENV}"
    return f"{BASE_URL}#apiEnv={API_ENV}"


# Max time (ms) to wait for API responses. HGVS variants need VEP
# normalization which can take 15s+, and SpliceAI model inference can
# take 10-30s for uncached variants.
TIMEOUT_MS = 60_000


class TestSpliceAILookupUI(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._pw = sync_playwright().start()
        cls.browser = cls._pw.chromium.launch(headless=True)

    @classmethod
    def tearDownClass(cls):
        cls.browser.close()
        cls._pw.stop()

    def setUp(self):
        # Use a fresh browser context per test (not just a new page on the
        # shared default context). Pages share localStorage / sessionStorage /
        # IndexedDB / cookies within a context, so test order would start to
        # matter as soon as the app — or any third-party script it loads —
        # begins using any of those stores.
        self.context = self.browser.new_context()
        self.page = self.context.new_page()
        self.page.set_default_timeout(TIMEOUT_MS)
        self.page.goto(_build_url(), wait_until="domcontentloaded")

    def tearDown(self):
        self.context.close()

    def _get_url_hash(self):
        """Return the URL fragment (everything after #), or empty string."""
        return self.page.url.split("#", 1)[-1] if "#" in self.page.url else ""

    def _wait_for_submit(self):
        """Wait for the submit button to enter and then exit its loading state."""
        self.page.wait_for_selector("#submit-button.loading", timeout=5_000)
        self.page.wait_for_selector("#submit-button:not(.loading)", timeout=TIMEOUT_MS)

    def _click_label(self, input_selector):
        """Click the Semantic UI label for a hidden checkbox/radio input.

        Semantic UI wraps inputs in <div class="ui checkbox"> and hides the
        real <input>.  Clicking the <label> sibling triggers the JS handler.
        """
        self.page.locator(input_selector).locator("xpath=..").locator("label").click()

    def _submit_variant(self, variant, hg="38"):
        """Enter a variant, select genome version, click submit, and wait for completion."""
        if hg == "37":
            self._click_label("input[name='hg'][value='37']")

        self.page.fill("#search-box", variant)
        self.page.click("#submit-button")
        self._wait_for_submit()

    # ------------------------------------------------------------------
    # Page load
    # ------------------------------------------------------------------

    def test_page_loads(self):
        """Title, search box, submit button, and genome version radios are present."""
        self.assertIn("SpliceAI", self.page.title())
        self.assertTrue(self.page.query_selector("#search-box"))
        self.assertTrue(self.page.query_selector("#submit-button"))
        self.assertTrue(self.page.query_selector("input[name='hg'][value='37']"))
        self.assertTrue(self.page.query_selector("input[name='hg'][value='38']"))

    # ------------------------------------------------------------------
    # Variant submission and result tables
    # ------------------------------------------------------------------

    def test_variant_shows_spliceai_results(self):
        """Submitting a coordinate variant shows the SpliceAI results table."""
        self._submit_variant("8-140300616-T-G")

        self.assertTrue(self.page.is_visible("#response-box"))
        self.assertTrue(self.page.is_visible("#spliceai-table"))

        rows = self.page.query_selector_all(".spliceai-result-row")
        self.assertGreater(len(rows), 0, "Expected at least one SpliceAI result row")

    def test_variant_shows_pangolin_results(self):
        """Submitting a coordinate variant shows the Pangolin results table."""
        self._submit_variant("8-140300616-T-G")

        self.assertTrue(self.page.is_visible("#pangolin-table"))

        rows = self.page.query_selector_all(".pangolin-result-row")
        self.assertGreater(len(rows), 0, "Expected at least one Pangolin result row")

    def test_spliceai_score_labels(self):
        """SpliceAI result rows contain the four expected score type labels."""
        self._submit_variant("8-140300616-T-G")

        table_text = self.page.inner_text("#spliceai-table")
        for label in ("Acceptor Loss", "Donor Loss", "Acceptor Gain", "Donor Gain"):
            self.assertIn(label, table_text, f"Missing label '{label}' in SpliceAI table")

    def test_pangolin_score_labels(self):
        """Pangolin result rows contain the two expected score type labels."""
        self._submit_variant("8-140300616-T-G")

        table_text = self.page.inner_text("#pangolin-table")
        for label in ("Splice Loss", "Splice Gain"):
            self.assertIn(label, table_text, f"Missing label '{label}' in Pangolin table")

    def test_insertion_variant(self):
        """Submitting an insertion variant shows results."""
        self._submit_variant("1-1042601-A-AGAGAG")

        self.assertTrue(self.page.is_visible("#response-box"))
        self.assertTrue(self.page.is_visible("#spliceai-table"))

        rows = self.page.query_selector_all(".spliceai-result-row")
        self.assertGreater(len(rows), 0, "Expected SpliceAI results for insertion variant")

    def test_deletion_variant(self):
        """Submitting a deletion variant shows results."""
        self._submit_variant("1-1042466-GGGC-G")

        self.assertTrue(self.page.is_visible("#response-box"))
        self.assertTrue(self.page.is_visible("#spliceai-table"))

        rows = self.page.query_selector_all(".spliceai-result-row")
        self.assertGreater(len(rows), 0, "Expected SpliceAI results for deletion variant")

    # ------------------------------------------------------------------
    # URL hash
    # ------------------------------------------------------------------

    def test_url_hash_updates_after_submit(self):
        """URL hash contains all expected params after submit."""
        self._submit_variant("8-140300616-T-G")

        url_hash = self._get_url_hash()
        for param in ("variant=", "hg=38", "bc=", "distance=", "mask=", "ra="):
            self.assertIn(param, url_hash, f"URL hash missing '{param}': {url_hash}")

    def test_url_hash_navigation(self):
        """Loading the page with a pre-populated hash auto-submits and shows results."""
        hash_url = _build_url("#variant=8-140300616-T-G&hg=38&bc=basic&distance=500&mask=0&ra=0")
        # Navigate away first so the next goto triggers a full page load
        # (same-origin hash-only changes don't re-fire $(document).ready).
        self.page.goto("about:blank")
        self.page.goto(hash_url, wait_until="domcontentloaded")

        # Use a longer timeout than the default. This is the one test in the
        # suite that cannot benefit from any earlier same-context warmup: it
        # navigates from about:blank, which forces the full GitHub-Pages load
        # plus a fresh /spliceai/ call against the configured Cloud Run
        # backend. A cold Cloud Run instance (~10-15s) plus VEP normalization
        # plus model inference can plausibly exceed the 60s default, and we
        # don't want a 60s flake to mask a real regression.
        self.page.wait_for_selector(".spliceai-result-row", timeout=TIMEOUT_MS * 2)
        self.assertTrue(self.page.is_visible("#spliceai-table"))

    def test_url_hash_restores_form_state(self):
        """Loading with hash params restores form state (hg, bc, distance, mask, ra).

        Uses `nosubmit=1` to skip the auto-submit so this test verifies state
        restoration only — independent of backend latency or availability.
        """
        hash_url = _build_url("#variant=8-140300616-T-G&hg=37&bc=comprehensive&distance=1000&mask=1&ra=1&nosubmit=1")
        self.page.goto("about:blank")
        self.page.goto(hash_url, wait_until="domcontentloaded")

        # Wait for the JS to apply hash params to form elements (happens in
        # applyUrlSettingsToFormElements; with nosubmit=1, no auto-submit fires).
        self.page.wait_for_function(
            "() => document.querySelector('#max-distance-input').value === '1000'",
            timeout=10_000,
        )

        self.assertTrue(self.page.is_checked("input[name='hg'][value='37']"))
        self.assertTrue(self.page.is_checked("input[name='gencode-gene-set'][value='comprehensive']"))
        self.assertEqual(self.page.input_value("#max-distance-input"), "1000")
        self.assertTrue(self.page.is_checked("input[name='mask']"))
        self.assertTrue(self.page.is_checked("input[name='show-ref-alt']"))

    # ------------------------------------------------------------------
    # Genome version
    # ------------------------------------------------------------------

    def test_genome_version_hg37(self):
        """Selecting hg37 and submitting uses hg=37 in the URL hash and shows results."""
        # Note: 8-140300616-T-G is an hg38 coordinate; we're testing the
        # hg37 UI workflow, not the biological accuracy of this position.
        self._submit_variant("8-140300616-T-G", hg="37")

        self.assertIn("hg=37", self._get_url_hash())
        self.assertTrue(self.page.is_visible("#response-box"))

    # ------------------------------------------------------------------
    # Form options
    # ------------------------------------------------------------------

    def test_max_distance_input(self):
        """Changing max distance to a non-default value is reflected in the URL hash."""
        self.page.fill("#max-distance-input", "1000")
        self._submit_variant("8-140300616-T-G")

        self.assertIn("distance=1000", self._get_url_hash())

    def test_mask_option(self):
        """Enabling the mask checkbox includes mask=1 in the URL hash."""
        self._click_label("input[name='mask']")
        self._submit_variant("8-140300616-T-G")

        self.assertIn("mask=1", self._get_url_hash())

    def test_gencode_comprehensive_option(self):
        """Selecting 'comprehensive' gencode set includes bc=comprehensive in the URL hash."""
        self._click_label("input[name='gencode-gene-set'][value='comprehensive']")
        self._submit_variant("8-140300616-T-G")

        self.assertIn("bc=comprehensive", self._get_url_hash())

    def test_ref_alt_score_columns(self):
        """Enabling REF/ALT scores checkbox shows the REF and ALT score columns."""
        self._click_label("input[name='show-ref-alt']")
        self._submit_variant("8-140300616-T-G")

        # The result table includes .ref-score-column / .alt-score-column <td>
        # elements only when the checkbox is checked at submit time.  The
        # header <th> elements always exist but may be hidden by default CSS.
        # Check that at least one result-row <td> with these classes exists.
        ref_tds = self.page.query_selector_all("#spliceai-table td.ref-score-column")
        alt_tds = self.page.query_selector_all("#spliceai-table td.alt-score-column")
        self.assertGreater(len(ref_tds), 0, "Expected REF score <td> elements in results")
        self.assertGreater(len(alt_tds), 0, "Expected ALT score <td> elements in results")

    # ------------------------------------------------------------------
    # Show more examples
    # ------------------------------------------------------------------

    def test_show_more_examples(self):
        """Clicking [show more examples] reveals hidden example rows."""
        hidden_rows = self.page.query_selector_all(".more-details2")
        self.assertGreater(len(hidden_rows), 0, "Expected .more-details2 rows to exist")
        for row in hidden_rows:
            self.assertFalse(row.is_visible())

        self.page.click("#more-details2-button")

        for row in hidden_rows:
            self.assertTrue(row.is_visible())

    # ------------------------------------------------------------------
    # Error handling
    # ------------------------------------------------------------------

    def test_invalid_variant_shows_error(self):
        """An unrecognizable variant string shows an error message."""
        self._submit_variant("INVALID-XYZ-VARIANT")

        self.assertTrue(self.page.is_visible("#error-box"))
        error_text = self.page.inner_text("#error-box")
        self.assertTrue(len(error_text.strip()) > 0, "Error box should contain a message")

    # ------------------------------------------------------------------
    # Transcript toggle
    # ------------------------------------------------------------------

    def test_transcript_toggle(self):
        """Main/all transcript buttons toggle visibility of non-main transcript rows."""
        self._submit_variant("1-930130-C-G")

        non_main_rows = self.page.query_selector_all(".spliceai-result-row.non-main-transcript")
        if not non_main_rows:
            self.skipTest("Variant 1-930130-C-G returned only main transcripts; cannot test toggle")

        # Click "All Transcripts" -- non-main rows should be visible
        self.page.click("#all-transcript-button")
        expect(self.page.locator(".spliceai-result-row.non-main-transcript").first).to_be_visible()

        # Click "Main Transcript" -- non-main rows should be hidden
        self.page.click("#main-transcript-button")
        expect(self.page.locator(".spliceai-result-row.non-main-transcript").first).to_be_hidden()

    # ------------------------------------------------------------------
    # Other predictors
    # ------------------------------------------------------------------

    def test_other_predictors_content(self):
        """The other-predictors table appears with predictor names and scores."""
        self._submit_variant("8-140300616-T-G")

        self.assertTrue(self.page.is_visible("#other-predictors-table"))
        table_text = self.page.inner_text("#other-predictors-table")
        # At least some standard predictors should appear
        found = [p for p in ("CADD", "REVEL", "AlphaMissense", "PrimateAI", "PhyloP")
                 if p.lower() in table_text.lower()]
        self.assertGreater(len(found), 0,
                           f"Expected predictor names in other-predictors table, got: {table_text[:300]}")

    # ------------------------------------------------------------------
    # Score color-coding
    # ------------------------------------------------------------------

    def test_score_color_coding(self):
        """Score cells have color-coded backgrounds based on score thresholds."""
        # Use a BRCA1 splice variant known to produce high delta scores
        self._submit_variant("17-41199659-C-T", hg="37")

        # Collect background colors from SpliceAI score cells. The score
        # column is the 4th <td> in each result row (after variant, gene,
        # delta-type).  getScoreStyle() applies inline styles:
        #   >= 0.8  ->  #fccfb8 (red/salmon)
        #   >= 0.5  ->  #fff19d (yellow)
        #   >= 0.2  ->  #cdffd7 (green)
        #   < 0.01  ->  color:#BBBBBB (gray text; no background)
        # Browsers normalise inline backgroundColor to rgb(r, g, b) when read
        # back via el.style.backgroundColor, so compare against the rgb form.
        known_rgb = {
            "rgb(252, 207, 184)",  # #fccfb8
            "rgb(255, 241, 157)",  # #fff19d
            "rgb(205, 255, 215)",  # #cdffd7
        }
        found_colors = set()
        for cell in self.page.query_selector_all("#spliceai-table td"):
            bg = cell.evaluate("el => el.style.backgroundColor")
            if bg:
                found_colors.add(bg)

        # At least one score cell should have a threshold-based background
        self.assertTrue(
            found_colors,
            "Expected at least one score cell with a color-coded background",
        )
        # And every background colour we see must be one of the palette
        # values — a cell with a rogue colour means getScoreStyle drifted.
        unexpected = found_colors - known_rgb
        self.assertFalse(
            unexpected,
            f"Unexpected score-cell background colors outside the known "
            f"palette: {unexpected} (expected a subset of {known_rgb})",
        )

    # ------------------------------------------------------------------
    # Insertion modal dialog
    # ------------------------------------------------------------------

    def test_insertion_modal_dialog(self):
        """Insertion variants can show a detailed per-base score modal."""
        # This insertion is known to trigger the SCORES_FOR_INSERTED_BASES
        # detail view (gain score at position 0).
        self._submit_variant("2-47790924-C-CAGTTG")

        # Look for the modal trigger icon (window maximize outline icon)
        modal_icons = self.page.query_selector_all("#spliceai-table .window.maximize.outline.icon")
        if not modal_icons:
            self.skipTest("No modal icons found for this insertion variant")

        modal_icons[0].click()
        # A Semantic UI modal should appear with a score table inside
        expect(self.page.locator(".ui.modal.visible")).to_be_visible(timeout=5_000)
        modal_text = self.page.inner_text(".ui.modal.visible")
        self.assertIn("REF acceptor score", modal_text)

    # ------------------------------------------------------------------
    # External links in results
    # ------------------------------------------------------------------

    def test_external_links_in_results(self):
        """Result rows contain external links to genome browsers and databases."""
        self._submit_variant("8-140300616-T-G")

        spliceai_html = self.page.inner_html("#spliceai-table")
        for domain in ("genome.ucsc.edu", "gnomad.broadinstitute.org",
                        "omim.org", "ensembl.org", "genecards.org"):
            self.assertIn(domain, spliceai_html, f"Missing link to {domain} in SpliceAI table")

    # ------------------------------------------------------------------
    # Enter key submission
    # ------------------------------------------------------------------

    def test_enter_key_submits_from_search_box(self):
        """Pressing Enter in the search box triggers form submission."""
        self.page.fill("#search-box", "8-140300616-T-G")
        self.page.press("#search-box", "Enter")

        self._wait_for_submit()
        self.assertTrue(self.page.is_visible("#response-box"))

    # ------------------------------------------------------------------
    # HGVS input
    # ------------------------------------------------------------------

    def test_hgvs_variant_input(self):
        """Submitting an HGVS variant (with VEP normalization) shows results for the correct gene."""
        self._submit_variant("NM_001089.3:c.875A>T")

        self.assertTrue(self.page.is_visible("#spliceai-table"))
        table_text = self.page.inner_text("#spliceai-table")
        self.assertIn("ABCA3", table_text, "Expected gene name ABCA3 in SpliceAI results")

    def test_hgvs_consequence_display(self):
        """HGVS variant shows a VEP consequence link in the results."""
        self._submit_variant("NM_001089.3:c.875A>T")

        spliceai_html = self.page.inner_html("#spliceai-table")
        # The consequence should link to Ensembl's predicted_data page
        self.assertIn("ensembl.org/info/genome/variation/prediction", spliceai_html,
                      "Expected VEP consequence link in results")

    def test_normalized_variant_display(self):
        """HGVS input shows a normalized coordinate with an arrow indicator."""
        self._submit_variant("NM_001089.3:c.875A>T")

        table_text = self.page.inner_text("#spliceai-table")
        # The normalization arrow is shown when the input differs from the
        # normalized chrom-pos-ref-alt form.
        self.assertIn("\u21d2", table_text,
                      "Expected normalization arrow (\u21d2) for HGVS input")

    # ------------------------------------------------------------------
    # Partial API failure
    # ------------------------------------------------------------------

    def test_partial_api_failure(self):
        """When SpliceAI returns an error, the error box shows non-empty content.

        Uses a variant deep in an unannotated intergenic region that should not
        produce SpliceAI scores. If a future change makes this variant succeed
        (e.g. annotation expansion), pick a different variant — the assertion
        below is intentionally strict about requiring a non-trivial error
        message so this test cannot silently degrade into a duplicate of
        test_variant_shows_spliceai_results.
        """
        self._submit_variant("1-10000-A-G")

        # Wait for loading to finish. The page must not still be showing the
        # spinner — an infinite-loading state should fail this test, not pass
        # it silently.
        self.page.wait_for_function(
            "!document.querySelector('#submit-button').classList.contains('loading')",
            timeout=30000,
        )

        # We expect the error path here. Require a non-trivial error message so
        # that "Network error" or a stray "[object Object]" doesn't silently
        # satisfy the assertion. Specific wording is verified by the
        # server-side tests; we just check that *some* meaningful error text
        # was rendered. Avoid substring-matching against a hardcoded list of
        # error fragments — that list inevitably rots whenever a server-side
        # message is rephrased and produces test failures unrelated to the UI.
        self.assertTrue(self.page.is_visible("#error-box"),
                        "Expected the error box to be visible for an unannotated variant")
        error_text = (self.page.inner_text("#error-box") or "").strip()
        self.assertGreaterEqual(
            len(error_text), 20,
            f"Error box should contain a substantive message; got: {error_text!r}",
        )

    # ------------------------------------------------------------------
    # Loading state
    # ------------------------------------------------------------------

    def test_loading_state_during_request(self):
        """Submit button gets 'loading' class during the API request."""
        self.page.fill("#search-box", "8-140300616-T-G")
        self.page.click("#submit-button")

        # wait_for_selector confirms the class is present; no extra evaluate needed.
        self.page.wait_for_selector("#submit-button.loading", timeout=5_000)

        # After the request completes, loading class should be removed.
        self.page.wait_for_selector("#submit-button:not(.loading)", timeout=TIMEOUT_MS)


if __name__ == "__main__":
    unittest.main()
