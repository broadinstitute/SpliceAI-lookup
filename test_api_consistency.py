"""
API-level score consistency tests for SpliceAI Lookup.

Queries the production SpliceAI and Pangolin APIs for test variants and
compares delta scores against a snapshot captured from the current
production deployment.

To capture the baseline (run ONCE before deploying changes):
    python3 test_api_consistency.py --capture

To run the consistency tests:
    python3 -m unittest test_api_consistency -v

The --capture step writes expected_scores.json. The test step loads it
and compares each variant's scores against the snapshot.

Note on the first run after a model re-pin: the pinned model commit is part of the
response cache key (MODEL_COMMIT in server.py), so re-pinning retires every cached
entry and each of these queries recomputes instead of being served from cache. A
computed request counts against the per-IP rate limit and a cached one does not, so
that first run costs up to 76 computed requests against a budget of 150 per 7 minutes.
It fits, and it also repopulates what it recomputed, so a later identical run is served
from cache and costs nothing again.
"""

import json
import os
import sys
import time
import unittest

import requests



# The GCP-assigned Cloud Run hostname slug changes if the service is recreated
# or moved between projects, so allow it to be overridden via
# SPLICEAI_API_URL_TEMPLATE=https://{tool}-{hg}-your-new-slug-uc.a.run.app
# (placeholders {tool} and {hg} are substituted per-endpoint).
#
# Set SPLICEAI_API_ENV=dev to point at the 'dev'-tagged Cloud Run revisions
# deployed by `build_and_deploy.py --dev` (no traffic). Ignored if
# SPLICEAI_API_URL_TEMPLATE is also set.
_DEFAULT_API_URL_TEMPLATE = "https://{tool}-{hg}-xwkwwwxdwq-uc.a.run.app"
if os.environ.get("SPLICEAI_API_ENV") == "dev":
    _DEFAULT_API_URL_TEMPLATE = "https://dev---{tool}-{hg}-xwkwwwxdwq-uc.a.run.app"
_API_URL_TEMPLATE = os.environ.get("SPLICEAI_API_URL_TEMPLATE", _DEFAULT_API_URL_TEMPLATE)
API_URLS = {
    (tool, hg): _API_URL_TEMPLATE.format(tool=tool, hg=hg)
    for tool in ("spliceai", "pangolin")
    for hg in ("37", "38")
}

TRANSCRIPT_PRIORITY = {"MS": 3, "MP": 2, "C": 1, "N": 0}

# Largest delta-score difference still treated as a match. Scores are compared after rounding
# to 2 decimals, so a model whose own output precision changes (bw2/SpliceAI and bw2/Pangolin
# raised theirs from 2 to 3 decimals in Aug 2026) shifts a value like 0.925 to either 0.92 or
# 0.93 depending on which side of the boundary it lands. That is a reporting difference, not a
# scoring regression, and this tolerance is what keeps it from being reported as one. Anything
# larger is a real change in what the model predicts.
SCORE_TOLERANCE = 0.01

# Substrings of the API error messages that mean "there is genuinely nothing to score at this
# variant". Those are a property of the variant, so --capture records them and the test skips
# that variant from then on. Every other failure -- a rate-limit rejection, a 5xx, a dropped
# connection -- is a property of the moment, and recording one would bake a permanent skip into
# the baseline for a variant that is fine. An earlier capture that tripped the rate limiter wrote
# 20 such entries into the baseline, which is what motivated this split; the baseline was
# re-captured afterwards and now has none.
CAPTURABLE_ERROR_MARKERS = (
    "did not return any scores",
    "unable to compute",
)

# Seconds between capture queries. The API restricts an IP that logs 150 computed requests
# within 7 minutes (block_ips in server.py), and a capture fires one request per
# (variant, hg, tool) with nothing in between, which is how those rate-limit entries got into
# the baseline in the first place.
CAPTURE_DELAY_SECONDS = 3


def scores_match(actual, expected):
    """Whether two rounded delta scores are within SCORE_TOLERANCE of each other.

    The difference is rounded before the comparison because binary floating point cannot hold
    2-decimal values exactly: abs(0.93 - 0.92) evaluates to 0.010000000000000009, which a plain
    `<= 0.01` would reject -- rejecting precisely the off-by-one-hundredth case this exists to
    accept.
    """
    return round(abs(actual - expected), 6) <= SCORE_TOLERANCE

EXPECTED_SCORES_PATH = os.path.join(os.path.dirname(__file__) or ".", "expected_scores.json")

# (variant, hg, label) -- one entry per (coord, assembly) pair.
#
# hg38 coords and their hg19 liftovers are listed as separate entries
# (Broad liftover service, bcftools plugin). Querying hg37 with an hg38 coord
# would miss the target gene by ~1-2 Mb and return no scores, so each pair
# gets its own entry with the assembly-appropriate coordinate.
#
# HGVS variants are resolved per-assembly by Ensembl VEP during --capture,
# so a single HGVS string gets one entry per assembly.
VARIANTS = [
    # -- From index.html examples --
    ("8-140300616-T-G",     "38", "index.html: simple coordinate example"),
    ("8-141310715-T-G",     "37", "index.html: simple coordinate example (hg19 liftover)"),
    ("6-31740453-G-T",      "38", "index.html: space-separated example"),
    ("6-31708230-G-T",      "37", "index.html: space-separated example (hg19 liftover)"),
    ("1-930130-C-G",        "38", "index.html: overlapping genes"),
    ("1-865510-C-G",        "37", "index.html: overlapping genes (hg19 liftover)"),
    ("1-1042601-A-AGAGAG",  "38", "index.html: insertion"),
    ("1-977981-A-AGAGAG",   "37", "index.html: insertion (hg19 liftover)"),
    ("1-1042466-GGGC-G",    "38", "index.html: deletion"),
    ("1-977846-GGGC-G",     "37", "index.html: deletion (hg19 liftover)"),

    # -- From Google Doc figures --
    ("17-58709857-A-C",     "38", "doc fig1: RAD51C exon skipping + partial exon deletion"),
    ("17-56787218-A-C",     "37", "doc fig1: RAD51C (hg19 liftover)"),
    ("22-28695127-C-T",     "38", "doc fig2: CHEK2 exon skipping"),
    ("22-29091115-C-T",     "37", "doc fig2: CHEK2 (hg19 liftover)"),
    ("16-23641111-T-C",     "38", "doc fig3: PALB2 partial exon deletion"),
    ("16-23652432-T-C",     "37", "doc fig3: PALB2 (hg19 liftover)"),
    ("17-58696687-T-A",     "38", "doc fig4: RAD51C partial intron retention"),
    ("17-56774048-T-A",     "37", "doc fig4: RAD51C partial intron retention (hg19 liftover)"),
    ("22-28741768-C-A",     "38", "doc fig5: CHEK2 non-coding partial exon deletion"),
    ("22-29137756-C-A",     "37", "doc fig5: CHEK2 non-coding partial exon deletion (hg19 liftover)"),
    ("13-32345247-T-G",     "38", "doc fig6: BRCA2 pseudoexon"),
    ("13-32919384-T-G",     "37", "doc fig6: BRCA2 pseudoexon (hg19 liftover)"),
    ("13-32363534-G-T",     "38", "doc fig8: BRCA2 exon skipping"),
    ("13-32937671-G-T",     "37", "doc fig8: BRCA2 exon skipping (hg19 liftover)"),

    # -- From changelog / code comments --
    ("2-47790924-C-CAGTTG", "38", "changelog: insertion warning example"),
    ("2-48018063-C-CAGTTG", "37", "changelog: insertion warning example (hg19 liftover)"),
    ("1-55039916-G-A",      "38", "index.html source: example variant comment"),
    ("1-55505589-G-A",      "37", "index.html source: example variant comment (hg19 liftover)"),

    # -- BRCA1 splice variants --
    # hg38 coords are themselves liftovers of hg19-native coords; the hg19
    # originals are in the BRCA1 hg19-native block below.
    ("17-43047638-C-G",     "38", "BRCA1 c.5467+5G>C (hg38, = hg19 liftover of 17-41199655)"),
    ("17-43047642-C-T",     "38", "BRCA1 c.5467+1G>A (hg38, = hg19 liftover of 17-41199659)"),

    # -- BRCA1 region on hg38 (ref allele verified via Ensembl REST API) --
    ("17-43045714-G-C",     "38", "BRCA1 region hg38, ref=G verified"),
    ("17-41197731-G-C",     "37", "BRCA1 region (hg19 liftover)"),

    # -- X chromosome (MECP2 region, clinically important for splicing) --
    ("X-154031414-C-T",     "38", "MECP2 region chrX hg38"),
    ("X-153296865-C-T",     "37", "MECP2 region chrX (hg19 liftover)"),

    # -- HGVS example from index.html, resolved to coordinate during --capture --
    # NM_000249.4:c.116G>A (MLH1). Ensembl VEP resolves per-assembly, so the
    # same HGVS string produces different chrom/pos per hg during --capture.
    ("HGVS:NM_000249.4:c.116G>A", "38", "MLH1 HGVS hg38, resolved during capture"),
    ("HGVS:NM_000249.4:c.116G>A", "37", "MLH1 HGVS hg37, resolved during capture"),

    # -- BRCA1 hg19-native coordinates (hg38 liftovers are in the BRCA1 block above) --
    ("17-41199655-C-G",     "37", "BRCA1 c.5467+5G>C (hg19 native coord)"),
    ("17-41199659-C-T",     "37", "BRCA1 c.5467+1G>A (hg19 native coord)"),
]


def query_api(tool, hg, variant):
    """Query a production API endpoint and return parsed JSON.

    Args:
        tool: "spliceai" or "pangolin"
        hg: "37" or "38"
        variant: chrom-pos-ref-alt format

    Returns:
        Parsed JSON response dict.
    """
    base_url = API_URLS[(tool, hg)]
    url = f"{base_url}/{tool}/?hg={hg}&distance=500&mask=0&variant={variant}&raw={variant}"
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    return response.json()


def normalize_hgvs(hgvs_string, hg="38"):
    """Resolve an HGVS string to a chrom-pos-ref-alt coordinate via the Ensembl VEP API.

    Uses the same Ensembl endpoint as the frontend (index.html normalizeVariant).

    Args:
        hgvs_string: e.g. "NM_000249.4:c.116G>A"
        hg: genome version ("37" or "38")

    Returns:
        Normalized variant string in chrom-pos-ref-alt format.
    """
    prefix = "grch37." if hg == "37" else ""
    url = f"https://{prefix}rest.ensembl.org/vep/human/hgvs/{hgvs_string}?content-type=application/json&vcf_string=1"
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    data = response.json()
    if isinstance(data, dict) and "error" in data:
        raise ValueError(f"VEP API error for {hgvs_string}: {data['error']}")
    if isinstance(data, list) and len(data) > 0:
        vcf_string = data[0].get("vcf_string")
        if vcf_string:
            return vcf_string
    raise ValueError(f"Could not extract variant from VEP response for {hgvs_string}")


def get_top_transcript(scores, tool):
    """Select the highest-priority transcript from a list of score dicts.

    Priority order: MS > MP > C > N. Tie-break by sum of |delta scores|.

    Args:
        scores: list of transcript score dicts from API response
        tool: "spliceai" or "pangolin"

    Returns:
        The highest-priority transcript score dict, or None.
    """
    if not scores:
        return None

    best = None
    best_priority = -1
    best_sum = -1.0

    for transcript in scores:
        priority = TRANSCRIPT_PRIORITY.get(transcript.get("t_priority", "N"), 0)
        if tool == "pangolin":
            score_sum = abs(float(transcript.get("DS_SL", 0))) + abs(float(transcript.get("DS_SG", 0)))
        else:
            score_sum = sum(
                abs(float(transcript.get(k, 0)))
                for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")
            )
        if priority > best_priority or (priority == best_priority and score_sum > best_sum):
            best = transcript
            best_priority = priority
            best_sum = score_sum

    return best


def extract_scores(transcript, tool):
    """Extract the relevant scores from a transcript dict for comparison.

    Args:
        transcript: transcript score dict
        tool: "spliceai" or "pangolin"

    Returns:
        Dict with g_name, t_id, t_priority, and delta scores.

    Raises:
        KeyError: if a required score key is missing from the transcript.
    """
    result = {
        "g_name": transcript.get("g_name", ""),
        "t_id": transcript.get("t_id", ""),
        "t_priority": transcript.get("t_priority", "N"),
    }
    score_keys = ("DS_SL", "DS_SG") if tool == "pangolin" else ("DS_AG", "DS_AL", "DS_DG", "DS_DL")
    for key in score_keys:
        if key not in transcript:
            raise KeyError(f"Missing score key '{key}' in transcript {result['t_id']}")
        result[key] = round(float(transcript[key]), 2)
    return result


def _capture_variant(expected, variant, hg, tool, label):
    """Capture scores for one (variant, hg, tool) combination.

    Raises on any failure that is a property of the moment rather than of the variant, so a
    capture that hits the rate limiter or a transport error fails loudly instead of writing a
    baseline that silently skips those variables forever. See CAPTURABLE_ERROR_MARKERS.
    """
    key = f"{variant}|{hg}|{tool}"
    actual_variant = variant
    if variant.startswith("HGVS:"):
        hgvs_string = variant[len("HGVS:"):]
        print(f"  Normalizing {hgvs_string} for hg{hg} ...")
        actual_variant = normalize_hgvs(hgvs_string, hg)
        print(f"    -> {actual_variant}")

    print(f"  Querying {tool} hg{hg} for {actual_variant} ({label}) ...")
    response = query_api(tool, hg, actual_variant)
    if response.get("error"):
        if not any(marker in response["error"].lower() for marker in CAPTURABLE_ERROR_MARKERS):
            raise RuntimeError(
                f"{key}: the API failed for a reason that is not a property of this variant, so "
                f"the capture is aborting rather than recording it: {response['error']}")
        print(f"    API error: {response['error']}")
        expected[key] = {"error": response["error"]}
        return

    top = get_top_transcript(response.get("scores", []), tool)
    if not top:
        print(f"    No transcripts returned")
        expected[key] = {"error": "no transcripts"}
        return

    scores = extract_scores(top, tool)
    if variant.startswith("HGVS:"):
        scores["resolved_variant"] = actual_variant
    expected[key] = scores
    print(f"    {scores}")


def capture_expected_scores():
    """Query all variants from production and write expected_scores.json."""
    expected = {}
    for variant, hg, label in VARIANTS:
        for tool in ("spliceai", "pangolin"):
            _capture_variant(expected, variant, hg, tool, label)
            time.sleep(CAPTURE_DELAY_SECONDS)

    with open(EXPECTED_SCORES_PATH, "w") as f:
        json.dump(expected, f, indent=2, sort_keys=True)

    errors = sum(1 for v in expected.values() if "error" in v)
    print(f"\nWrote {len(expected)} entries to {EXPECTED_SCORES_PATH} "
          f"({len(expected) - errors} succeeded, {errors} errors)")


def load_expected_scores():
    """Load expected_scores.json, returning the parsed dict."""
    if not os.path.exists(EXPECTED_SCORES_PATH):
        raise FileNotFoundError(
            f"{EXPECTED_SCORES_PATH} not found. Run 'python3 test_api_consistency.py --capture' first."
        )
    with open(EXPECTED_SCORES_PATH) as f:
        return json.load(f)


class TestAPIConsistency(unittest.TestCase):
    """Verify production API scores match the expected snapshot."""

    @classmethod
    def setUpClass(cls):
        cls.expected = load_expected_scores()

    def _check_variant(self, variant, hg, tool):
        """Query the API and compare scores against the expected snapshot."""
        key = f"{variant}|{hg}|{tool}"
        if key not in self.expected:
            # Fail rather than skip. A key missing from the baseline means the two have drifted
            # apart -- VARIANTS gained an entry, or the baseline was captured from a different
            # list -- and skipping made that invisible: the baseline once covered only 44 of the
            # 76 generated cases while the run still reported success.
            self.fail(f"No expected scores for {key}. The baseline no longer covers VARIANTS; "
                      f"re-capture it with 'python3 test_api_consistency.py --capture'.")

        expected = self.expected[key]
        if "error" in expected:
            self.skipTest(f"Expected error for {key}: {expected['error']}")

        # For HGVS variants, use the resolved coordinate
        actual_variant = expected.get("resolved_variant", variant)
        if actual_variant.startswith("HGVS:"):
            self.skipTest(f"HGVS variant {variant} was not resolved during capture")

        response = query_api(tool, hg, actual_variant)
        self.assertNotIn("error", response, f"API error for {actual_variant}: {response.get('error')}")

        top = get_top_transcript(response.get("scores", []), tool)
        self.assertIsNotNone(top, f"No transcripts returned for {actual_variant}")

        actual = extract_scores(top, tool)

        # Compare transcript identity (base ID without version suffix)
        self.assertEqual(
            actual["t_id"].split(".")[0],
            expected["t_id"].split(".")[0],
            f"{key}: transcript ID mismatch",
        )
        self.assertEqual(actual["g_name"], expected["g_name"], f"{key}: gene name mismatch")
        self.assertEqual(actual["t_priority"], expected["t_priority"], f"{key}: priority mismatch")

        # Compare delta scores at 2 decimal places, within SCORE_TOLERANCE
        score_keys = ("DS_SL", "DS_SG") if tool == "pangolin" else ("DS_AG", "DS_AL", "DS_DG", "DS_DL")
        for sk in score_keys:
            self.assertTrue(
                scores_match(actual[sk], expected[sk]),
                f"{key}: {sk} mismatch: got {actual[sk]}, expected {expected[sk]} "
                f"(differs by {abs(actual[sk] - expected[sk]):.3f}, tolerance {SCORE_TOLERANCE})",
            )


# Dynamically generate one test method per (variant, hg, tool) combination.
def _make_test(variant, hg, tool):
    def test_fn(self):
        self._check_variant(variant, hg, tool)
    safe_name = variant.replace("-", "_").replace(":", "_").replace(".", "_")
    test_fn.__name__ = f"test_{tool}_{safe_name}_hg{hg}"
    return test_fn


for _variant, _hg, _label in VARIANTS:
    for _tool in ("spliceai", "pangolin"):
        _fn = _make_test(_variant, _hg, _tool)
        setattr(TestAPIConsistency, _fn.__name__, _fn)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capture", action="store_true",
                        help="Capture expected scores from production APIs and write to expected_scores.json")
    args, remaining = parser.parse_known_args()

    if args.capture:
        capture_expected_scores()
    else:
        unittest.main(argv=[sys.argv[0]] + remaining)
