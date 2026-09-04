#!/usr/bin/env python3
"""Cloud Run production health monitor for spliceai-lookup-412920.

Loops forever, printing every --interval minutes (default 30):
  - error log breakdown by signature
  - response-code totals (2xx/3xx/4xx/5xx)
  - CPU and memory utilization (p95/p99)
  - instance count against each service's max-instances ceiling
  - request latency (p50/p95/p99) with sample counts
  - GeneBe and Ensembl VEP API latency, measured by probing them directly (see probe_external_apis)
  - container cold-start count and startup latency (p50/p95/p99)
  - response-cache hit rate over the past 24 hours / 7 days / 30 days, split by tool and genome build

If a baseline window is provided, latency p50/p95/p99 are also compared
to that window so regressions stand out.

Examples:
  ./monitor_google_cloud_run_latency_stats_and_errors.py
  ./monitor_google_cloud_run_latency_stats_and_errors.py --interval 60 --window-hours 6
  ./monitor_google_cloud_run_latency_stats_and_errors.py --once \
      --baseline-end 2026-05-07T04:00:00Z --baseline-days 7
"""

import argparse
import collections
import http.client
import json
import os
import select
import statistics
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone, timedelta

import psycopg2

from google.cloud import monitoring_v3
from google.cloud import bigquery


PROJECT = "spliceai-lookup-412920"
REGION = "us-central1"
# The un-suffixed services serve the basic gene set (their original names, kept so the
# published API URLs stay valid); the "-comprehensive" ones serve the comprehensive gene set.
SERVICES = ["liftover",
            "spliceai-37", "spliceai-38", "pangolin-37", "pangolin-38",
            "spliceai-37-comprehensive", "spliceai-38-comprehensive",
            "pangolin-37-comprehensive", "pangolin-38-comprehensive"]

# Cache populated by `~/.claude/skills/analyze-gcloud-costs/scripts/cost_analysis.py`.
# If absent, the cost section is skipped with a warning instead of failing.
BILLING_CACHE_DIR = os.path.expanduser("~/.cache/analyze-gcloud-costs")

# The Cloud SQL instance server.py logs every scoring request to, read by cache_hit_counts().
# Same instance connect_to_db.sh points at, and the same password file build_and_deploy.py reads.
DB_INSTANCE = "spliceai-lookup-db"
DB_NAME = "spliceai-lookup-db"
DB_USER = "postgres"
DB_PASSWORD_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".pgpass")

# server.py's log() writes exactly one of these per scoring request that gets past validation
# and the rate limiter, so the two together are the denominator of the cache hit rate.
CACHE_EVENT_NAMES = tuple(f"{tool}:{outcome}"
                          for tool in ("spliceai", "pangolin")
                          for outcome in ("from-cache", "computed"))
CACHE_HIT_OUTCOME = "from-cache"

# The external-API table is the one table that measures nothing historical: it calls the APIs
# while the snapshot is printing, so it has no window to report.
PROBE_WINDOW_LABEL = "n/a, live probes"

# log.genome holds the "hg" request parameter, which is "37" or "38"; the site calls those
# builds hg19 and hg38, so the report does too.
GENOME_LABELS = {"37": "hg19", "38": "hg38"}

# Per-request timeouts for the external-API probes below, matching GENEBE_TIMEOUT_MS and
# ENSEMBL_TIMEOUT_MS in index.html so a probe gives up exactly when the page would.
GENEBE_TIMEOUT_SEC = 15
ENSEMBL_TIMEOUT_SEC = 90


def genebe_response_is_usable(body):
    """True if a GeneBe variant-relaxed response carries an annotation the page can use.

    Mirrors the `annotationResponse.variants` check in annotateVariantWithGeneBe (index.html),
    so a 200 that answers nothing counts as a failure here just as it does in the browser. The
    `ok` that the page tests alongside it is not a GeneBe field: makeRequest sets it from the
    HTTP status, which here is already covered by urlopen raising HTTPError on a non-2xx.
    """
    return bool(isinstance(body, dict) and body.get("variants"))


def ensembl_response_is_usable(body):
    """True if an Ensembl VEP response carries the coordinates the page needs.

    Mirrors the `ensemblApiResponseJson[0].vcf_string` check in normalizeVariant (index.html).
    """
    return bool(isinstance(body, list) and body and isinstance(body[0], dict) and body[0].get("vcf_string"))


# The external APIs index.html calls straight from the user's browser. Cloud Monitoring never
# sees these requests, so probing the APIs from here is the only way to know what they cost a
# search. One entry per (API, query shape) the page actually sends; see probe_external_apis.
#
# The probe variants come from the page's own Examples table (index.html), except the hg19
# coordinate, which is the GRCh37 position of the same variant as chr8-140300616-T-G (the liftover
# test_api_consistency.py and test_ui.py use; on hg19 the hg38 position itself has a different
# reference allele).
GENEBE_URL_PREFIX = "https://api.genebe.net/cloud/api-public/v1/variant-relaxed"
ENSEMBL_HGVS_EXAMPLE = "NM_000249.4(MLH1):c.116G>A"
# What genomicHgvs(index.html) writes for the hg38 coordinate example, which is what
# annotateVariantWithEnsembl sends for a plain coordinate search on hg38.
ENSEMBL_GENOMIC_HGVS_EXAMPLE = "chr8:g.140300616T>G"
EXTERNAL_API_PROBES = [
    # (api, query shape, url, timeout_sec, response-usable predicate)
    ("GeneBe", "hg38 coords",
     f"{GENEBE_URL_PREFIX}?variant={urllib.parse.quote('chr8-140300616-T-G')}&genome=hg38",
     GENEBE_TIMEOUT_SEC, genebe_response_is_usable),
    ("GeneBe", "hg19 coords",
     f"{GENEBE_URL_PREFIX}?variant={urllib.parse.quote('8-141310715-T-G')}&genome=hg19",
     GENEBE_TIMEOUT_SEC, genebe_response_is_usable),
    ("GeneBe", "hg38 HGVS",
     f"{GENEBE_URL_PREFIX}?variant={urllib.parse.quote(ENSEMBL_HGVS_EXAMPLE)}&genome=hg38",
     GENEBE_TIMEOUT_SEC, genebe_response_is_usable),
    # Two shapes reach Ensembl. Transcript HGVS goes there to be converted to coordinates, on
    # both builds. Parsed coordinates go there too, but only on hg38 and only for consequences:
    # annotateVariantWithEnsembl races GeneBe for those, writing the coordinates back out as
    # genomic HGVS first (index.html's genomicHgvs). That second shape is the one a plain
    # coordinate search pays for, which is the page's most common hg38 query.
    ("Ensembl VEP", "hg38 g. HGVS",
     f"https://rest.ensembl.org/vep/human/hgvs/{urllib.parse.quote(ENSEMBL_GENOMIC_HGVS_EXAMPLE)}"
     "?content-type=application/json&vcf_string=1",
     ENSEMBL_TIMEOUT_SEC, ensembl_response_is_usable),
    ("Ensembl VEP", "hg38 HGVS",
     f"https://rest.ensembl.org/vep/human/hgvs/{urllib.parse.quote(ENSEMBL_HGVS_EXAMPLE)}"
     "?content-type=application/json&vcf_string=1",
     ENSEMBL_TIMEOUT_SEC, ensembl_response_is_usable),
    ("Ensembl VEP", "hg19 HGVS",
     f"https://grch37.rest.ensembl.org/vep/human/hgvs/{urllib.parse.quote(ENSEMBL_HGVS_EXAMPLE)}"
     "?content-type=application/json&vcf_string=1",
     ENSEMBL_TIMEOUT_SEC, ensembl_response_is_usable),
]


def parse_iso(s):
    return datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)


def describe_production_services():
    """Return ({service: revision at 100% traffic}, {service: max-instances ceiling}).

    Both are read from one `gcloud run services describe` per service, because the same JSON
    carries them: a second round of describes for the ceiling would double the calls and could
    see a different state part-way through a deploy, leaving the table's limit disagreeing with
    the revision its counts came from.

    The revision is any traffic entry with `percent == 100`. The ceiling is the
    autoscaling.knative.dev/maxScale annotation on the service's template, which Cloud Run
    leaves unset when none was ever applied -- that means the account default, not "unlimited",
    so an absent annotation is reported as None and instance_counts()'s table says so rather
    than inventing a number to compare against.

    A service that can't be described is warned about and skipped rather than raising: the
    deploy workflow runs one job per tool/genome, so a partial or in-progress rollout can
    legitimately leave some of the services in SERVICES absent, and one missing service must
    not blank out the error counts and latencies of the ones that are up.
    """
    revisions = {}
    limits = {}
    for svc in SERVICES:
        proc = subprocess.run([
            "gcloud", "run", "services", "describe", svc,
            f"--project={PROJECT}", f"--region={REGION}", "--format=json",
        ], capture_output=True, text=True)
        if proc.returncode != 0:
            print(f"  WARNING: skipping {svc}: {proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else 'gcloud describe failed'}")
            continue
        described = json.loads(proc.stdout)
        for entry in described["status"].get("traffic", []):
            if entry.get("percent", 0) != 100:
                continue
            revisions[svc] = entry["revisionName"]
            break
        max_scale = (described.get("spec", {}).get("template", {}).get("metadata", {})
                     .get("annotations", {}).get("autoscaling.knative.dev/maxScale"))
        if max_scale is not None:
            limits[svc] = int(max_scale)
    return revisions, limits


def revisions_with_traffic(client, start, end):
    """Return {service_name: {revision_name: request_count}} for the window.

    Every revision that answered a request, production or tagged, so callers can widen the
    revision filter past the one revision that happens to be at 100% right now.
    """
    span = max(60, int((end - start).total_seconds()))
    out = {}
    for ts in client.list_time_series(request={
        "name": f"projects/{PROJECT}",
        "filter": 'metric.type="run.googleapis.com/request_count"',
        "interval": monitoring_v3.TimeInterval(end_time=end, start_time=start),
        "view": monitoring_v3.ListTimeSeriesRequest.TimeSeriesView.FULL,
        "aggregation": monitoring_v3.Aggregation(
            alignment_period={"seconds": span},
            per_series_aligner=monitoring_v3.Aggregation.Aligner.ALIGN_SUM,
            cross_series_reducer=monitoring_v3.Aggregation.Reducer.REDUCE_SUM,
            group_by_fields=["resource.label.service_name", "resource.label.revision_name"],
        ),
    }):
        svc = ts.resource.labels.get("service_name", "?")
        rev = ts.resource.labels.get("revision_name", "?")
        for p in ts.points:
            out.setdefault(svc, {})[rev] = out.get(svc, {}).get(rev, 0) + p.value.int64_value
    return out


def tagged_revisions(start):
    """Return the set of revision names that served a tagged (non-production) URL since `start`.

    Cloud Run gives each tagged revision its own hostname of the form
    `<tag>---<service>-<hash>-<region>.a.run.app`, so a request whose URL carries the `---`
    separator went to a tag (dev, profile, ...) rather than to production.

    Deciding this from the traffic rather than from the service's current tag list is what makes
    it correct across a window: a revision that was the `dev` target earlier in the window can
    have been retagged since, and its dev traffic would then be counted as production.

    The query is narrow by construction (only tag hostnames match), so the 1000-entry cap is not
    a practical concern the way it is for gcloud_errors().
    """
    out = subprocess.run([
        "gcloud", "logging", "read",
        f'resource.type="cloud_run_revision" AND httpRequest.requestUrl:"---" AND '
        f'timestamp>="{start.strftime("%Y-%m-%dT%H:%M:%SZ")}"',
        f"--project={PROJECT}",
        "--format=json", "--limit=1000",
    ], capture_output=True, text=True, check=True)
    return {
        e.get("resource", {}).get("labels", {}).get("revision_name")
        for e in json.loads(out.stdout or "[]")
    } - {None}


def production_revisions(client, start, end):
    """Return ({service: [revision, ...]}, {service: current_100_percent_revision}, {service: max-instances}).

    The first mapping is every revision that served production traffic during the window, which
    is what the metric and log queries must filter on. Filtering on the single revision at 100%
    instead drops the whole pre-deploy part of the window: on 2026-09-01 a mid-window redeploy
    left spliceai-37 reporting 325 of the 1205 requests it actually served, and with them the
    dozens of 5xx it returned during that morning's container-start failures, so the service
    read as 0.00% 5xx while it was failing.

    A revision that only ever answered on a tag hostname is excluded, which is the dev/test
    traffic the revision filter existed to keep out in the first place.
    """
    current, limits = describe_production_services()
    served = revisions_with_traffic(client, start, end)
    tagged = tagged_revisions(start)
    out = {}
    for svc in SERVICES:
        revs = {r for r in served.get(svc, {}) if r not in tagged}
        # Keep the current revision even with no requests yet: the CPU, memory and startup
        # metrics still have something to say about an instance that is up but idle.
        if svc in current:
            revs.add(current[svc])
        if revs:
            out[svc] = sorted(revs)
    return out, current, limits


def uptime_check_accepted_codes():
    """Return {service_name: set(status_code)} to treat as that service's own uptime probes.

    Read from the live config rather than hard-coded, so the report is correct on both sides of
    the migration to /uptime/ and keeps following the checks if an accepted status changes again.

    An accepted 200 is excluded: real traffic returns 200 too, so suppressing it would hide the
    successes this report exists to count (liftover's check probes / and accepts 200). Any other
    accepted status is suppressed, which is exact once a check accepts the 204 that server.py's
    /uptime/ endpoint answers with, because nothing else in this API returns 204. It is not exact
    while a check still probes the parameterless /spliceai/ or /pangolin/ and accepts the 400 the
    server answers there: a genuine client 400 is then pooled in with ~1,700 probes/day/service
    and lands in the ignored count instead of the 4xx total. Repointing the checks at /uptime/ is
    what ends that; until then the ignored count printed below is where those 400s go.
    """
    proc = subprocess.run([
        "gcloud", "monitoring", "uptime", "list-configs",
        f"--project={PROJECT}", "--format=json",
    ], capture_output=True, text=True)
    if proc.returncode != 0:
        print("  WARNING: could not read uptime check configs; uptime-check status codes will "
              "be counted as errors below.")
        return {}
    out = {}
    for cfg in json.loads(proc.stdout or "[]"):
        host = cfg.get("monitoredResource", {}).get("labels", {}).get("host", "")
        # "spliceai-37" also prefixes "spliceai-37-comprehensive", so the longest match wins.
        svc = max((s for s in SERVICES if host.startswith(f"{s}-")), key=len, default=None)
        if not svc:
            continue
        codes = {
            c["statusValue"]
            for c in cfg.get("httpCheck", {}).get("acceptedResponseStatusCodes", [])
            if "statusValue" in c
        }
        out.setdefault(svc, set()).update(c for c in codes if c != 200)
    return out


def revision_filter_clause(revisions):
    # Cloud Monitoring filter syntax: singular `resource.label.<KEY>` (verified working).
    # Cloud Logging uses plural `resource.labels.<KEY>` — see gcloud_errors() below.
    return "(" + " OR ".join(f'resource.label.revision_name="{r}"' for r in revisions) + ")"


def gcloud_errors(start, revisions=None):
    """Fetch ERROR-severity Cloud Run logs since `start` and group by service+signature.

    If `revisions` is provided, only errors on those revisions are counted.
    """
    f = (f'resource.type="cloud_run_revision" AND severity>=ERROR AND '
         f'timestamp>="{start.strftime("%Y-%m-%dT%H:%M:%SZ")}"')
    if revisions:
        rev_clause = " OR ".join(f'resource.labels.revision_name="{r}"' for r in revisions)
        f += f" AND ({rev_clause})"
    out = subprocess.run([
        "gcloud", "logging", "read", f,
        f"--project={PROJECT}",
        "--format=json", "--limit=1000",
    ], capture_output=True, text=True, check=True)
    entries = json.loads(out.stdout or "[]")
    by_sig = collections.Counter()
    for e in entries:
        svc = e.get("resource", {}).get("labels", {}).get("service_name", "?")
        msg = (e.get("textPayload") or e.get("jsonPayload", {}).get("message") or "").strip()
        if "SystemError" in msg:
            sig = "SystemError(Fasta)"
        elif "M does not exist" in msg:
            sig = "KeyError chrM"
        elif "KeyError" in msg:
            sig = "KeyError other"
        elif "malformed or connection" in msg:
            sig = "malformed/connection"
        elif "Traceback" in msg:
            sig = "Traceback (other)"
        elif not msg:
            sig = "(empty payload — likely platform-level)"
        else:
            sig = msg.split("\n")[0][:80]
        by_sig[(svc, sig)] += 1
    return by_sig, len(entries) >= 1000


def percentiles(client, metric, start, end, revisions=None, extra_filter=None):
    """Return {service: {p50, p95, p99}} for a DELTA+DISTRIBUTION metric over the window.

    If `revisions` is provided, only those revisions contribute (filters out dev/test traffic).
    `extra_filter` is appended to the Cloud Monitoring filter, for metrics that carry a label
    worth excluding -- see the uptime-probe exclusion on request_latencies in snapshot().
    """
    span = max(60, int((end - start).total_seconds()))
    f = f'metric.type="{metric}"'
    if revisions:
        f += f" AND {revision_filter_clause(revisions)}"
    if extra_filter:
        f += f" AND {extra_filter}"
    out = {}
    # ALIGN_DELTA + REDUCE_PERCENTILE_X computes the true percentile from the pooled
    # distribution across all series sharing service_name. ALIGN_PERCENTILE_X + REDUCE_MEAN
    # would average per-series percentiles — wrong for request_latencies, which is split
    # by response_code_class (so the mean of per-class p95s is not the service-level p95).
    for name, red in [
        ("p50", monitoring_v3.Aggregation.Reducer.REDUCE_PERCENTILE_50),
        ("p95", monitoring_v3.Aggregation.Reducer.REDUCE_PERCENTILE_95),
        ("p99", monitoring_v3.Aggregation.Reducer.REDUCE_PERCENTILE_99),
    ]:
        for ts in client.list_time_series(request={
            "name": f"projects/{PROJECT}",
            "filter": f,
            "interval": monitoring_v3.TimeInterval(end_time=end, start_time=start),
            "view": monitoring_v3.ListTimeSeriesRequest.TimeSeriesView.FULL,
            "aggregation": monitoring_v3.Aggregation(
                alignment_period={"seconds": span},
                per_series_aligner=monitoring_v3.Aggregation.Aligner.ALIGN_DELTA,
                cross_series_reducer=red,
                group_by_fields=["resource.label.service_name"],
            ),
        }):
            svc = ts.resource.labels.get("service_name", "?")
            for p in ts.points:
                out.setdefault(svc, {})[name] = p.value.double_value
    return out


def sample_count(client, metric, start, end, revisions=None, extra_filter=None):
    """Return {service: total_sample_count} for a DELTA+DISTRIBUTION metric over the window.

    Sums the per-bucket counts across the distribution so we know how many raw observations
    contributed to each percentile in `percentiles()` above, so it takes the same
    `extra_filter` -- a count over a different population would not describe those percentiles.
    """
    span = max(60, int((end - start).total_seconds()))
    f = f'metric.type="{metric}"'
    if revisions:
        f += f" AND {revision_filter_clause(revisions)}"
    if extra_filter:
        f += f" AND {extra_filter}"
    out = {}
    for ts in client.list_time_series(request={
        "name": f"projects/{PROJECT}",
        "filter": f,
        "interval": monitoring_v3.TimeInterval(end_time=end, start_time=start),
        "view": monitoring_v3.ListTimeSeriesRequest.TimeSeriesView.FULL,
        "aggregation": monitoring_v3.Aggregation(
            alignment_period={"seconds": span},
            per_series_aligner=monitoring_v3.Aggregation.Aligner.ALIGN_DELTA,
            cross_series_reducer=monitoring_v3.Aggregation.Reducer.REDUCE_SUM,
            group_by_fields=["resource.label.service_name"],
        ),
    }):
        svc = ts.resource.labels.get("service_name", "?")
        for p in ts.points:
            out[svc] = out.get(svc, 0) + int(p.value.distribution_value.count)
    return out


def instance_counts(client, start, end, revisions=None):
    """Return {service: [instance count per minute, ...]} over the window.

    ALIGN_MAX within each series and then REDUCE_SUM across them, because instance_count is
    split by a `state` label (active/idle) as well as by revision while max-instances caps the
    total, so answering "did it reach the ceiling" means adding those back together. The max
    rather than the mean within each minute keeps a brief scale-out from being averaged away,
    which is the event this is read for.

    Note for the caller: Cloud Run reports nothing while a service is scaled to zero, so the
    returned list is the minutes the service was actually running, not the whole window. A
    percentage computed from it is a share of the service's running time -- which is the useful
    denominator here, since idle minutes cannot be near a ceiling.
    """
    f = f'metric.type="run.googleapis.com/container/instance_count"'
    if revisions:
        f += f" AND {revision_filter_clause(revisions)}"
    out = {}
    for ts in client.list_time_series(request={
        "name": f"projects/{PROJECT}",
        "filter": f,
        "interval": monitoring_v3.TimeInterval(end_time=end, start_time=start),
        "view": monitoring_v3.ListTimeSeriesRequest.TimeSeriesView.FULL,
        "aggregation": monitoring_v3.Aggregation(
            alignment_period={"seconds": 60},
            per_series_aligner=monitoring_v3.Aggregation.Aligner.ALIGN_MAX,
            cross_series_reducer=monitoring_v3.Aggregation.Reducer.REDUCE_SUM,
            group_by_fields=["resource.label.service_name"],
        ),
    }):
        svc = ts.resource.labels.get("service_name", "?")
        for p in ts.points:
            out.setdefault(svc, []).append(p.value.double_value or p.value.int64_value)
    return out


def request_counts(client, start, end, revisions=None):
    """Return {service: {response_code_int: count}} sums over the window.

    Grouping by raw `response_code` (rather than `response_code_class`) lets the
    caller derive class totals AND show a per-code breakdown — useful when one
    code (e.g. 404 from probe traffic) dominates the class total and would
    otherwise mask real 400/429/etc. errors.
    """
    span = max(60, int((end - start).total_seconds()))
    f = 'metric.type="run.googleapis.com/request_count"'
    if revisions:
        f += f" AND {revision_filter_clause(revisions)}"
    out = {}
    for ts in client.list_time_series(request={
        "name": f"projects/{PROJECT}",
        "filter": f,
        "interval": monitoring_v3.TimeInterval(end_time=end, start_time=start),
        "view": monitoring_v3.ListTimeSeriesRequest.TimeSeriesView.FULL,
        "aggregation": monitoring_v3.Aggregation(
            alignment_period={"seconds": span},
            per_series_aligner=monitoring_v3.Aggregation.Aligner.ALIGN_SUM,
            cross_series_reducer=monitoring_v3.Aggregation.Reducer.REDUCE_SUM,
            group_by_fields=["resource.label.service_name", "metric.label.response_code"],
        ),
    }):
        svc = ts.resource.labels.get("service_name", "?")
        try:
            code = int(ts.metric.labels.get("response_code", ""))
        except ValueError:
            continue
        for p in ts.points:
            out.setdefault(svc, {})[code] = out.get(svc, {}).get(code, 0) + p.value.int64_value
    return out


def fmt_pct(x):
    return f"{x*100:.1f}%" if x is not None else "?"


def fmt_s(x):
    return f"{x/1000:.2f}s" if x is not None else "?"


def fmt_hit_rate(hits, total):
    """Return "NN.N% (hits/total)", with the rate padded so a column of these lines up.

    The counts vary in width from one row to the next, so the rate is what has to be padded
    (and the column left-aligned) for the percentages to read down the column as a column.
    """
    return f"{hits/total*100:5.1f}% ({hits:,}/{total:,})" if total else f"{'-':>6}"


def discover_billing_export():
    """Return the fully-qualified BigQuery billing-export table for PROJECT, or None.

    Reads the cache populated by `/analyze-gcloud-costs`. If the cache file or
    billing-account lookup is unavailable, returns None so the caller can skip
    the cost section instead of failing the whole snapshot.
    """
    try:
        out = subprocess.run([
            "gcloud", "billing", "projects", "describe", PROJECT,
            "--format=value(billingAccountName)",
        ], capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    if not out:
        return None
    cache_file = os.path.join(BILLING_CACHE_DIR, f"{out.replace('billingAccounts/', '')}.json")
    if not os.path.exists(cache_file):
        return None
    with open(cache_file) as f:
        e = json.load(f)
    return f"{e['bq_project']}.{e['dataset']}.{e['standard_table']}"


def daily_costs(bq_client, table, days):
    """Return (daily_totals, top_skus) for PROJECT over the last `days` days.

    daily_totals: list of (date, usd) ordered by date ascending.
    top_skus:     list of (service, sku, usd) ordered by usd descending, top 10.
    Net cost = gross usage cost + credits (credits are negative).
    """
    params = [
        bigquery.ScalarQueryParameter("project_id", "STRING", PROJECT),
        bigquery.ScalarQueryParameter("days", "INT64", days),
    ]
    # Pad partition filter by 5 days for late-arriving billing data. The usage_start_time
    # comparison uses `>` (not `>=`) so `--cost-days N` returns exactly N rows: today plus
    # the prior N-1 days. Today's row is partial — that's intentional, so the chart shows
    # spending-so-far rather than dropping to zero at the right edge.
    partition_clause = (
        "_PARTITIONTIME >= TIMESTAMP_SUB(CURRENT_TIMESTAMP(), INTERVAL (@days + 5) DAY) "
        "AND DATE(usage_start_time) > DATE_SUB(CURRENT_DATE(), INTERVAL @days DAY) "
        "AND project.id = @project_id"
    )
    daily = [(r.day, float(r.net_usd or 0)) for r in bq_client.query(f"""
        SELECT
          DATE(usage_start_time) AS day,
          ROUND(SUM(cost) + SUM(IFNULL((SELECT SUM(c.amount) FROM UNNEST(credits) c), 0)), 4) AS net_usd
        FROM `{table}`
        WHERE {partition_clause}
        GROUP BY day
        ORDER BY day
    """, job_config=bigquery.QueryJobConfig(query_parameters=params)).result()]
    skus = [(r.service, r.sku, float(r.net_usd or 0)) for r in bq_client.query(f"""
        SELECT
          service.description AS service,
          sku.description AS sku,
          ROUND(SUM(cost) + SUM(IFNULL((SELECT SUM(c.amount) FROM UNNEST(credits) c), 0)), 2) AS net_usd
        FROM `{table}`
        WHERE {partition_clause}
        GROUP BY service, sku
        HAVING net_usd >= 0.01
        ORDER BY net_usd DESC
        LIMIT 10
    """, job_config=bigquery.QueryJobConfig(query_parameters=params)).result()]
    return daily, skus


def print_cost_chart(daily, skus, bar_width=40):
    """Print an ASCII bar chart of daily net cost plus the top-SKU breakdown.

    Bar length is scaled so the largest day fills `bar_width` cells. Uses the
    eight Unicode left-block fractions (U+2589..U+2588) for sub-cell precision
    so small days stay visible without a `log` scale distorting comparisons.
    """
    if not daily:
        print("  (no billing data returned for the window)")
        return
    blocks = " ▏▎▍▌▋▊▉█"
    max_usd = max(usd for _, usd in daily) or 1.0
    total = sum(usd for _, usd in daily)
    for day, usd in daily:
        units = bar_width * 8 * usd / max_usd
        full, frac = divmod(round(units), 8)
        bar = "█" * full + (blocks[frac] if frac else "")
        print(f"  {day.strftime('%a %b %d')}  {bar:<{bar_width}}  ${usd:>7,.2f}")
    print(f"\n  Total over {len(daily)} day(s): ${total:,.2f}   "
          f"(avg ${total/len(daily):,.2f}/day)")

    if skus:
        print()
        print("  Top SKUs:")
        rows = [["$usd", "service", "sku"]]
        for service, sku, usd in skus:
            rows.append([f"${usd:,.2f}", service, sku[:60] + ("…" if len(sku) > 60 else "")])
        print_table(rows, aligns=['r', 'l', 'l'], indent="    ")


def discover_db_connection_params():
    """Return psycopg2 connect kwargs for the log database, or None if it can't be reached.

    The host comes from SPLICEAI_LOOKUP_DB_HOST when that is set (the variable
    build_and_deploy.py reads out of .env), and otherwise from the instance's PRIMARY public IP
    as gcloud reports it, so this works on a machine that has no .env. Returns None rather than
    raising when the password file is missing or gcloud can't describe the instance, so the
    caller skips the cache-hit section instead of failing the whole snapshot.
    """
    if not os.path.exists(DB_PASSWORD_FILE):
        return None
    host = os.environ.get("SPLICEAI_LOOKUP_DB_HOST")
    if not host:
        proc = subprocess.run([
            "gcloud", "sql", "instances", "describe", DB_INSTANCE,
            f"--project={PROJECT}", "--format=json",
        ], capture_output=True, text=True)
        if proc.returncode != 0:
            return None
        host = next((ip["ipAddress"] for ip in json.loads(proc.stdout).get("ipAddresses", [])
                     if ip.get("type") == "PRIMARY"), None)
    if not host:
        return None
    with open(DB_PASSWORD_FILE) as f:
        password = f.read().strip()
    return {"host": host, "dbname": DB_NAME, "user": DB_USER, "password": password,
            "connect_timeout": 20}


def cache_hit_windows(now):
    """Return [(label, start), ...] for the cache-hit report, as naive UTC datetimes.

    Naive rather than timezone-aware because log.logtime is a `TIMESTAMP` (no zone) written by
    the database's own now(), and the instance runs in UTC. Comparing that column against an
    aware value would make the boundary depend on the session's TimeZone setting instead of
    being the plain timestamp comparison it looks like.

    All three end at `now`. The shortest was the calendar UTC day until that turned out to be
    read as a full day whatever the hour: run at 00:33Z it counted 109 requests over its 33
    minutes, where the 24 hours before it held about 19,000.
    """
    naive_now = now.replace(tzinfo=None)
    return [
        ("past 24 hours", naive_now - timedelta(hours=24)),
        ("past 7 days", naive_now - timedelta(days=7)),
        ("past 30 days", naive_now - timedelta(days=30)),
    ]


def cache_hit_counts(connect_params, windows):
    """Count cache hits and total scoring requests per tool, genome build and window.

    Every window is counted in a single pass with one FILTER clause each, because log.logtime
    carries no index (see _SCHEMA_DDL_STATEMENTS in server.py): a query per window would be a
    sequential scan of the whole table per window.

    Args:
        connect_params (dict): psycopg2.connect kwargs, from discover_db_connection_params().
        windows (list): (label, start_datetime) pairs, from cache_hit_windows().

    Returns:
        tuple: (hits, totals, earliest), where hits and totals are Counters keyed by
            (tool, genome, window_label), and earliest is the oldest logtime the query saw
            (None if it matched nothing) -- which says how much of the widest window the log
            table still covers.
    """
    filters = ", ".join(f"count(*) FILTER (WHERE logtime >= %s) AS w{i}"
                        for i in range(len(windows)))
    placeholders = ", ".join(["%s"] * len(CACHE_EVENT_NAMES))
    params = ([start for _, start in windows]
              + [min(start for _, start in windows)]
              + list(CACHE_EVENT_NAMES))
    hits, totals = collections.Counter(), collections.Counter()
    earliest = None
    conn = psycopg2.connect(**connect_params)
    try:
        # psycopg2's connection context manager ends the transaction but leaves the connection
        # open, so the close below is what actually releases it -- and this runs every snapshot.
        with conn.cursor() as cursor:
            cursor.execute(f"""
                SELECT split_part(event_name, ':', 1) AS tool,
                       genome,
                       split_part(event_name, ':', 2) AS outcome,
                       min(logtime) AS earliest,
                       {filters}
                FROM log
                WHERE logtime >= %s AND event_name IN ({placeholders})
                GROUP BY 1, 2, 3
            """, params)
            for tool, genome, outcome, group_earliest, *counts in cursor.fetchall():
                earliest = group_earliest if earliest is None else min(earliest, group_earliest)
                for (label, _), count in zip(windows, counts):
                    totals[(tool, genome, label)] += count
                    if outcome == CACHE_HIT_OUTCOME:
                        hits[(tool, genome, label)] += count
    finally:
        conn.close()
    return hits, totals, earliest


def print_cache_hit_rates(hits, totals, windows, earliest):
    """Print the cache hit rate per tool and genome build, one column per window."""
    if not totals:
        print("  (no scoring requests logged over these windows)")
        return
    labels = [label for label, _ in windows]
    rows = [["tool", "genome"] + labels]
    for tool, genome in sorted({(tool, genome) for tool, genome, _ in totals}):
        rows.append([tool, GENOME_LABELS.get(genome, f"hg={genome}")] + [
            fmt_hit_rate(hits[(tool, genome, label)], totals[(tool, genome, label)])
            for label in labels
        ])
    rows.append(["all", "all"] + [
        fmt_hit_rate(sum(c for (_, _, l), c in hits.items() if l == label),
                     sum(c for (_, _, l), c in totals.items() if l == label))
        for label in labels
    ])
    print_table(rows, aligns=['l', 'l'] + ['l'] * len(labels))
    widest_label, widest_start = min(windows, key=lambda w: w[1])
    # A window can outrun what the log table still holds (it currently starts at 2026-07-15).
    # Without this, the widest column would keep its header and quietly describe a shorter
    # period than the header claims.
    if earliest is not None and earliest - widest_start > timedelta(days=1):
        print(f"  NOTE: the oldest matching log row is from {earliest.strftime('%Y-%m-%d %H:%MZ')}, "
              f"so \"{widest_label}\" only reaches back to then.")


def probe_external_apis(samples, probes=EXTERNAL_API_PROBES):
    """Time each external-API query shape in `probes` by calling it `samples` times.

    Requests are issued serially, one shape fully finished before the next starts, because
    GeneBe throttles concurrent requests (measured while comparing it to Ensembl on 2026-08-16:
    calls that failed inside a 5-worker run succeeded when replayed alone).

    A call is only counted as successful if the response is one the page could actually use, so
    an HTTP 200 carrying no annotation is a failure here exactly as it is in the browser.

    Args:
        samples (int): number of requests to make per query shape.
        probes (list): (api, shape, url, timeout_sec, is_usable) tuples.

    Returns:
        list: one (api, shape, results) triple per query shape, where results is a list of
            (elapsed_seconds, outcome) and outcome is "ok" or a short failure description.
    """
    out = []
    for api, shape, url, timeout, is_usable in probes:
        results = []
        for _ in range(samples):
            # Identify the probe rather than sending urllib's default "Python-urllib/3.x",
            # which some APIs rate-limit or reject outright.
            request = urllib.request.Request(url, headers={
                "User-Agent": "SpliceAI-lookup-monitor/1.0 (+https://spliceailookup.broadinstitute.org)",
                "Accept": "application/json",
            })
            started = time.monotonic()
            try:
                with urllib.request.urlopen(request, timeout=timeout) as response:
                    outcome = "ok" if is_usable(json.loads(response.read())) else "unusable response"
            except urllib.error.HTTPError as e:
                outcome = f"HTTP {e.code}"
            except urllib.error.URLError as e:
                # A timeout during connect or read surfaces wrapped in URLError, while one
                # during the TLS handshake can be raised bare, hence the separate clause below.
                outcome = (f"timeout >{timeout:g}s" if isinstance(e.reason, TimeoutError)
                           else f"URLError ({e.reason})")
            except TimeoutError:
                outcome = f"timeout >{timeout:g}s"
            except (OSError, http.client.HTTPException) as e:
                # urllib only wraps errors from sending the request in URLError. Anything raised
                # while reading the response (RemoteDisconnected, ConnectionResetError,
                # IncompleteRead) comes through bare, and one flaky call should show up in the
                # failures column rather than abort the whole snapshot.
                outcome = type(e).__name__
            except (json.JSONDecodeError, UnicodeDecodeError):
                outcome = "unparseable response"
            results.append((time.monotonic() - started, outcome))
        out.append((api, shape, results))
    return out


def print_external_api_latencies(probe_results):
    """Print the per-query-shape latency table for the external APIs.

    min/median/max cover the successful calls only, so a shape that is entirely failing shows
    "-" rather than the time its failures took to come back. The failures are summarized in the
    last column with their own timing, since a 503 after 10s and a 503 after 0.2s say very
    different things about where the API is broken.
    """
    rows = [["api", "query shape", "n", "ok", "min", "med", "max", "failures"]]
    for api, shape, results in probe_results:
        ok = sorted(sec for sec, outcome in results if outcome == "ok")
        failures = collections.Counter(outcome for _, outcome in results if outcome != "ok")
        notes = []
        for outcome, count in failures.most_common():
            elapsed = [sec for sec, o in results if o == outcome]
            notes.append(f"{count}x {outcome} (~{statistics.median(elapsed):.2f}s)")
        rows.append([
            api,
            shape,
            str(len(results)),
            str(len(ok)),
            fmt_s(ok[0] * 1000) if ok else "-",
            fmt_s(statistics.median(ok) * 1000) if ok else "-",
            fmt_s(ok[-1] * 1000) if ok else "-",
            ", ".join(notes),
        ])
    print_table(rows, aligns=['l', 'l', 'r', 'r', 'r', 'r', 'r', 'l'])


def print_section_header(title, window):
    """Print a table's `=== title (window: ...) ===` line.

    Every table spells out its own window because they do not all share one: most follow
    --window-hours, the cache table has its own three and the cost chart follows --cost-days.
    Without this the window named once in the snapshot header reads as if it covered every
    table under it.

    Args:
        title (str): what the table shows, without the `===` markers.
        window (str): the period it covers, e.g. "last 12h" or "n/a, live probes".
    """
    print(f"=== {title} (window: {window}) ===")


def print_table(rows, aligns=None, indent="  ", gap="  "):
    """Print rows aligned by max column widths.

    rows: list of lists of cell strings (first row treated as header — same alignment).
    aligns: list of 'l' or 'r' per column (default: all 'l').
    """
    if not rows:
        return
    n_cols = len(rows[0])
    aligns = aligns or ['l'] * n_cols
    widths = [max(len(row[i]) for row in rows) for i in range(n_cols)]
    for row in rows:
        cells = [
            row[i].ljust(widths[i]) if aligns[i] == 'l' else row[i].rjust(widths[i])
            for i in range(n_cols)
        ]
        print((indent + gap.join(cells)).rstrip())


def snapshot(client, args, bq_client=None, billing_table=None, db_connect_params=None):
    now = datetime.now(timezone.utc)
    start = now - timedelta(hours=args.window_hours)
    window_label = f"last {args.window_hours:g}h"

    if args.all_revisions:
        prod_revs = None
        max_instances = {}
        rev_label = "(all revisions, including dev/test traffic)"
    else:
        prod_map, current_map, max_instances = production_revisions(client, start, now)
        if not prod_map:
            # Either nothing served production traffic and no service has a revision at 100%
            # (e.g. all simultaneously rolling out), or every describe call failed. Fall back
            # explicitly: an empty list would otherwise be falsy in the per-query `if revisions:`
            # guards and silently query all revisions under a misleading label. Say plainly that
            # dev/test traffic is now mixed in, since these numbers are no longer production's.
            prod_revs = None
            rev_label = ("WARNING: could not identify a production revision for ANY service, so the "
                         "figures below cover ALL revisions and include dev/test traffic.")
        else:
            prod_revs = sorted({r for revs in prod_map.values() for r in revs})
            lines = ["production revisions only (dev/test tags excluded):"]
            for svc in sorted(prod_map):
                revs = prod_map[svc]
                # Name the current revision first, then say how many older ones the window also
                # covers. More than one means a deploy landed mid-window; before this was fixed
                # everything before that deploy was silently dropped from every figure below.
                head = current_map.get(svc, revs[-1])
                extra = [r for r in revs if r != head]
                note = f"  (+{len(extra)} earlier this window: {', '.join(extra)})" if extra else ""
                lines.append(f"  {svc:<26} {head}{note}")
            rev_label = "\n".join(lines)

    print("=" * 100)
    print(f"Snapshot at {now.strftime('%Y-%m-%d %H:%M:%SZ')}    "
          f"window: {start.strftime('%Y-%m-%d %H:%M:%SZ')} -> {now.strftime('%H:%M:%SZ')} ({args.window_hours:g}h)")
    print(rev_label)
    print("=" * 100)
    print()

    print_section_header("Errors by signature", window_label)
    errs, truncated = gcloud_errors(start, revisions=prod_revs)
    if not errs:
        print("  none")
    else:
        for (svc, sig), c in errs.most_common():
            print(f"  {c:3d}  {svc:<14}  {sig}")
    if truncated:
        print("  WARNING: error log query hit 1000-entry cap — older errors in window are truncated.")
    print()

    print_section_header("Response codes (3xx, 404s and each service's uptime-check status "
                         "ignored as probe/redirect noise)", window_label)
    codes = request_counts(client, start, now, revisions=prod_revs)
    accepted = uptime_check_accepted_codes()
    for svc in SERVICES:
        all_codes = codes.get(svc, {})
        # 404s are Google's TsunamiSecurityScanner walking /etc/passwd, /login and .jsp upload
        # paths; a client that mistypes an endpoint lands here too, so the count is still printed
        # in the ignored note rather than thrown away. The rest is whatever this service's own
        # check calls a pass — see uptime_check_accepted_codes for why that is read live and what
        # it still costs before the /uptime/ migration. 405 is counted normally: not every one is
        # the scanner, so it is not safe to drop on the assumption that it is noise.
        noise = {404} | accepted.get(svc, set())
        by_code = {code: c for code, c in all_codes.items() if code not in noise and code // 100 != 3}
        classes = collections.Counter()
        for code, c in by_code.items():
            classes[f"{code // 100}xx"] += c
        total = sum(classes.values())
        rate = classes["5xx"] / total * 100 if total else 0
        rate_2xx = classes["2xx"] / total * 100 if total else 0
        ignored_items = sorted(
            ((code, c) for code, c in all_codes.items() if code in noise and c),
            key=lambda x: -x[1],
        )
        ignored_note = ("; +" + ", ".join(f"{c} {code}s" for code, c in ignored_items) + " ignored"
                        if ignored_items else "")
        print(f"  {svc:<14}  2xx={classes['2xx']:<5} ({rate_2xx:5.1f}%) "
              f"4xx={classes['4xx']:<5} 5xx={classes['5xx']:<3}  "
              f"({rate:.2f}% 5xx of {total}{ignored_note})")
        for cls in ("4xx", "5xx"):
            items = sorted(
                ((code, c) for code, c in by_code.items() if code // 100 == int(cls[0]) and c),
                key=lambda x: -x[1],
            )
            if items:
                print(f"                  {cls}: " + ", ".join(f"{code}={c}" for code, c in items))
    print()

    print_section_header("CPU / Memory utilization", window_label)
    cpu_metric = "run.googleapis.com/container/cpu/utilizations"
    mem_metric = "run.googleapis.com/container/memory/utilizations"
    cpu = percentiles(client, cpu_metric, start, now, revisions=prod_revs)
    mem = percentiles(client, mem_metric, start, now, revisions=prod_revs)
    cpu_n = sample_count(client, cpu_metric, start, now, revisions=prod_revs)
    mem_n = sample_count(client, mem_metric, start, now, revisions=prod_revs)
    rows = [["service", "CPU n", "CPU p95", "CPU p99", "Mem n", "Mem p95", "Mem p99"]]
    for svc in SERVICES:
        cv = cpu.get(svc, {}); mv = mem.get(svc, {})
        rows.append([
            svc,
            str(cpu_n.get(svc, 0)),
            fmt_pct(cv.get('p95')),
            fmt_pct(cv.get('p99')),
            str(mem_n.get(svc, 0)),
            fmt_pct(mv.get('p95')),
            fmt_pct(mv.get('p99')),
        ])
    print_table(rows, aligns=['l', 'r', 'r', 'r', 'r', 'r', 'r'])
    print()

    # Whether the autoscaler ever ran out of room. A service pinned at its ceiling queues
    # requests inside its instances instead of scaling out, which shows up as latency rather
    # than as an error, so nothing else in this report would name it. "mins" is the minutes the
    # service was running (Cloud Run reports nothing while scaled to zero), and the percentages
    # are shares of that, not of the window.
    print_section_header("Instances vs the max-instances ceiling", window_label)
    inst = instance_counts(client, start, now, revisions=prod_revs)
    rows = [["service", "limit", "mins", "peak", "at limit", "within 1", "mean"]]
    for svc in SERVICES:
        vals = inst.get(svc, [])
        limit = max_instances.get(svc)
        limit_s = str(limit) if limit else "?"
        if not vals:
            rows.append([svc, limit_s, "0", "?", "?", "?", "?"])
            continue
        n = len(vals)
        if limit:
            at = sum(1 for v in vals if v >= limit)
            at_s = f"{at} ({at / n * 100:.1f}%)"
            # "within 1" only says something once limit-1 is at least 2: a running service
            # always has an instance, so on liftover's limit of 2 it would report every minute
            # as near the ceiling and read as saturation that is not there.
            if limit >= 3:
                near = sum(1 for v in vals if v >= limit - 1)
                near_s = f"{near} ({near / n * 100:.1f}%)"
            else:
                near_s = "-"
        else:
            # Either --all-revisions skipped the describe calls, or the service carries no
            # maxScale annotation. Report the counts that stand on their own rather than
            # comparing against a ceiling this run does not know.
            at_s = near_s = "-"
        rows.append([svc, limit_s, str(n), f"{max(vals):.0f}", at_s, near_s, f"{sum(vals) / n:.2f}"])
    print_table(rows, aligns=['l', 'r', 'r', 'r', 'r', 'r', 'r'])
    print()

    lat_metric = "run.googleapis.com/request_latencies"
    # Drop the uptime probes here too, for the same reason the response-code table drops them:
    # each scoring service takes roughly 1,700 a day, which on the quiet comprehensive services
    # is more than their entire real traffic, and a probe returns without doing any scoring work.
    # Left in, the p50 would describe the probe rather than the requests this table is read for.
    # The statuses come from the same live check configs (see uptime_check_accepted_codes); the
    # filter is one query for every service, so it is the union rather than each service's own
    # status -- harmless while every scoring check accepts the same one.
    probe_codes = sorted({code for codes in accepted.values() for code in codes})
    lat_filter = " AND ".join(f'metric.label.response_code != "{code}"' for code in probe_codes)
    lat = percentiles(client, lat_metric, start, now, revisions=prod_revs, extra_filter=lat_filter)
    lat_n = sample_count(client, lat_metric, start, now, revisions=prod_revs, extra_filter=lat_filter)
    if args.baseline_end:
        baseline_end = parse_iso(args.baseline_end)
        baseline_start = baseline_end - timedelta(days=args.baseline_days)
        # Baseline window predates current revisions; query unfiltered to capture pre-deploy traffic.
        baseline_lat = percentiles(client, lat_metric, baseline_start, baseline_end,
                                   extra_filter=lat_filter)
        print_section_header("Latency (p50/p95/p99 in s)",
                             f"{window_label} vs a {args.baseline_days:g}d baseline "
                             f"ending {args.baseline_end}")
    else:
        baseline_lat = None
        print_section_header("Latency (p50/p95/p99 in s)", window_label)

    rows = [["service", "n", "p50", "p95", "p99"]]
    for svc in SERVICES:
        cells = [svc, str(lat_n.get(svc, 0))]
        for k in ("p50", "p95", "p99"):
            post = lat.get(svc, {}).get(k)
            if baseline_lat is None:
                cells.append(fmt_s(post))
            else:
                pre = baseline_lat.get(svc, {}).get(k)
                if pre is None or post is None or pre == 0:
                    cells.append(f"{fmt_s(pre)} → {fmt_s(post)}")
                else:
                    cells.append(f"{fmt_s(pre)} → {fmt_s(post)} ({(post-pre)/pre*100:+.0f}%)")
        rows.append(cells)
    print_table(rows, aligns=['l', 'r', 'r', 'r', 'r'])
    print()

    # These calls are made by the user's browser, not by our Cloud Run services, so they appear
    # in none of the metrics above even though a slow GeneBe or Ensembl delays every search that
    # needs one. The numbers below are live probes from this machine rather than real user
    # traffic, so they measure the API rather than what any particular user experienced.
    if args.probe_samples < 1:
        print_section_header("External API latency (skipped, --probe-samples 0)", PROBE_WINDOW_LABEL)
    else:
        print_section_header(f"External API latency ({args.probe_samples} probe"
                             f"{'' if args.probe_samples == 1 else 's'} per query shape, from this machine)",
                             PROBE_WINDOW_LABEL)
        print_external_api_latencies(probe_external_apis(args.probe_samples))
    print()

    # Container startup latency distribution — count = number of cold starts in window,
    # percentiles = how long each new instance took to become ready to serve requests.
    print_section_header("Container startup (cold starts; p50/p95/p99 in s)", window_label)
    startup_metric = "run.googleapis.com/container/startup_latencies"
    startup = percentiles(client, startup_metric, start, now, revisions=prod_revs)
    startup_n = sample_count(client, startup_metric, start, now, revisions=prod_revs)
    rows = [["service", "n", "p50", "p95", "p99"]]
    for svc in SERVICES:
        sv = startup.get(svc, {})
        rows.append([
            svc,
            str(startup_n.get(svc, 0)),
            fmt_s(sv.get('p50')),
            fmt_s(sv.get('p95')),
            fmt_s(sv.get('p99')),
        ])
    print_table(rows, aligns=['l', 'r', 'r', 'r', 'r'])
    print()

    # A cache hit never reaches the model, so it costs a fraction of a miss -- this rate is what
    # says how much of the traffic the scoring services actually have to compute. It comes from
    # the `log` table server.py writes, not from Cloud Monitoring, which cannot tell the two
    # apart. So unlike every table above it covers its own fixed windows rather than
    # --window-hours, and it includes dev/tagged traffic, since the table records no revision.
    # A `force=1` request counts as a miss, which is what it is: the lookup is skipped and the
    # model runs. Basic and comprehensive gene sets are pooled here; the table's `bc` column
    # separates them for anyone who needs that breakdown.

    # Computed before the header rather than inside the else below so the header names the
    # windows even when the database is unreachable and no table follows it.
    windows = cache_hit_windows(now)
    print_section_header("Cache hit rate (share of scoring requests answered from the response cache)",
                         ", ".join(label for label, _ in windows))
    if db_connect_params is None:
        print(f"  (skipped, no database connection: needs {DB_PASSWORD_FILE}, plus either "
              f"SPLICEAI_LOOKUP_DB_HOST or a working "
              f"`gcloud sql instances describe {DB_INSTANCE}`)")
    else:
        try:
            hits, totals, earliest = cache_hit_counts(db_connect_params, windows)
        except psycopg2.Error as e:
            # Degrade like the cost section rather than raising: an unreachable database (the
            # instance's authorized-networks list not covering this machine, say) must not cost
            # the snapshot the cost chart printed after it.
            print(f"  (skipped, could not read the log table: {str(e).strip()})")
        else:
            print_cache_hit_rates(hits, totals, windows, earliest)
    print()

    print_section_header("Project cost (net of credits)", f"last {args.cost_days:g} days")
    if bq_client is None or billing_table is None:
        print("  (skipped — billing-export discovery cache not found at "
              f"{BILLING_CACHE_DIR}; run /analyze-gcloud-costs once to populate it)")
    else:
        daily, skus = daily_costs(bq_client, billing_table, int(args.cost_days))
        print_cost_chart(daily, skus)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--interval", type=float, default=30.0,
                        help="Minutes between snapshots (default: 30). Use --once for a single snapshot.")
    parser.add_argument("--once", action="store_true",
                        help="Print one snapshot and exit (default is to loop forever).")
    parser.add_argument("--window-hours", type=float, default=2.0,
                        help="Length of each snapshot's window in hours (default: 2). Covers every "
                             "table except the external-API probes, the cache hit rate (its own "
                             "24h/7d/30d windows) and the cost chart (--cost-days).")
    parser.add_argument("--baseline-end",
                        help="Optional ISO-8601 (YYYY-MM-DDTHH:MM:SSZ) end of a baseline window for latency comparison.")
    parser.add_argument("--baseline-days", type=float, default=7.0,
                        help="Length of the baseline window in days (default: 7)")
    parser.add_argument("--all-revisions", action="store_true",
                        help="Aggregate metrics across all revisions of each service "
                             "(default: filter to the revisions that served production traffic during "
                             "the window, so dev/test traffic against `dev---*` URLs doesn't "
                             "contaminate metrics).")
    parser.add_argument("--cost-days", type=float, default=14.0,
                        help="Number of days of daily-cost history to chart (default: 14).")
    parser.add_argument("--probe-samples", type=int, default=3,
                        help="Requests to send to each external-API query shape (GeneBe and Ensembl VEP) "
                             "per snapshot (default: 3). Use 0 to skip the probes: they are serial, so "
                             "with both APIs down a snapshot spends several minutes waiting on timeouts.")
    args = parser.parse_args()

    client = monitoring_v3.MetricServiceClient()
    billing_table = discover_billing_export()
    # Reuse a single client across iterations so we don't repeat ADC + project discovery
    # every 30 minutes. The BQ project hosts the billing export (not PROJECT itself).
    bq_client = bigquery.Client(project=billing_table.split(".")[0]) if billing_table else None
    # Discovered once for the same reason: the instance IP costs a gcloud call to look up, and
    # it does not change between snapshots. The connection itself is opened per snapshot.
    db_connect_params = discover_db_connection_params()
    while True:
        print("Processing...")
        try:
            snapshot(client, args, bq_client=bq_client, billing_table=billing_table,
                     db_connect_params=db_connect_params)
        except Exception as e:
            import traceback
            print(f"\n[snapshot failed: {type(e).__name__}: {e} — retrying next interval]")
            traceback.print_exc()
        if args.once:
            return
        next_at = datetime.now(timezone.utc) + timedelta(minutes=args.interval)
        print(f"\n[next snapshot at {next_at.strftime('%H:%M:%SZ')} — Enter to run now, Ctrl-C to stop]\n", flush=True)
        try:
            deadline = time.monotonic() + args.interval * 60
            while True:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    break
                # select() on stdin lets Enter break the wait early without
                # putting the terminal into raw mode. stdin stays line-buffered,
                # so it only becomes readable after a newline is typed.
                ready, _, _ = select.select([sys.stdin], [], [], remaining)
                if ready:
                    sys.stdin.readline()
                    break
        except KeyboardInterrupt:
            return


if __name__ == "__main__":
    main()
