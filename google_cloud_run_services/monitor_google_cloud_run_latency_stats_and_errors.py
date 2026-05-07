#!/usr/bin/env python3
"""Cloud Run production health monitor for spliceai-lookup-412920.

Loops forever, printing every --interval minutes (default 30):
  - error log breakdown by signature
  - response-code totals (2xx/3xx/4xx/5xx)
  - CPU and memory utilization (p95/p99)
  - request latency (p50/p95/p99) with sample counts

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
import json
import subprocess
import time
from datetime import datetime, timezone, timedelta

from google.cloud import monitoring_v3


PROJECT = "spliceai-lookup-412920"
REGION = "us-central1"
SERVICES = ["liftover", "spliceai-37", "spliceai-38", "pangolin-37", "pangolin-38"]


def parse_iso(s):
    return datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)


def get_production_revisions():
    """Return {service_name: revision_name} for the revision currently serving 100% production traffic.

    Picks any traffic entry with `percent == 100`, including a `--tag`-decorated one
    (e.g. after `--to-tags dev=100` promotes a dev revision to production).
    """
    out = {}
    for svc in SERVICES:
        proc = subprocess.run([
            "gcloud", "run", "services", "describe", svc,
            f"--project={PROJECT}", f"--region={REGION}", "--format=json",
        ], capture_output=True, text=True, check=True)
        for entry in json.loads(proc.stdout)["status"].get("traffic", []):
            if entry.get("percent", 0) != 100:
                continue
            out[svc] = entry["revisionName"]
            break
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


def percentiles(client, metric, start, end, revisions=None):
    """Return {service: {p50, p95, p99}} for a DELTA+DISTRIBUTION metric over the window.

    If `revisions` is provided, only those revisions contribute (filters out dev/test traffic).
    """
    span = max(60, int((end - start).total_seconds()))
    f = f'metric.type="{metric}"'
    if revisions:
        f += f" AND {revision_filter_clause(revisions)}"
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


def request_counts(client, start, end, revisions=None):
    """Return {service: {2xx,3xx,4xx,5xx}} sums over the window."""
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
            group_by_fields=["resource.label.service_name", "metric.label.response_code_class"],
        ),
    }):
        svc = ts.resource.labels.get("service_name", "?")
        cls = ts.metric.labels.get("response_code_class", "?")
        for p in ts.points:
            out.setdefault(svc, {})[cls] = out.get(svc, {}).get(cls, 0) + p.value.int64_value
    return out


def fmt_pct(x):
    return f"{x*100:5.1f}%" if x is not None else "  ?  "


def fmt_s(x):
    return f"{x/1000:6.2f}s" if x is not None else "  ?  "


def snapshot(client, args):
    now = datetime.now(timezone.utc)
    start = now - timedelta(hours=args.window_hours)

    if args.all_revisions:
        prod_revs = None
        rev_label = "(all revisions, including dev/test traffic)"
    else:
        prod_map = get_production_revisions()
        if not prod_map:
            # No service has a revision at 100% traffic (e.g. all simultaneously rolling out).
            # Fall back explicitly: an empty list would otherwise be falsy in the per-query
            # `if revisions:` guards and silently query all revisions under a misleading label.
            prod_revs = None
            rev_label = "(all revisions — no service had a 100% production revision)"
        else:
            prod_revs = list(prod_map.values())
            rev_label = "production revisions only:  " + ", ".join(f"{s}={r}" for s, r in sorted(prod_map.items()))

    print("=" * 100)
    print(f"Snapshot at {now.strftime('%Y-%m-%d %H:%M:%SZ')}    "
          f"window: {start.strftime('%Y-%m-%d %H:%M:%SZ')} -> {now.strftime('%H:%M:%SZ')} ({args.window_hours:g}h)")
    print(rev_label)
    print("=" * 100)
    print()

    print("=== Errors by signature (window) ===")
    errs, truncated = gcloud_errors(start, revisions=prod_revs)
    if not errs:
        print("  none")
    else:
        for (svc, sig), c in errs.most_common():
            print(f"  {c:3d}  {svc:<14}  {sig}")
    if truncated:
        print("  WARNING: error log query hit 1000-entry cap — older errors in window are truncated.")
    print()

    print("=== Response codes ===")
    codes = request_counts(client, start, now, revisions=prod_revs)
    totals = {svc: sum(codes.get(svc, {}).get(k, 0) for k in ("2xx", "3xx", "4xx", "5xx")) for svc in SERVICES}
    for svc in SERVICES:
        c = codes.get(svc, {})
        total = totals[svc]
        rate = c.get("5xx", 0) / total * 100 if total else 0
        print(f"  {svc:<14}  2xx={c.get('2xx',0):<5} 3xx={c.get('3xx',0):<3} 4xx={c.get('4xx',0):<5} 5xx={c.get('5xx',0):<3}  ({rate:.2f}% 5xx of {total})")
    print()

    print("=== CPU / Memory utilization ===")
    cpu = percentiles(client, "run.googleapis.com/container/cpu/utilizations", start, now, revisions=prod_revs)
    mem = percentiles(client, "run.googleapis.com/container/memory/utilizations", start, now, revisions=prod_revs)
    print(f"  {'service':<14}  {'CPU p95':<8} {'CPU p99':<8}  {'Mem p95':<8} {'Mem p99':<8}")
    for svc in SERVICES:
        cv = cpu.get(svc, {}); mv = mem.get(svc, {})
        print(f"  {svc:<14}  {fmt_pct(cv.get('p95'))}  {fmt_pct(cv.get('p99'))}   {fmt_pct(mv.get('p95'))}  {fmt_pct(mv.get('p99'))}")
    print()

    lat = percentiles(client, "run.googleapis.com/request_latencies", start, now, revisions=prod_revs)
    if args.baseline_end:
        baseline_end = parse_iso(args.baseline_end)
        baseline_start = baseline_end - timedelta(days=args.baseline_days)
        # Baseline window predates current revisions; query unfiltered to capture pre-deploy traffic.
        baseline_lat = percentiles(client, "run.googleapis.com/request_latencies", baseline_start, baseline_end)
        print(f"=== Latency (p50/p95/p99 in s) — vs {args.baseline_days:g}d baseline ending {args.baseline_end} ===")
    else:
        baseline_lat = None
        print("=== Latency (p50/p95/p99 in s) ===")

    header = f"  {'service':<14}  {'n':<5}  {'p50':<22}  {'p95':<22}  {'p99':<22}"
    print(header)
    for svc in SERVICES:
        n = totals[svc]
        cells = [f"{svc:<14}", f"{n:<5}"]
        for k in ("p50", "p95", "p99"):
            post = lat.get(svc, {}).get(k)
            if baseline_lat is None:
                cells.append(f"{fmt_s(post)}")
            else:
                pre = baseline_lat.get(svc, {}).get(k)
                if pre is None or post is None or pre == 0:
                    cells.append(f"{fmt_s(pre)} → {fmt_s(post)}")
                else:
                    cells.append(f"{fmt_s(pre)} → {fmt_s(post)} ({(post-pre)/pre*100:+4.0f}%)")
        print("  " + "  ".join(cells))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--interval", type=float, default=30.0,
                        help="Minutes between snapshots (default: 30). Use --once for a single snapshot.")
    parser.add_argument("--once", action="store_true",
                        help="Print one snapshot and exit (default is to loop forever).")
    parser.add_argument("--window-hours", type=float, default=2.0,
                        help="Length of each snapshot's window in hours (default: 2)")
    parser.add_argument("--baseline-end",
                        help="Optional ISO-8601 (YYYY-MM-DDTHH:MM:SSZ) end of a baseline window for latency comparison.")
    parser.add_argument("--baseline-days", type=float, default=7.0,
                        help="Length of the baseline window in days (default: 7)")
    parser.add_argument("--all-revisions", action="store_true",
                        help="Aggregate metrics across all revisions of each service "
                             "(default: filter to the revision currently serving 100%% production traffic, "
                             "so dev/test traffic against `dev---*` URLs doesn't contaminate metrics).")
    args = parser.parse_args()

    client = monitoring_v3.MetricServiceClient()
    while True:
        try:
            snapshot(client, args)
        except Exception as e:
            import traceback
            print(f"\n[snapshot failed: {type(e).__name__}: {e} — retrying next interval]")
            traceback.print_exc()
        if args.once:
            return
        next_at = datetime.now(timezone.utc) + timedelta(minutes=args.interval)
        print(f"\n[next snapshot at {next_at.strftime('%H:%M:%SZ')} — Ctrl-C to stop]\n", flush=True)
        try:
            time.sleep(args.interval * 60)
        except KeyboardInterrupt:
            return


if __name__ == "__main__":
    main()
