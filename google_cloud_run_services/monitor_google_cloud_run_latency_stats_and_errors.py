#!/usr/bin/env python3
"""Cloud Run production health monitor for spliceai-lookup-412920.

Loops forever, printing every --interval minutes (default 30):
  - error log breakdown by signature
  - response-code totals (2xx/3xx/4xx/5xx)
  - CPU and memory utilization (p95/p99)
  - request latency (p50/p95/p99) with sample counts
  - container cold-start count and startup latency (p50/p95/p99)

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
import os
import select
import subprocess
import sys
import time
from datetime import datetime, timezone, timedelta

from google.cloud import monitoring_v3
from google.cloud import bigquery


PROJECT = "spliceai-lookup-412920"
REGION = "us-central1"
SERVICES = ["liftover", "spliceai-37", "spliceai-38", "pangolin-37", "pangolin-38"]

# Cache populated by `~/.claude/skills/analyze-gcloud-costs/scripts/cost_analysis.py`.
# If absent, the cost section is skipped with a warning instead of failing.
BILLING_CACHE_DIR = os.path.expanduser("~/.cache/analyze-gcloud-costs")


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


def sample_count(client, metric, start, end, revisions=None):
    """Return {service: total_sample_count} for a DELTA+DISTRIBUTION metric over the window.

    Sums the per-bucket counts across the distribution so we know how many raw observations
    contributed to each percentile in `percentiles()` above.
    """
    span = max(60, int((end - start).total_seconds()))
    f = f'metric.type="{metric}"'
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
            per_series_aligner=monitoring_v3.Aggregation.Aligner.ALIGN_DELTA,
            cross_series_reducer=monitoring_v3.Aggregation.Reducer.REDUCE_SUM,
            group_by_fields=["resource.label.service_name"],
        ),
    }):
        svc = ts.resource.labels.get("service_name", "?")
        for p in ts.points:
            out[svc] = out.get(svc, 0) + int(p.value.distribution_value.count)
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


def snapshot(client, args, bq_client=None, billing_table=None):
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

    print("=== Response codes (3xx and 404s ignored as probe/redirect noise) ===")
    codes = request_counts(client, start, now, revisions=prod_revs)
    for svc in SERVICES:
        all_codes = codes.get(svc, {})
        by_code = {code: c for code, c in all_codes.items() if code != 404 and code // 100 != 3}
        classes = collections.Counter()
        for code, c in by_code.items():
            classes[f"{code // 100}xx"] += c
        total = sum(classes.values())
        rate = classes["5xx"] / total * 100 if total else 0
        ignored = all_codes.get(404, 0)
        ignored_note = f"; +{ignored} 404s ignored" if ignored else ""
        print(f"  {svc:<14}  2xx={classes['2xx']:<5} "
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

    print("=== CPU / Memory utilization ===")
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

    lat_metric = "run.googleapis.com/request_latencies"
    lat = percentiles(client, lat_metric, start, now, revisions=prod_revs)
    lat_n = sample_count(client, lat_metric, start, now, revisions=prod_revs)
    if args.baseline_end:
        baseline_end = parse_iso(args.baseline_end)
        baseline_start = baseline_end - timedelta(days=args.baseline_days)
        # Baseline window predates current revisions; query unfiltered to capture pre-deploy traffic.
        baseline_lat = percentiles(client, lat_metric, baseline_start, baseline_end)
        print(f"=== Latency (p50/p95/p99 in s) — vs {args.baseline_days:g}d baseline ending {args.baseline_end} ===")
    else:
        baseline_lat = None
        print("=== Latency (p50/p95/p99 in s) ===")

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

    # Container startup latency distribution — count = number of cold starts in window,
    # percentiles = how long each new instance took to become ready to serve requests.
    print("=== Container startup (cold starts; p50/p95/p99 in s) ===")
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

    print(f"=== Project cost over the last {args.cost_days:g} days (net of credits) ===")
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
                        help="Length of each snapshot's window in hours (default: 2)")
    parser.add_argument("--baseline-end",
                        help="Optional ISO-8601 (YYYY-MM-DDTHH:MM:SSZ) end of a baseline window for latency comparison.")
    parser.add_argument("--baseline-days", type=float, default=7.0,
                        help="Length of the baseline window in days (default: 7)")
    parser.add_argument("--all-revisions", action="store_true",
                        help="Aggregate metrics across all revisions of each service "
                             "(default: filter to the revision currently serving 100%% production traffic, "
                             "so dev/test traffic against `dev---*` URLs doesn't contaminate metrics).")
    parser.add_argument("--cost-days", type=float, default=14.0,
                        help="Number of days of daily-cost history to chart (default: 14).")
    args = parser.parse_args()

    client = monitoring_v3.MetricServiceClient()
    billing_table = discover_billing_export()
    # Reuse a single client across iterations so we don't repeat ADC + project discovery
    # every 30 minutes. The BQ project hosts the billing export (not PROJECT itself).
    bq_client = bigquery.Client(project=billing_table.split(".")[0]) if billing_table else None
    while True:
        print("Processing...")
        try:
            snapshot(client, args, bq_client=bq_client, billing_table=billing_table)
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
