"""Use the cache to check if any scores have changed since the last time the scores were computed."""

import collections
import configargparse
import datetime
import json
import os
import pandas as pd
import psycopg2
import re
import requests
import tqdm
import time

from contextlib import contextmanager

# Mirror server.py's get_pangolin_scores transcript priority + score-sum logic.
# Earlier versions of this script imported sai10k_select_transcript, but that
# function tie-breaks on the SpliceAI delta keys (DS_AG/DS_AL/DS_DG/DS_DL),
# which are absent from Pangolin cache rows — the tie-break sum was always 0
# and selection collapsed to "first transcript at the highest priority",
# disagreeing with the production response and causing spurious "Transcript
# ids don't match" skips that masked real consistency regressions.
PANGOLIN_TRANSCRIPT_PRIORITY_ORDER = {"MS": 4, "MP": 3, "C": 2, "N": 1}

def _select_pangolin_transcript(scores_list):
    if not scores_list:
        return None
    best = None
    best_priority = -1
    best_sum = -1.0
    for s in scores_list:
        priority = PANGOLIN_TRANSCRIPT_PRIORITY_ORDER.get(s.get("t_priority", "N"), 0)
        score_sum = abs(float(s.get("DS_SL", 0))) + abs(float(s.get("DS_SG", 0)))
        if priority > best_priority or (priority == best_priority and score_sum > best_sum):
            best = s
            best_priority = priority
            best_sum = score_sum
    return best

@contextmanager
def get_db_cursor(conn):
    """Get a database cursor"""
    cursor = conn.cursor()
    try:
        yield cursor
        conn.commit()
    finally:
        cursor.close()

def run_sql(conn, sql_query, *params):
    with get_db_cursor(conn) as cursor:
        cursor.execute(sql_query, *params)
        try:
            results = cursor.fetchall()
        except psycopg2.ProgrammingError:
            # No result set (e.g. from a DELETE/INSERT/UPDATE); caller just needs [].
            results = []
    return results

def main():
    # Body was previously at module top level, which meant `python3 -m unittest
    # test_score_consistency` (or any unittest discovery run) would issue
    # production-DB queries and DELETE rows during import. Wrap in a main()
    # gated by __main__ so the file is safe to import.
    p = configargparse.ArgParser(default_config_files=["~/.spliceai_lookup_db_config"])
    p.add_argument("--ip", required=True)
    p.add_argument("--user", required=True)
    p.add_argument("--password", required=True)
    p.add_argument("--db", default="spliceai-lookup-db")
    p.add_argument("-n", type=int, help="number of rows to query", default=1000)
    p.add_argument("-p", "--show-progress-bar", action="store_true")
    p.add_argument(
        "--delete-own-logs",
        action="store_true",
        help="After each API call, DELETE log rows this script generated since "
             "it started running (matched on the caller's public IP from "
             "checkip.dyndns.com AND logtime >= the script's start time). "
             "Off by default because an unexpected IP match — e.g. if checkip "
             "returns a shared/corporate NAT address — would delete unrelated "
             "production log rows from the same run window.",
    )
    args, _ = p.parse_known_args()

    # Capture start time BEFORE any HTTP calls so the DELETE bound is strictly
    # earlier than any log row this script could have produced.
    run_start_ts = datetime.datetime.now()

    myip_match = None
    if args.delete_own_logs:
        myip = requests.get("http://checkip.dyndns.com", timeout=10).text
        myip_match = re.search(r'Address: (\d+\.\d+\.\d+\.\d+)', myip)
        if myip_match:
            myip_match = myip_match.group(1)

    days_ago = 30
    conn = psycopg2.connect(f"dbname='{args.db}' user='{args.user}' host='{args.ip}' password='{args.password}'")
    #query = f"SELECT key, value, accessed FROM cache WHERE accessed < now() - INTERVAL '{days_ago} days' ORDER BY accessed ASC"
    #query = f"SELECT key, value, accessed FROM cache WHERE key LIKE 'pangolin%hg38%' AND accessed > now() - INTERVAL '{days_ago} days' ORDER BY accessed ASC"
    query = f"SELECT key, value, accessed FROM cache WHERE key LIKE 'pangolin%hg38%' ORDER BY accessed ASC"
    df = pd.read_sql_query(query, conn)
    print(f"Retrieved {len(df):,d} records from cache that were last accessed less than {days_ago} days ago.")
    if args.n:
        keep_every_kth_record = len(df)//args.n
        if keep_every_kth_record > 1:
            df = df[df.index % keep_every_kth_record == 0]
            print(f"Kept {len(df):,d} records after applying -n {args.n} arg")

    counter = collections.Counter()
    iterator = zip(df.key, df.value, df.accessed)
    if args.show_progress_bar:
        iterator = tqdm.tqdm(iterator, total=len(df), unit=" variants", unit_scale=True)

    # server.py appends "__{DEPLOYMENT}" to the key of every non-prod deployment, because dev
    # revisions share this database with production (see get_splicing_scores_cache_key). The
    # endpoint is chosen once, from SPLICEAI_API_ENV, so a run that did not match the two would
    # recompute dev-namespaced entries against production, or production entries against dev,
    # and report the difference between two deployments as a scoring inconsistency.
    want_dev_namespace = os.environ.get("SPLICEAI_API_ENV") == "dev"

    for i, (cache_key, cache_value, last_accessed) in enumerate(iterator):
        if cache_key.endswith("__dev") != want_dev_namespace:
            counter[f"skipped: not the {'dev' if want_dev_namespace else 'prod'} namespace"] += 1
            continue
        print(f"{i+1:3,d}: Processing", cache_key, "which was last accessed on", last_accessed)
        data = json.loads(cache_value)

        if not data.get("scores"):
            print("ERROR: No scores found in cached value. Skipping...")
            continue

        tool = data["source"].split(":")[0]
        hg = data["genomeVersion"]
        distance = data["distance"]
        # Versioned components follow the mask in the key (the cache version, and depending on the
        # tool and the deployment some of sai10k / model / gencode / dev), so cache_key[-1] isn't
        # the mask digit. Pull the mask out of the `__m{0,1}__` field by name instead.
        mask_match = re.search(r"__m([01])__", cache_key)
        assert mask_match, f"Could not parse mask from cache key: {cache_key}"
        mask = mask_match.group(1)
        # The gene set is part of the key too, and it has to be sent back explicitly: each
        # service now pins one gene set, and a request that omits bc is answered with whatever
        # the contacted service pins (server.py's GENE_SET default). Without this, every
        # __comprehensive__ row would be recomputed by the basic service and reported as a
        # mismatch. The server 307-redirects bc=comprehensive to the -comprehensive sibling,
        # and requests follows that, so the un-suffixed base_url below still reaches the right
        # service.
        # Anchored to the mask field so it cannot match anything inside the variant, and the
        # trailing "__" is optional because the field is last in a legacy Pangolin key: before
        # the CACHE_VERSION component was added, the key ended at the gene set for every tool
        # whose suffix was empty, which is Pangolin (the sai10k suffix is spliceai-only). Those
        # rows are exactly the ones that must not silently fall back to basic -- a legacy
        # "__comprehensive" entry was computed with the comprehensive annotations and has to be
        # rechecked against them.
        gene_set_match = re.search(r"__m[01]__(basic|comprehensive)(?:__|$)", cache_key)
        gene_set = gene_set_match.group(1) if gene_set_match else None
        variant = data["variant"]

        # get json response
        # time requests
        start_time = time.time()
        # Hostname slug is GCP-assigned and may change on service recreate;
        # allow override via SPLICEAI_API_URL_TEMPLATE (same env var as
        # test_api_consistency.py). SPLICEAI_API_ENV=dev points at the
        # 'dev'-tagged revisions from build_and_deploy.py --dev.
        default_template = "https://{tool}-{hg}-xwkwwwxdwq-uc.a.run.app"
        if os.environ.get("SPLICEAI_API_ENV") == "dev":
            default_template = "https://dev---{tool}-{hg}-xwkwwwxdwq-uc.a.run.app"
        url_template = os.environ.get("SPLICEAI_API_URL_TEMPLATE", default_template)
        base_url = url_template.format(tool=tool, hg=hg)
        bc_arg = f"&bc={gene_set}" if gene_set else ""
        url = f"{base_url}/{tool}/?hg={hg}&distance={distance}&mask={mask}{bc_arg}&variant={variant}&raw={variant}"
        # print(url)
        try:
            response_json = requests.get(f"{url}&force=1").json()
        except Exception as e:
            print(f"ERROR: {e} when retrieving {url}  Skipping...")
            continue

        if not response_json.get("scores"):
            print(f"ERROR: {url} response doesn't contain scores: {response_json}. Skipping...")
            continue

        elapsed_time = time.time() - start_time
        # Pick the canonical transcript the same way get_pangolin_scores does
        # (MS > MP > C > N priority, tie-broken by sum of |DS_SL| + |DS_SG|).
        # Sorting alphabetically by t_id could select a non-canonical transcript
        # on one side and a different one on the other, producing spurious
        # mismatches and masking real ones.
        response_scores = _select_pangolin_transcript(response_json["scores"])
        cached_scores = _select_pangolin_transcript(data["scores"])
        if response_scores is None or cached_scores is None:
            print(f"ERROR: no transcripts in response or cache for {variant}. Skipping...")
            continue

        if not response_scores.get("t_id") or response_scores.get("t_id") != cached_scores.get("t_id"):
            print("Transcript ids don't match:", response_scores.get("t_id"), "vs", cached_scores.get("t_id"),
                  ". Skipping...")
            continue

        if not response_scores.get("g_id") or response_scores.get("g_id") != cached_scores.get("g_id"):
            print("Gene ids don't match:", response_scores.get("g_id"), "vs", cached_scores.get("g_id"),
                  ". Skipping...")
            continue

        counter[f"   {tool}"] += 1
        counter[f"  hg{hg}"] += 1
        counter[f" m{mask}"] += 1

        missing_keys = set()
        mismatched_values = set()
        values_to_print = set()
        for k, v1 in cached_scores.items():
            if k in ("t_refseq_ids", "t_id", "g_id", "g_name"):
                # differences in gene ids are not important
                continue

            # ALL_NON_ZERO_SCORES is written to the cache but stripped from the HTTP response
            # (server.py keeps it in the dict for the cache, then pops it alongside NAME and
            # STRAND before serializing), so comparing it here reports it missing for every
            # single entry -- which also means a genuine missing key could never be spotted.
            # EXON_STARTS / EXON_ENDS are lists, and the values below go into a set, which a
            # list cannot enter: without this the first Pangolin entry carrying exon
            # coordinates aborts the whole run with TypeError: unhashable type: 'list'.
            if k in ("ALL_NON_ZERO_SCORES", "EXON_STARTS", "EXON_ENDS"):
                continue

            if k not in response_scores:
                missing_keys.add(k)
                continue

            v2 = response_scores[k]
            try:
                diff = float(v1) - float(v2)
            except:
                diff = "?"

            values_to_print.add((k, v1, v2, diff))
            if v1 != v2:
                mismatched_values.add((k, v1, v2, diff))
                continue

        if missing_keys:
            print(f"ERROR: {cache_key} which was last accessed on {last_accessed} is missing keys: {missing_keys}. Response: {json.dumps(response_json, indent=1)}")

        if mismatched_values:
            counter["ERROR: mismatched_values"] += 1
            print(f"ERROR: {cache_key} which was last accessed on {last_accessed} has mismatched values for keys: "
                  f"{', '.join(sorted([t[0] for t in mismatched_values]))} "
                  f"with max delta_score_diff="
                  f"{max([abs(t[3]) for t in mismatched_values if t[0].startswith('DS')] or [None])} "
                  f"and max raw_score_diff="
                  f"{max([abs(t[3]) for t in mismatched_values if t[0].startswith('S')] or [None])} ")

            for k, v1, v2, diff in sorted(values_to_print):
                print(f"    {k}:  {v1}  vs  {v2}   diff: {diff}")

            #print(f"	Cache: {json.dumps(data, indent=1)}")
            #print(f"	Response: {json.dumps(response_json, indent=1)}")

            df = pd.read_sql_query("SELECT * FROM log WHERE variant=%s", conn, params=(variant,))
            print(f"        Log:")
            print(df.to_string(index=False))

        print(f"{i+1:3,d}: Done with", cache_key, f"elapsed_time={elapsed_time:.1f}s")

    # One DELETE after the loop instead of one per iteration. The
    # `logtime >= run_start_ts` clause already scopes the delete to this run, so
    # repeating it inside the loop just issues N redundant DELETEs against
    # production for no benefit.
    if myip_match:
        print(f"Deleting logs for ip {myip_match} since {run_start_ts}")
        run_sql(
            conn,
            "DELETE FROM log WHERE ip=%s AND logtime >= %s",
            (myip_match, run_start_ts),
        )

    conn.close()

    print(f"Done")

    print("Stats:")
    for key, value in sorted(counter.items()):
        print(f"{value:10,d}    {key}")


if __name__ == "__main__":
    main()
