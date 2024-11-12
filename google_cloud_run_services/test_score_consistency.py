"""Use the cache to check if any scores have changed since the last time the scores were computed."""

import collections
import configargparse
import json
import os
import pandas as pd
import psycopg2
import requests
import tqdm
import time

p = configargparse.ArgParser(default_config_files=["~/.spliceai_lookup_db_config"])
p.add_argument("--ip", required=True)
p.add_argument("--user", required=True)
p.add_argument("--password", required=True)
p.add_argument("--db", default="spliceai-lookup-db")
p.add_argument("-n", type=int, help="number of rows to query", default=1000)
p.add_argument("-p", "--show-progress-bar", action="store_true")
args, _ = p.parse_known_args()

days_ago = 3
conn = psycopg2.connect(f"dbname='{args.db}' user='{args.user}' host='{args.ip}' password='{args.password}'")
#query = f"SELECT key, value, accessed FROM cache WHERE accessed < now() - INTERVAL '{days_ago} days' ORDER BY accessed ASC"
query = f"SELECT key, value, accessed FROM cache WHERE key LIKE 'pangolin%hg38%' AND accessed > now() - INTERVAL '{days_ago} months' ORDER BY accessed ASC"
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

for i, (cache_key, cache_value, last_accessed) in enumerate(iterator):
    print(f"{i+1:3,d}: Processing", cache_key, "which was last accessed on", last_accessed)
    data = json.loads(cache_value)

    if not data.get("scores"):
        print("ERROR: No scores found in cached value. Skipping...")
        continue

    tool = data["source"].split(":")[0]
    hg = data["genomeVersion"]
    distance = data["distance"]
    cache_key = cache_key.replace("__basic", "").replace("__comprehensive", "")
    assert cache_key[-2:] in ("m1", "m0")
    mask = cache_key[-1]
    variant = data["variant"]

    # get json response
    # time requests
    start_time = time.time()
    url = f"https://{tool}-{hg}-xwkwwwxdwq-uc.a.run.app/{tool}/?hg={hg}&distance={distance}&mask={mask}&variant={variant}&raw={variant}"
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
    response_json["scores"] = list(sorted(response_json["scores"], key=lambda s: s.get("t_id")))
    response_scores = response_json["scores"][0]

    data["scores"] = list(sorted(data["scores"], key=lambda s: s.get("t_id")))
    cached_scores = data["scores"][0]

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

        df = pd.read_sql_query(f"SELECT * FROM log WHERE variant='{variant}'", conn)
        print(f"        Log:")
        print(df.to_string(index=False))
                
    print(f"{i+1:3,d}: Done with", cache_key, f"elapsed_time={elapsed_time:.1f}s")

conn.close()

print(f"Done")

print("Stats:")
for key, value in sorted(counter.items()):
    print(f"{value:10,d}    {key}")

#%%
