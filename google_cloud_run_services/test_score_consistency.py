"""Use the cache to check if any scores have changed since the last time the scores were computed."""

import collections
import configargparse
import json
import os
import pandas as pd
import psycopg2
import requests
import tqdm

p = configargparse.ArgParser(default_config_files=["~/.spliceai_lookup_db_config"])
p.add_argument("--ip", required=True)
p.add_argument("--user", required=True)
p.add_argument("--password", required=True)
p.add_argument("--db", default="spliceai-lookup-db")
p.add_argument("-n", type=int, help="number of rows to query", default=1000)
args, _ = p.parse_known_args()


conn = psycopg2.connect(f"dbname='{args.db}' user='{args.user}' host='{args.ip}' password='{args.password}'")
df = pd.read_sql_query(f"SELECT key, value FROM cache LIMIT {args.n}", conn)
counter = collections.Counter()
for cache_key, cache_value in tqdm.tqdm(zip(df.key, df.value), total=len(df), unit=" variants", unit_scale=True):
	print("Processing", cache_key)
	data = json.loads(cache_value)
	hg = data["genomeVersion"]
	tool = data["source"].split(":")[0]
	distance = data["distance"]
	assert cache_key[-2:] in ("m1", "m0")
	mask = cache_key[-1]
	variant = data["variant"]
	cached_scores = data["scores"]

	counter[f"   {tool}"] += 1
	counter[f"  hg{hg}"] += 1
	counter[f" m{mask}"] += 1

	# get json response
	url = f"https://{tool}-{hg}-xwkwwwxdwq-uc.a.run.app/{tool}/?hg={hg}&distance={distance}&mask={mask}&variant={variant}&raw={variant}"
	#print(url)
	response_json = requests.get(f"{url}&force=1").json()
	response_scores = response_json["scores"]
	cached_scores = cached_scores[0]
	response_scores = response_scores[0]
	missing_keys = set()
	mismatched_values = set()
	for k, v in cached_scores.items():
		if k not in response_scores:
			missing_keys.add(k)
			continue
		if v != response_scores[k]:
			mismatched_values.add(k)
			continue

	if missing_keys:
		print(f"ERROR: {cache_key} response is missing keys: {missing_keys}. Response: {json.dumps(response_json, indent=1)}")
	if mismatched_values:
		print(f"ERROR: {cache_key} response has mismatched values for keys: {mismatched_values}")
		print(f"	Cache: {json.dumps(data, indent=1)}")
		print(f"	Response: {json.dumps(response_json, indent=1)}")

	#print("Done with", cache_key)
conn.close()

print(f"Done")

print("Stats:")
for key, value in sorted(counter.items()):
	print(f"{value:10,d}    {key}")

#%%
