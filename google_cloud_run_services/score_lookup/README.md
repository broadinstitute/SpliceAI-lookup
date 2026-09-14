This folder contains the `score-lookup` [Cloud Run function](https://cloud.google.com/functions/docs/concepts/overview) used by the "Other scores" section of [spliceai-lookup.broadinstitute.org](https://spliceai-lookup.broadinstitute.org). For one variant it returns the AlphaMissense, PrimateAI-3D, PromoterAI, and AlphaGenome Variant Impact (AVI) scores from the tabix-indexed tables in `gs://spliceai-lookup-reference-data`.

The page used to query those tables directly from the browser, which downloaded several MiB of table indexes per session and was billed as internet data transfer. The function runs in the same region as the bucket and sends back only a few hundred bytes of JSON.

### API

```
GET https://score-lookup-xwkwwwxdwq-uc.a.run.app/?hg=38&chrom=8&pos=140300616&ref=T&alt=G
```

```json
{
  "variant": {"hg": "38", "chrom": "chr8", "pos": 140300616, "ref": "T", "alt": "G"},
  "scores": {"alphagenome_avi": {"rawScore": 1.59, "phredScore": 30.48}},
  "errors": []
}
```

* `hg` is `37` or `38`, `chrom` may include the `chr` prefix, and `ref` / `alt` are the variant's alleles.
* `scores` includes only the scores found for the variant: `alphamissense`, `primateai3d`, `promoterai`, and (hg38 only) `alphagenome_avi`.
* `errors` lists any table whose query failed. The other tables' scores are still returned, and the response isn't cached.
* Invalid parameters return a 400 with an `error` message.

### Development

* `npm install` then `npm test` runs the tests. The tests query the real tables unless `SKIP_NETWORK_TESTS=1` is set.
* `npm start` runs the function locally on port 8080.
* `python3 deploy.py --dev` deploys a no-traffic revision tagged `dev`, which the dev site uses. `python3 deploy.py --promote` then moves all traffic to it. `python3 deploy.py` deploys straight to production.
