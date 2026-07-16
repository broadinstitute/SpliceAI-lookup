This repo contains: 
- client-side code for [spliceailookup.broadinstitute.org](https://spliceailookup.broadinstitute.org/) - contained within the [index.html](index.html) file and hosted via GitHub Pages.
- server-side code for SpliceAI and Pangolin REST APIs - contained within the [google_cloud_run_services/](google_cloud_run_services/) subdirectory and hosted on Google Cloud Run. 
- the original Flask web app that previously powered spliceailookup.broadinstitute.org - contained within the [original_flask_app/](original_flask_app/) subdirectory. 

---

#### SpliceAI, Pangolin APIs


<b>NOTE:</b> These APIs are intended for interactive use only, and do not support more than several requests per user per minute. More frequent queries will trigger a "rate limit" error in the response. To process large batches of variants, please set up and query your own local instance of the API server. This is easy to do using the publicly available docker images (see below for details). Alternatively, you can intall and run the underlying SpliceAI and/or Pangolin models directly on your local infrastructure. Their source code is available @ [https://github.com/bw2/SpliceAI](https://github.com/bw2/SpliceAI) and [https://github.com/bw2/Pangolin](https://github.com/bw2/Pangolin). <br />
<br />

The SpliceAI and Pangolin APIs have different base urls for different genome versions:

`https://spliceai-37-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=37&variant=` - SpliceAI for variants on GRCh37<br />
`https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=38&variant=` - SpliceAI for variants on GRCh38<br />
`https://pangolin-37-xwkwwwxdwq-uc.a.run.app/pangolin/?hg=37&variant=` - Pangolin for variants on GRCh37<br />
`https://pangolin-38-xwkwwwxdwq-uc.a.run.app/pangolin/?hg=38&variant=` - Pangolin for variants on GRCh38 <br />

To query the API, append your variant of interest in `chrom-pos-ref-alt` format to the appropriate base url above.

For example, to get SpliceAI scores for `chr8-140300616-T-G`:<br>

*[https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=38&variant=chr8-140300616-T-G](https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=38&variant=chr8-140300616-T-G)*
  
or to get Pangolin scores while also setting the `distance` and `mask` parameters:<br>

*[https://pangolin-38-xwkwwwxdwq-uc.a.run.app/pangolin/?hg=38&variant=chr8-140300616-T-G&distance=1000&mask=1](https://pangolin-38-xwkwwwxdwq-uc.a.run.app/pangolin/?hg=38&variant=chr8-140300616-T-G&distance=1000&mask=1)*

#### API parameters

Parameter descriptions:  

- **variant** (required) a variant in the format "chrom-pos-ref-alt"  
- **hg** (required) can be 37 or 38  
- **distance** (optional) distance parameter of SpliceAI model (default: 50)   
- **mask** (optional) can be 0 which means raw scores or 1 which means masked scores (default: 0). 
Splicing changes corresponding to strengthening annotated splice sites and weakening unannotated splice sites are typically much less pathogenic than weakening annotated splice sites and
strengthening unannotated splice sites. When this parameter is = 1 (masked), the delta scores of such splicing changes are set to 0. SpliceAI developers recommend using raw (0) for alternative splicing analysis and masked (1) for variant interpretation.  


---
#### Running Your Own Local API Server

If you have [docker](https://docs.docker.com/engine/install/) installed, you can easily start your own SpliceAI-lookup API server by running one of these commands (depending on which model and genome version you want to query):

```
docker run -p 8080:8080 docker.io/weisburd/spliceai-38:latest
docker run -p 8080:8080 docker.io/weisburd/spliceai-37:latest
docker run -p 8080:8080 docker.io/weisburd/pangolin-38:latest
docker run -p 8080:8080 docker.io/weisburd/pangolin-37:latest
```

On startup the container loads the model and prints some TensorFlow / model-loading warnings (e.g. `WARNING:absl:No training configuration found...`); these can be ignored. Once it is listening on port 8080, you can query it. For example, if you started the `spliceai-38` image, you can open http://localhost:8080/spliceai/?hg=38&variant=chr8-140300616-T-G in your browser. Each per-genome/per-tool image only answers requests for its own tool and `hg` value. 

Optional environment variables (pass with `docker run -e NAME=value ...`):
- `DISABLE_RATE_LIMIT=1` — explicitly turn off per-IP rate limiting. Only relevant if you connect your own database (see below); without a database, rate limiting is already disabled.
- `DATABASE_ENABLED=1` / `0` — force the database on or off. Defaults to on when `DB_PASSWORD` is set and off otherwise.
- `DB_HOST`, `DB_PORT`, `DB_NAME`, `DB_USER`, `DB_PASSWORD` — connection settings for an optional PostgreSQL database (see below).


##### Optionally attaching your own PostgreSQL database

A database is not required, but attaching one to a local instance adds response caching, so repeated queries for the same variant return instantly. To attach a local PostgreSQL:

1. Create an empty database, e.g. `createdb spliceai-lookup-db`.
2. Point the container at it. From a Docker container, the host's PostgreSQL is reachable at `host.docker.internal` (on Docker Desktop for Mac/Windows; on Linux add `--add-host=host.docker.internal:host-gateway`):
   ```
   docker run -p 8080:8080 \
     -e DATABASE_ENABLED=1 \
     -e DB_HOST=host.docker.internal \
     -e DB_PORT=5432 \
     -e DB_NAME=spliceai-lookup-db \
     -e DB_USER=postgres \
     -e DB_PASSWORD=yourpassword \
     docker.io/weisburd/spliceai-38:latest
   ```
   If your PostgreSQL uses passwordless (`trust`) auth, omit `-e DB_PASSWORD` and keep `-e DATABASE_ENABLED=1`.

The server creates the tables it needs automatically on the first request (`cache`, `log`, `restricted_ips`, `whitelist_ips`) — no manual schema setup is required.

The `transcripts_hg37`/`transcripts_hg38` tables used for SAI-10k transcript-structure enrichment are *not* created automatically; without them SAI-10k falls back to the bundled annotations. To populate them, see the `update_transcript_tables` command in [build_and_deploy.py](https://github.com/broadinstitute/SpliceAI-lookup/blob/master/google_cloud_run_services/build_and_deploy.py).

If you would like to run your own API instance on Google Cloud instead of locally, see the [build_and_deploy.py](https://github.com/broadinstitute/SpliceAI-lookup/blob/master/google_cloud_run_services/build_and_deploy.py#L224-L238) script which we use to deploy and update the SpliceAI-lookup API on [Google Cloud Run](https://cloud.google.com/run?hl=en). Submit a GitHub issue if you have any questions.

---
#### For Developers

For how to make, test, and/or contribute changes, see [Developer Guidelines](CONTRIBUTING.md).


