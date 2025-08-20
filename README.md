This repo contains: 
- client-side code for [spliceailookup.broadinstitute.org](https://spliceailookup.broadinstitute.org/) - contained within the [index.html](index.html) file and hosted via GitHub Pages.
- server-side code for SpliceAI and Pangolin REST APIs - contained within the [google_cloud_run_services/](google_cloud_run_services/) subdirectory and hosted on Google Cloud Run. 

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
When it starts, it will print:  
```
 * Serving Flask app 'server'
 * Debug mode: on
```   

Let's say you ran the `spliceai-38` instance. You should then be able to query it by, for example, opening http://localhost:8080/spliceai/?hg=38&variant=chr8-140300616-T-G in your browser.
The docker container will initially print:   
```
ERROR: Unable to connect to SQL database...
WARNING:absl:No training configuration found...
WARNING:tensorflow:...
```
but these messages can be ignored, and subsequent queries will run faster.


If you would like to run your own API instance on Google Cloud instead of locally, see the [build_and_deploy.py](https://github.com/broadinstitute/SpliceAI-lookup/blob/master/google_cloud_run_services/build_and_deploy.py#L224-L238) script which we use to deploy and update the SpliceAI-lookup API on [Google Cloud Run](https://cloud.google.com/run?hl=en). Submit a GitHub issue if you have any questions.

---
#### Code Overview For Developers

The [spliceailookup.broadinstitute.org](https://spliceailookup.broadinstitute.org) front-end is contained within [index.html](index.html). It uses ES6 javascript with [Semantic UI](https://semantic-ui.com) and [jQuery](https://en.wikipedia.org/wiki/JQuery). Also, it uses a [custom version of igv.js](https://github.com/bw2/igv.js) that includes new track types for visualizing the SpliceAI & Pangolin scores. The new server-side code is in the [google_cloud_run_services/](google_cloud_run_services/) subdirectory and includes Dockerfiles for building API server images, as well as the [build_and_deploy.py](https://github.com/broadinstitute/SpliceAI-lookup/blob/master/google_cloud_run_services/build_and_deploy.py#L224-L238) script for deploying SpliceAI and Pangolin API services to [Google Cloud Run](https://cloud.google.com/run?hl=en). 
The API server logic is in [google_cloud_run_services/server.py](https://github.com/broadinstitute/SpliceAI-lookup/blob/master/google_cloud_run_services/server.py) and uses the [Flask](https://flask.palletsprojects.com/en/3.0.x) library.


