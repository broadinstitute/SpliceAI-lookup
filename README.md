This repo contains: 
- client-side code for [spliceailookup.broadinstitute.org](https://spliceailookup.broadinstitute.org/) - contained within the [index.html](index.html) file and hosted via GitHub Pages.
- server-side code for SpliceAI and Pangolin REST APIs - contained within the [google_cloud_run_services/](google_cloud_run_services/) subdirectory and hosted on Google Cloud Run. 

---
#### API Examples:

*[/spliceai/?hg=38&distance=50&variant=chr8-140300616-T-G](http://spliceailookup-api.broadinstitute.org/spliceai/?hg=38&variant=chr8-140300616-T-G)*
  
Get SpliceAI scores for the given variant.   

- **variant** (required) a variant in the format "chrom-pos-ref-alt"
- **hg** (required) can be 37 or 38
- **distance** (optional) distance parameter of SpliceAI model (default: 50)  
- **mask** (optional) can be 0 which means raw scores or 1 which means masked scores (default: 0). 
Splicing changes corresponding to strengthening annotated splice sites and weakening unannotated splice sites are typically much less pathogenic than weakening annotated splice sites and
strengthening unannotated splice sites. When this parameter is = 1 (masked), the delta scores of such splicing changes are set to 0. SpliceAI developers recommend using raw (0) for alternative splicing analysis and masked (1) for variant interpretation.  

*[/pangolin/?hg=38&distance=50&variant=chr8-140300616-T-G](http://spliceailookup-api.broadinstitute.org/pangolin/?hg=38&variant=chr8-140300616-T-G)*

Get Pangolin scores for the given variant.

- **variant** (required) a variant in the format "chrom-pos-ref-alt"
- **hg** (required) can be 37 or 38
- **distance** (optional) distance parameter of SpliceAI model (default: 50)
- **mask** (optional) can be 0 which means raw scores or 1 which means masked scores (default: 0).
  Splicing changes corresponding to strengthening annotated splice sites and weakening unannotated splice sites are typically much less pathogenic than weakening annotated splice sites and
  strengthening unannotated splice sites. When this parameter is = 1 (masked), the delta scores of such splicing changes are set to 0. SpliceAI developers recommend using raw (0) for alternative splicing analysis and masked (1) for variant interpretation.

---
#### Local Install

The steps below describe how to install a SpliceAI API server on your local infrastructure.
The details will vary depending on your OS, etc. If you run into issues, please submit them
to the [issue tracker](https://github.com/broadinstitute/SpliceAI-lookup/issues).

1. Install pytorch as described in the [Pangolin installation docs](https://github.com/tkzeng/Pangolin#installation)
1. Install the modified versions of SpliceAI and Pangolin from [https://github.com/bw2/SpliceAI](https://github.com/bw2/SpliceAI) and [https://github.com/bw2/Pangolin](https://github.com/bw2/Pangolin)
1. Install and start a [redis](https://redis.io/) server. It's used to cache previously computed API server responses so that they don't have to be computed again.
1. Download reference fasta files: [hg19.fa](https://storage.cloud.google.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta) and [hg38.fa](https://storage.cloud.google.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta)
1. Download [annotation files](https://spliceailookup-api.broadinstitute.org/annotations) into your local ./annotations directory, or generate them yourself using the steps in the [annotations README](https://github.com/broadinstitute/SpliceAI-lookup/blob/master/annotations/README.md).
1. Start a SpliceAI API server on localhost port 8080. To modify server options, edit the `start_local_server.sh` script:

```
$ git clone git@github.com:broadinstitute/SpliceAI-lookup.git  # clone this repo  
$ cd SpliceAI-lookup  
$ python3 -m pip install -r requirements.txt  # install python dependencies  
$ ./start_local_server.sh  
```

The server uses ~1.5 Gb RAM per server thread.

---
#### For Developers

The [spliceailookup.broadinstitute.org](https://spliceailookup.broadinstitute.org) front-end is contained within [index.html](index.html). It uses ES6 javascript with [Semantic UI](https://semantic-ui.com) and [jQuery](https://en.wikipedia.org/wiki/JQuery). Also, it uses a [custom version of igv.js](https://github.com/bw2/igv.js) that includes new track types for visualizing the SpliceAI & Pangolin scores. The original server-side code is in [server.py](server.py) and uses the [Flask](https://flask.palletsprojects.com/en/3.0.x) library. It is designed to run on a plain Linux or MacOS machine. The new server-side code is in the [google_cloud_run_services/](google_cloud_run_services/) subdirectory and includes Dockerfiles and scripts for deploying SpliceAI and Pangolin API services to Google Cloud Run. 


