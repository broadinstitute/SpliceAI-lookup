
#### SpliceAI Lookup API

This API serves as the backend for [spliceailookup.broadinstitute.org](http://spliceailookup.broadinstitute.org).

The source code for the server and web UI is available @ [github.com/broadinstitute/SpliceAI-lookup](https://github.com/broadinstitute/SpliceAI-lookup) and is maintained by the [TGG](https://the-tgg.org/).   

For more details on **SpliceAI**, see [Jaganathan et al, Cell 2019 in press.](https://doi.org/10.1016/j.cell.2018.12.015) and [github.com/Illumina/SpliceAI](https://github.com/Illumina/SpliceAI).

As of February, 2023, this server also computes [Pangolin](https://github.com/tkzeng/Pangolin) scores.

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
1. Install and start a [redis](https://redis.io/) server. It's used to cache previously computed API server responses so that they don't have to be computed again.
1. Download reference fasta files: [hg19.fa](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) and [hg38.fa](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
1. Download [annotation files](https://spliceailookup-api.broadinstitute.org/annotations) into your local ./annotations directory.
1. Optionally download pre-computed scores .vcf.gz and .vcf.gz.tbi files from [Illumina Basespace](https://basespace.illumina.com/s/otSPW8hnhaZR)   
1. Start a SpliceAI API server on localhost port 8080. To modify server options, edit the `start_local_server.sh` script:

```
$ git clone git@github.com:broadinstitute/SpliceAI-lookup.git  # clone this repo  
$ cd SpliceAI-lookup  
$ python3 -m pip install -r requirements.txt  # install python dependencies  
$ ./start_local_server.sh  
```

The server uses ~1.5 Gb RAM per server thread.

