
**SpliceAI Lookup API** 

This API serves as the backend for [spliceailookup.broadinstitute.org](http://spliceailookup.broadinstitute.org). 
It works by using precomputed scores provided by Illumina for SNVs and small InDels, and 
by running the [SpliceAI model](https://github.com/Illumina/SpliceAI) for all other queries - such as larger InDels or non-default **distance** values.  

The server and web UI code is available @ [github.com/broadinstitute/SpliceAI-lookup](https://github.com/broadinstitute/SpliceAI-lookup) and is maintained by the [TGG](https://the-tgg.org/).   

For more details on **SpliceAI**, see [Jaganathan et al, Cell 2019 in press.](https://doi.org/10.1016/j.cell.2018.12.015) and [github.com/Illumina/SpliceAI](https://github.com/Illumina/SpliceAI).

---
**API Example:**

*[/spliceai/?hg=38&distance=50&variant=chr8-140300616-T-G](http://34.71.77.14/spliceai/?hg=38&variant=chr8-140300616-T-G)*
  
Returns SpliceAI scores for the given variant.   

- **variant** (required) a variant in the format "chrom-pos-ref-alt"
- **hg** (required) can be 37 or 38
- **distance** (optional) distance parameter of SpliceAI model (default: 50)    
- **mask** (optional) can be 0 which means raw scores or 1 which means masked scores (default: 0). 
Splicing changes corresponding to strengthening annotated splice sites and weakening unannotated splice sites are typically much less pathogenic than weakening annotated splice sites and
strengthening unannotated splice sites. When this parameter is = 1 (masked), the delta scores of such splicing changes are set to 0. SpliceAI developers recommend using raw (0) for alternative splicing analysis and masked (1) for variant interpretation.    
