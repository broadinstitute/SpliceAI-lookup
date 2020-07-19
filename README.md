
**This server provides the following APIs:**

----
*[/spliceai/?hg=38&distance=50&variants=chr8-140300615-C-G,chr8-140300616-T-G](/spliceai/?hg=38&variants=chr8-140300615-C-G,chr8-140300616-T-G)*
  
Returns spliceai scores for the given variants.  
The reponse is a json list that's the same length as the input "variants" list and 
contains the SpliceAI scores or error message for each variant.  
NOTE: The server takes ~ 2 seconds per variant.

- **hg** (required) can be 37 or 38
- **variants** (required) a comma-spearated list of one or variants in the format "chrom-pos-ref-alt" or "chrom:pos ref&gt;alt" or "chrom pos ref alt"
- **distance** (optional) distance parameter of SpliceAI model (default: 50)    


TODO:
- move variant parsing, ucsc links, etc. to client
- add firestore tables for precomputed variants
- cloud function
