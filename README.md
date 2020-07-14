
**This server provides the following APIs:**

----
*[/spliceai/?hg=38&variants=chr8-140300615-T-G,chr8-140300615-C-G](/spliceai/?hg=38&variants=chr8-140300615-T-G,chr8-140300615-C-G)*
  
Returns spliceai scores for the given variants.  
The reponse is a json list that's the same length as the input "variants" list and 
contains the SpliceAI scores or error message for each variant.  
NOTE: The server takes ~ 2 seconds per variant.

- **hg** can be 37 or 38    
- **variants** can have the format "chrom-pos-ref-alt" or "chrom:pos ref&gt;alt" or "chrom pos ref alt" <br />
