#%%
import argparse
import gzip
import os
import pandas as pd
import re

from collections import defaultdict

#%%

test_args = ["./annotations/gencode.v36lift37.annotation.gtf.gz"]
test_args = ["./annotations/gencode.v36.annotation.gtf.gz"]
test_args = None

#%%

p = argparse.ArgumentParser(description="""This script takes a Gencode .gtf.gz file
    and outputs an annotation file which can be passed to SpliceAI instead of 
    the default SpliceAI annotations which are still on Gencode v24. 
""")

p.add_argument("gtf_gz_path", help="Path of gene annotations file in GTF format")
args = p.parse_args(test_args)

print(f"Parsing {args.gtf_gz_path}")

# NOTE: the UCSC .bed file is in a format that would have been easier to use here, but even the UCSC comprehensive .bed
# file is missing some of the genes and transcripts from the comprehensive Gencode .gtf file, and additionally contains
# some unnecessary extra genes from super-contigs (which are present only in the Gencode "all" .gtf superset file
# called gencode.v36.chr_patch_hapl_scaff.annotation.gtf.gz)

#%%

def parse_gencode_file(gencode_gtf_path):
    with gzip.open(gencode_gtf_path, "rt") as gencode_gtf:
        for line in gencode_gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "exon":
                continue

            annotation_source = fields[1]
            chrom = fields[0]
            start_1based = int(fields[3])
            end_1based = int(fields[4])

            # don't emit the same exon coords more than once
            exon_coords = (chrom, start_1based, end_1based)

            meta_fields = {}
            for meta_field in fields[8].strip("; ").split(";"):
                key, value = meta_field.strip().replace('"', '').split()
                meta_fields[key] = value

            strand = fields[6]

            yield {
                "chrom": chrom,
                "start_1based": start_1based,
                "end_1based": end_1based,
                "annotation_source": annotation_source,
                "strand": strand,
                "gene_id": meta_fields["gene_id"].split(".")[0],
                "transcript_id": meta_fields["transcript_id"].split(".")[0],
                "gene_name": meta_fields["gene_name"],
                "gene_type": meta_fields["transcript_type"],
                "transcript_type":  meta_fields["transcript_type"],
            }


#%%

# aggregate gtf exon records into buckets keyed by (chrom, gene name, strand)
all_exons = defaultdict(set)
for record in parse_gencode_file(args.gtf_gz_path):
    key = (record["chrom"], record["gene_name"], record["strand"])
    all_exons[key].add((int(record['start_1based']), int(record['end_1based'])))

#%%

# reformat the aggregated records into a list which can be turned into a pandas DataFrame
output_records = []
for key in sorted(all_exons.keys()):
    exons_set = all_exons[key]

    chrom, gene_name, strand = key

    tx_start_0based = min([start_1based - 1 for start_1based, _ in exons_set])
    tx_end_1based = max([end_1based for _, end_1based in exons_set])
    exon_starts_0based = sorted([start_1based - 1 for start_1based, _ in exons_set])
    exon_ends_1based = sorted([end_1based for _, end_1based in exons_set])
    output_records.append({
        "#NAME": gene_name,
        "CHROM": chrom,
        "STRAND": strand,
        "TX_START": str(tx_start_0based),
        "TX_END": str(tx_end_1based),
        "EXON_START": ",".join([str(s) for s in exon_starts_0based]) + ",",
        "EXON_END": ",".join([str(s) for s in exon_ends_1based]) + ",",
    })


#%%

# generate output table
output_df = pd.DataFrame(output_records)
output_df = output_df[["#NAME", "CHROM", "STRAND", "TX_START", "TX_END", "EXON_START", "EXON_END"]]

#output_df[output_df['#NAME'] == "FGF16"]
output_path = re.sub(".gtf.gz$", "", os.path.basename(args.gtf_gz_path)) + ".txt"

output_df.to_csv(output_path, index=False, sep="\t")

print(f"Wrote {os.path.abspath(output_path)}")

#%%
