#%%

import argparse
from collections import defaultdict, Counter
import gzip
from intervaltree import IntervalTree, Interval
import os
import pandas as pd
from pprint import pprint
import re
import spliceai

from annotations.get_ensembl_db_info import get_canonical_transcripts, CURRENT_ENSEMBL_DATABASE

#%%
official_annotations_gene_names = set()
for genome_version in "37", "38":
    official_annotations_path = os.path.join(os.path.dirname(spliceai.__file__), f"annotations/grch{genome_version}.txt")
    with open(official_annotations_path, "rt") as f:
        official_annotations_gene_names.update({line.split("\t")[0] for line in f})


#%%

test_args = None
#test_args = ["./annotations/gencode.v24.annotation.gtf.gz"]
#test_args = ["./annotations/gencode.v36lift37.annotation.gtf.gz"]
#test_args = ["./annotations/gencode.v36.annotation.gtf.gz"]

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

TRANSCRIPT_TYPES_BY_PRIORITY = [
    'protein_coding',
    'translated_processed_pseudogene',
    'translated_unprocessed_pseudogene',
    'lncRNA',
    'lincRNA',
    'ribozyme',
    'polymorphic_pseudogene',
    'retained_intron',
    'sense_intronic',
    'non_stop_decay',
    'nonsense_mediated_decay',
    'sense_overlapping',
    'antisense',
    'processed_transcript',
    'rRNA_pseudogene',
    'unitary_pseudogene',
    'transcribed_processed_pseudogene',
    'transcribed_unitary_pseudogene',
    'transcribed_unprocessed_pseudogene',
    'unprocessed_pseudogene',
    'processed_pseudogene',
    'pseudogene',
    'vault_RNA',
    'rRNA',
    'snRNA',
    'sRNA',
    'scaRNA',
    'snoRNA',
    'scRNA',
    'Mt_tRNA',
    'Mt_rRNA',
    'IG_V_gene',
    'IG_C_gene',
    'IG_J_gene',
    'TR_C_gene',
    'TR_J_gene',
    'TR_V_gene',
    'TR_D_gene',
    'IG_D_gene',
    'IG_V_pseudogene',
    'TR_V_pseudogene',
    'IG_C_pseudogene',
    'TR_J_pseudogene',
    'IG_J_pseudogene',
    'IG_pseudogene',
    'TEC',
    'miRNA',
    'misc_RNA',
]


#%%

def parse_gencode_file(gencode_gtf_path):
    gene_type_summary_counter = Counter()
    gene_type_counter = Counter()
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

            meta_fields = {}
            for meta_field in fields[8].strip("; ").split(";"):
                key, value = meta_field.strip().replace('"', '').split()
                meta_fields[key] = value

            strand = fields[6]

            gene_type_counter[meta_fields["transcript_type"]] += 1
            if meta_fields["gene_name"] in official_annotations_gene_names:
                priority = "primary"
                gene_type_summary_counter["Including gene name from default SpliceAI annotation file"] += 1
            elif meta_fields["transcript_type"] == "protein_coding":
                priority = "primary"
                gene_type_summary_counter["Including new protein-coding gene"] += 1
            else:
                priority = "secondary"
                gene_type_summary_counter["Including new other gene"] += 1

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
                "priority": priority,
            }


    print("Summary:")
    pprint(gene_type_summary_counter)
    #print("Gene Types:")
    #pprint(gene_type_counter)
    #pprint(list(gene_type_counter.keys()))


#%%

# aggregate gtf exon records into buckets keyed by (chrom, gene name, strand)
all_exons_by_priority = {
    "primary": defaultdict(lambda: defaultdict(set)),
    "secondary": defaultdict(lambda: defaultdict(set)),
}

print(f"Getting canonical transcripts from {CURRENT_ENSEMBL_DATABASE}")
gene_id_to_canonical_transcript_id = get_canonical_transcripts()

for record in parse_gencode_file(args.gtf_gz_path):
    priority = record["priority"]
    transcript_type = record["transcript_type"]

    is_canonical_transcript = "no"
    if record["gene_id"] not in gene_id_to_canonical_transcript_id:
        print(f"WARNING: no canonical transcript for " + record["gene_id"])
    elif record["transcript_id"] == gene_id_to_canonical_transcript_id[record["gene_id"]]:
        is_canonical_transcript = "yes"

    name = "---".join([record["gene_name"], record["gene_id"], record["transcript_id"], is_canonical_transcript])
    key = (record["chrom"], name, record["strand"])

    all_exons_by_priority[priority][transcript_type][key].add((int(record['start_1based']), int(record['end_1based'])))


#%%

def transcript_type_order(transcript_type):
    try:
        return TRANSCRIPT_TYPES_BY_PRIORITY.index(transcript_type)
    except ValueError:
        return len(TRANSCRIPT_TYPES_BY_PRIORITY) + 1


# reformat the aggregated records into a list which can be turned into a pandas DataFrame
output_records = []
interval_trees = defaultdict(IntervalTree)
for priority in all_exons_by_priority:
    for transcript_type in sorted(all_exons_by_priority[priority].keys(), key=transcript_type_order):
        current_exon_sets = all_exons_by_priority[priority][transcript_type]
        for key in sorted(current_exon_sets.keys()):
            exons_set = current_exon_sets[key]

            chrom, gene_name, strand = key

            tx_start_0based = min([start_1based - 1 for start_1based, _ in exons_set])
            tx_end_1based = max([end_1based for _, end_1based in exons_set])

            # check for overlap with previously-added transcripts
            if priority != "primary" and interval_trees[chrom].overlaps(Interval(tx_start_0based, tx_end_1based)):
                # skip any secondary transcripts that overlap already-added primary transcripts
                overlapping_genes = sorted(set([i.data for i in interval_trees[chrom][tx_start_0based:tx_end_1based]]))
                print(f"Skipping {priority} {transcript_type} gene {gene_name} since it overlaps {len(overlapping_genes)} gene(s): " +
                    ", ".join(overlapping_genes[:5]) + ("..." if len(overlapping_genes) > 5 else ""))
                continue

            interval_trees[chrom].add(Interval(tx_start_0based, tx_end_1based + 0.1, gene_name))

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

#%%
#output_df[output_df['#NAME'] == "FGF16"]
output_path = re.sub(".gtf.gz$", "", os.path.basename(args.gtf_gz_path)) + ".txt.gz"

output_df.to_csv(output_path, index=False, sep="\t")

print(f"Wrote {len(output_df)} records to {os.path.abspath(output_path)}")

#%%
