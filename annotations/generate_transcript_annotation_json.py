#%%

"""This script creates a json file with annotation fields for each transcript id, including
whether it's a "MANE select" transcript, canonical transcript, etc.
"""

import argparse
import collections
from collections import defaultdict, Counter
import gzip
from intervaltree import IntervalTree, Interval
import os
import pandas as pd
import re
import sys
sys.path.append(os.path.expanduser("~/code/SpliceAI"))          # code from https://github.com/bw2/SpliceAI
sys.path.append(os.path.expanduser("~/code/annotation-utils"))  # code from https://github.com/bw2/annotation-utils

# import from https://github.com/bw2/annotation-utils
from annotations.get_ensembl_db_info import get_gene_id_to_canonical_transcript_id, CURRENT_ENSEMBL_DATABASE
from annotations.get_MANE_table import get_MANE_ensembl_transcript_table

# cd to the dir that contains this script
os.chdir(os.path.expanduser("~/code/SpliceAI-lookup/annotations"))

#%%

test_args = None
#test_args = ["./annotations/gencode.v24.annotation.gtf.gz"]
#test_args = ["./annotations/gencode.v36lift37.annotation.gtf.gz"]
#test_args = ["./annotations/gencode.v36.annotation.gtf.gz"]

#%%

test_args = ["./gencode.v44.basic.annotation.gtf.gz"]
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

# update the ENSG to RefSeq id map
os.system("./download_ENSG_to_RefSeq_mapping.sh")

ENST_to_refseq_map = collections.defaultdict(set)
with open("./ENST_to_RefSeq_map.txt", "rt") as f:
    for line in f:
        ENST_id, refseq_id = line.strip().split("\t")
        ENST_to_refseq_map[ENST_id].add(refseq_id)

#%%

# Retrieve all MANE_select transcripts
MANE_df = get_MANE_ensembl_transcript_table("transcript")
MANE_df = MANE_df[MANE_df.tag.isin({"MANE_Plus_Clinical", "MANE_Select"})]

MANE_select_gene_ids_to_transcript_ids_without_version_number = {
    row.gene_id.split(".")[0]: row.transcript_id.split(".")[0]
    for _, row in MANE_df[MANE_df.tag == "Mane_Select"][["gene_id", "transcript_id"]].iterrows()
}

MANE_plus_clinical_gene_ids_to_transcript_ids_without_version_number = {
    row.gene_id.split(".")[0]: row.transcript_id.split(".")[0]
    for _, row in MANE_df[MANE_df.tag == "MANE_Plus_Clinical"][["gene_id", "transcript_id"]].iterrows()
}

#%%

def parse_gencode_file(gencode_gtf_path):
    transcript_type_summary_counter = Counter()
    transcript_type_counter = Counter()
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

            transcript_type_counter[meta_fields["transcript_type"]] += 1
            if meta_fields["transcript_type"] == "protein_coding":
                priority = "primary"
                transcript_type_summary_counter["Including new protein-coding gene"] += 1
            else:
                priority = "secondary"
                transcript_type_summary_counter["Including new other gene"] += 1

            yield {
                "chrom": chrom,
                "start_1based": start_1based,
                "end_1based": end_1based,
                "annotation_source": annotation_source,
                "strand": strand,
                "gene_id": meta_fields["gene_id"],
                "transcript_id": meta_fields["transcript_id"],
                "gene_name": meta_fields["gene_name"],
                "gene_type": meta_fields["transcript_type"],
                "transcript_type":  meta_fields["transcript_type"],
                "priority": priority,
            }

    print("Exon counts per transcript types:")
    for k, v in sorted(transcript_type_counter.items(), key=lambda i: -i[1]):
        print(f"{v:10d}: {k}")

    print("Summary:")
    for k, v in sorted(transcript_type_summary_counter.items(), key=lambda i: -i[1]):
        print(f"{v:10d}: {k}")

    #pprint(list(gene_type_counter.keys()))


#%%

# aggregate gtf exon records into buckets keyed by (chrom, gene name, strand)
all_exons_by_priority = {
    "primary": defaultdict(lambda: defaultdict(set)),
    "secondary": defaultdict(lambda: defaultdict(set)),
}

print(f"Getting canonical transcripts from {CURRENT_ENSEMBL_DATABASE}")
gene_id_to_canonical_transcript_id = get_gene_id_to_canonical_transcript_id()

genes_without_canonical_transcripts = 0
genes_without_MANE_select_transcripts = 0
total_transcripts = 0

gene_set_all = set()
gene_set_without_MANE_select = set()

transcript_id_to_composite_name = {}
for record in parse_gencode_file(os.path.expanduser(args.gtf_gz_path)):
    priority = record["priority"]
    transcript_type = record["transcript_type"]

    gene_id_without_version = record["gene_id"].split(".")[0]

    total_transcripts += 1
    gene_set_all.add(gene_id_without_version)

    transcript_id_without_version = record["transcript_id"].split(".")[0]

    # set the is_main_transcript variable
    is_main_transcript = "N"  # possible values are "MS" (for MANE select), "MP" (for MANE plus), "C" (for canonical) and "N" (for none of the above)
    if is_main_transcript == "N":
        if gene_id_without_version not in MANE_select_gene_ids_to_transcript_ids_without_version_number:
            #print("WARNING: no MANE select transcript for " + record["gene_id"])
            genes_without_MANE_select_transcripts += 1
            gene_set_without_MANE_select.add(gene_id_without_version)
        elif transcript_id_without_version == MANE_select_gene_ids_to_transcript_ids_without_version_number[gene_id_without_version]:
            is_main_transcript = "MS"

    if is_main_transcript == "N":
        if gene_id_without_version not in MANE_plus_clinical_gene_ids_to_transcript_ids_without_version_number:
            pass
        elif transcript_id_without_version == MANE_plus_clinical_gene_ids_to_transcript_ids_without_version_number[gene_id_without_version]:
            is_main_transcript = "MP"

    if is_main_transcript == "N":
        if gene_id_without_version not in gene_id_to_canonical_transcript_id:
            print(f"WARNING: no canonical transcript for " + record["gene_id"])
            genes_without_canonical_transcripts += 1
        elif transcript_id_without_version == gene_id_to_canonical_transcript_id[gene_id_without_version]:
            is_main_transcript = "C"


    refseq_transcript_ids_set = ENST_to_refseq_map[transcript_id_without_version]
    name = "---".join([record["gene_name"], record["gene_id"], record["transcript_id"], is_main_transcript, record["transcript_type"], ",".join(refseq_transcript_ids_set)])
    key = (record["chrom"], name, record["strand"])
    all_exons_by_priority[priority][transcript_type][key].add((int(record['start_1based']), int(record['end_1based'])))
    transcript_id_to_composite_name[transcript_id_without_version] = name


print(f"{genes_without_canonical_transcripts:,d} out of {total_transcripts:,d} transcripts are from genes without canonical transcripts")
print(f"{genes_without_MANE_select_transcripts:,d} out of {total_transcripts:,d} transcripts are from genes without MANE select transcripts")
print(f"{len(gene_set_without_MANE_select):,d} out of {len(gene_set_all):,d} genes don't have a MANE select transcript")

#%%

def transcript_type_order(transcript_type):
    try:
        return TRANSCRIPT_TYPES_BY_PRIORITY.index(transcript_type)
    except ValueError:
        return len(TRANSCRIPT_TYPES_BY_PRIORITY) + 1


# reformat the aggregated records into a list which can be turned into a pandas DataFrame
output_records = []
used_transcript_type_counter = Counter()
if SKIP_SECONDARY_TRANSCRIPTS:
    interval_trees = defaultdict(IntervalTree)
    skipped_transcript_type_counter = Counter()
for priority in all_exons_by_priority:
    for transcript_type in sorted(all_exons_by_priority[priority].keys(), key=transcript_type_order):
        current_exon_sets = all_exons_by_priority[priority][transcript_type]
        for key in sorted(current_exon_sets.keys()):
            exons_set = current_exon_sets[key]

            chrom, gene_name, strand = key

            tx_start_0based = min([start_1based - 1 for start_1based, _ in exons_set])
            tx_end_1based = max([end_1based for _, end_1based in exons_set])

            # check for overlap with previously-added transcripts
            if SKIP_SECONDARY_TRANSCRIPTS:
                if priority != "primary" and interval_trees[chrom].overlaps(Interval(tx_start_0based, tx_end_1based)):
                    # skip any secondary transcripts that overlap already-added primary transcripts
                    overlapping_genes = sorted(set([i.data for i in interval_trees[chrom][tx_start_0based:tx_end_1based]]))
                    #print(f"Skipping {priority} {transcript_type} gene {gene_name} since it overlaps {len(overlapping_genes)} gene(s): " +
                    #    ", ".join(overlapping_genes[:5]) + ("..." if len(overlapping_genes) > 5 else ""))
                    skipped_transcript_type_counter[transcript_type] += 1
                    continue
                interval_trees[chrom].add(Interval(tx_start_0based, tx_end_1based + 0.1, gene_name))

            used_transcript_type_counter[transcript_type] += 1

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


print("Transcript types counter:")
for k, v in sorted(used_transcript_type_counter.items(), key=lambda i: -i[1]):
    print(f"{v:10d}: {k}")

if SKIP_SECONDARY_TRANSCRIPTS:
    print("Skipped transcript types counter:")
    for k, v in sorted(skipped_transcript_type_counter.items(), key=lambda i: -i[1]):
        print(f"{v:10d}: {k}")

#%%

# generate output table for SpliceAI

output_df = pd.DataFrame(output_records)
output_df = output_df[["#NAME", "CHROM", "STRAND", "TX_START", "TX_END", "EXON_START", "EXON_END"]]

#output_df[output_df['#NAME'] == "FGF16"]
output_path = re.sub(".gtf.gz$", "", os.path.basename(args.gtf_gz_path)) + ".txt.gz"

output_df.to_csv(output_path, index=False, sep="\t")

print(f"Wrote {len(output_df):,d} records to {os.path.abspath(output_path)}")

#%%


# generate a GTF for pangolin where the gene name has been replaced with the composite name
pangolin_gtf_output_path = re.sub(".gtf.gz$", "", os.path.basename(args.gtf_gz_path)) + ".for_pangolin.gtf.gz"
with gzip.open(args.gtf_gz_path, "rt") as f, gzip.open(pangolin_gtf_output_path, "wt") as f2:
    for line in f:
        if line.startswith("#"):
            print(line.strip())
            continue
        fields = line.strip().split("\t")
        if fields[2] != "transcript":
            continue

        meta_fields = {}
        for meta_field in fields[8].strip("; ").split(";"):
            key, value = meta_field.strip().replace('"', '').split()
            meta_fields[key] = value

        transcript_id_without_version = meta_fields["transcript_id"].split(".")[0]
        if transcript_id_without_version not in transcript_id_to_composite_name:
            print(f"WARNING: no composite name for {transcript_id_without_version}")
            continue

        meta_fields["gene_name"] = transcript_id_to_composite_name[transcript_id_without_version]
        meta_fields["gene_id"] = transcript_id_to_composite_name[transcript_id_without_version]
        meta_fields["transcript_id"] = transcript_id_to_composite_name[transcript_id_without_version]

        fields[8] = "; ".join([f'{k} "{v}"' for k, v in meta_fields.items()])
        f2.write("\t".join(fields) + "\n")

print(f"Done writing to {pangolin_gtf_output_path}")
#%%