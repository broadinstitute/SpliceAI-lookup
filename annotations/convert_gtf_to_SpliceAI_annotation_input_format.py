#%%

import argparse
import collections
import gzip
import os
import pandas as pd
import re


def parse_exons_from_gencode_gtf(gencode_gtf_path):
    transcript_type_summary_counter = collections.Counter()
    transcript_type_counter = collections.Counter()
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
            }

    print("Exon counts per transcript types:")
    for k, v in sorted(transcript_type_counter.items(), key=lambda i: -i[1]):
        print(f"{v:10d}: {k}")

    print("Summary:")
    for k, v in sorted(transcript_type_summary_counter.items(), key=lambda i: -i[1]):
        print(f"{v:10d}: {k}")


#%%

def main():
    p = argparse.ArgumentParser(description="""This script takes a Gencode .gtf.gz file
        and outputs an annotation file which can be passed to SpliceAI instead of 
        the default SpliceAI annotations which are still on Gencode v24.""")

    p.add_argument("gtf_gz_path", help="Path of gene annotations file in GTF format")
    args = p.parse_args()

    print(f"Parsing {args.gtf_gz_path}")

    transcript_id_to_exons = collections.defaultdict(set)
    for record in parse_exons_from_gencode_gtf(os.path.expanduser(args.gtf_gz_path)):
        transcript_id_without_version = record["transcript_id"].split(".")[0]
        key = (transcript_id_without_version, record["strand"], record["chrom"])
        exon_2_tuple = (record["start_1based"], record["end_1based"])
        if exon_2_tuple in transcript_id_to_exons[key]:
            raise ValueError(f"Duplicate exon: {exon_2_tuple} in transcript {key}")
        transcript_id_to_exons[key].add(exon_2_tuple)

    # reformat the aggregated records into a list which can be turned into a pandas DataFrame
    output_records = []
    for (transcript_id_without_version, strand, chrom), exon_set in transcript_id_to_exons.items():
        tx_start_0based = min([start_1based - 1 for start_1based, _ in exon_set])
        tx_end_1based = max([end_1based for _, end_1based in exon_set])

        exon_list = sorted(list(exon_set))
        exon_starts_0based = [start_1based - 1 for start_1based, _ in exon_list]
        exon_ends_1based = [end_1based for _, end_1based in exon_list]

        output_records.append({
            "#NAME": transcript_id_without_version,
            "CHROM": chrom,
            "STRAND": strand,
            "TX_START": str(tx_start_0based),
            "TX_END": str(tx_end_1based),
            "EXON_START": ",".join([str(s) for s in exon_starts_0based]) + ",",
            "EXON_END": ",".join([str(s) for s in exon_ends_1based]) + ",",
        })

    output_df = pd.DataFrame(output_records)
    output_df = output_df[["#NAME", "CHROM", "STRAND", "TX_START", "TX_END", "EXON_START", "EXON_END"]]
    output_path = re.sub(".gtf.gz$", "", os.path.basename(args.gtf_gz_path)) + ".txt.gz"
    output_df.to_csv(output_path, index=False, sep="\t")

    print(f"Wrote {len(output_df):,d} records to {os.path.abspath(output_path)}")


if __name__ == "__main__":
    main()
