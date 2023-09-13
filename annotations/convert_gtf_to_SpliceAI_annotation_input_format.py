#%%

import argparse
import collections
import os
import pandas as pd
import re

from bw2_annotation_utils.gtf_utils import parse_gtf


def main():
    p = argparse.ArgumentParser(description="""This script takes a Gencode .gtf.gz file
        and outputs an annotation file which can be passed to SpliceAI instead of 
        the default SpliceAI annotations which are still on Gencode v24.""")

    p.add_argument("gtf_gz_path", help="Path of gene annotations file in GTF format")
    args = p.parse_args()

    print(f"Parsing {args.gtf_gz_path}")

    transcript_id_to_exons = collections.defaultdict(set)
    for record in parse_gtf(os.path.expanduser(args.gtf_gz_path), "exon"):
        key = (record["transcript_id"], record["strand"], record["chrom"])
        exon_2_tuple = (record["start"], record["end"])
        if exon_2_tuple in transcript_id_to_exons[key]:
            raise ValueError(f"Duplicate exon: {exon_2_tuple} in transcript {key}")
        transcript_id_to_exons[key].add(exon_2_tuple)

    # reformat the aggregated records into a list which can be turned into a pandas DataFrame
    output_records = []
    for (transcript_id, strand, chrom), exon_set in transcript_id_to_exons.items():
        tx_start_0based = min([start_1based - 1 for start_1based, _ in exon_set])
        tx_end_1based = max([end_1based for _, end_1based in exon_set])

        exon_list = sorted(list(exon_set))
        exon_starts_0based = [start_1based - 1 for start_1based, _ in exon_list]
        exon_ends_1based = [end_1based for _, end_1based in exon_list]

        output_records.append({
            "#NAME": transcript_id,
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
