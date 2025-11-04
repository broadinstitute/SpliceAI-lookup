#%%

import argparse
import collections
import gzip
import json
import os
import pandas as pd
import re

from annotation_utils.gtf_utils import parse_gtf


def main():
    p = argparse.ArgumentParser(description="""This script takes a Gencode .gtf.gz file
        and outputs an annotation file which can be passed to SpliceAI instead of 
        the default SpliceAI annotations which are still on Gencode v24.""")

    p.add_argument("-a", "--annotation-json-path", required=True, help="Path of the transcript annotations JSON file "
                   "created by generate_transcript_annotation_json.py")
    p.add_argument("--gtf-id-field", default="transcript_id", choices=["transcript_id", "gene_id"])
    p.add_argument("gtf_gz_path", help="Path of gene annotations file in GTF format")
    args = p.parse_args()

    for path in args.annotation_json_path, args.gtf_gz_path:
        if not os.path.exists(path):
            p.error(f"File not found: {path}")

    fopen = gzip.open if args.gtf_gz_path.endswith("gz") else open
    with fopen(args.annotation_json_path, "rt") as f:
        transcript_annotations = json.load(f)

    print(f"Parsing {args.gtf_gz_path}")
    gtf_id_to_exons = collections.defaultdict(set)
    for record in parse_gtf(os.path.expanduser(args.gtf_gz_path), "exon"):
        key = (record[args.gtf_id_field], record["strand"], record["chrom"])
        exon_tuple = (record["start"], record["end"])
        if exon_tuple in gtf_id_to_exons[key]:
            raise ValueError(f"Duplicate exon: {exon_tuple} in transcript {key}")
        gtf_id_to_exons[key].add(exon_tuple)

    output_records = []
    # SpliceAI predictions (prior to 'masking') depend only on transcript chrom/start/end/strand. Often, transcripts
    # within the same gene have the same chrom/start/end/strand and differ only in their internal exon structure.
    # We can discard these redundant transcripts (while making sure to keep all MANE Select and canonical transcripts).
    output_records_transcript_keys = set()
    maybe_output_records = []
    for (gtf_id, strand, chrom), exon_set in gtf_id_to_exons.items():
        tx_start_0based = min([start_1based - 1 for start_1based, _ in exon_set])
        tx_end_1based = max([end_1based for _, end_1based in exon_set])
        gtf_id_without_version = gtf_id.split(".")[0]
        if gtf_id_without_version not in transcript_annotations:
            print(f"WARNING: transcript {gtf_id_without_version} not found in {args.annotation_json_path}")
            continue
        transcript_annotation = transcript_annotations[gtf_id_without_version]
        if transcript_annotation['t_priority'] == "N":
            output_list = maybe_output_records
        else:
            output_list = output_records
            transcript_key = (chrom, strand, str(tx_start_0based), str(tx_end_1based))
            output_records_transcript_keys.add(transcript_key)

        # if it's a MANE Select, MANE Plus Clinical or canonical transcript
        exon_list = sorted(list(exon_set))
        exon_starts_0based = [start_1based - 1 for start_1based, _ in exon_list]
        exon_ends_1based = [end_1based for _, end_1based in exon_list]

        # reformat the records into a list which can be turned into a pandas DataFrame
        output_list.append({
            "#NAME": gtf_id,
            "CHROM": chrom,
            "STRAND": strand,
            "TX_START": str(tx_start_0based),
            "TX_END": str(tx_end_1based),
            "EXON_START": ",".join([str(s) for s in exon_starts_0based]) + ",",
            "EXON_END": ",".join([str(s) for s in exon_ends_1based]) + ",",
        })

    transcripts_kept_counter1 = len(output_records)
    transcripts_kept_counter2 = 0
    for output_record in maybe_output_records:
        transcript_key = (output_record["CHROM"], output_record["STRAND"], output_record["TX_START"], output_record["TX_END"])
        if transcript_key not in output_records_transcript_keys:
            # if this transcript has a chrom/start/end/strand that hasn't been seen before, add it to the output
            output_records_transcript_keys.add(transcript_key)
            output_records.append(output_record)
            transcripts_kept_counter2 += 1

    assert transcripts_kept_counter1 + transcripts_kept_counter2 == len(output_records)
    print(f"Kept {transcripts_kept_counter1:,d} transcripts which were MANE Select, MANE Plus Clinical or canonical.")
    print(f"Kept {transcripts_kept_counter2:,d} additional transcripts with unique transcript start/stop coords.")
    print(f"Discarded {len(maybe_output_records) - transcripts_kept_counter2:,d} out of {len(gtf_id_to_exons):,d} "
          f"({(len(maybe_output_records) - transcripts_kept_counter2) / len(gtf_id_to_exons):.1%}) transcripts "
          f"because they were redundant.")

    output_df = pd.DataFrame(output_records)
    output_df = output_df[["#NAME", "CHROM", "STRAND", "TX_START", "TX_END", "EXON_START", "EXON_END"]]
    output_path = re.sub(".gtf.gz$", "", os.path.basename(args.gtf_gz_path)) + ".txt.gz"
    output_df.to_csv(output_path, index=False, sep="\t")

    print(f"Wrote {len(output_df):,d} records to {os.path.abspath(output_path)}")


if __name__ == "__main__":
    main()
