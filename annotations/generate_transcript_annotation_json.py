#%%

"""This script creates a json file with annotation fields for each transcript id, including
whether it's a "MANE select" transcript, canonical transcript, etc.
"""

import argparse
import gzip
import json
import os
import re

from annotation_utils.get_ensembl_db_info import get_gene_id_to_canonical_transcript_id, \
    get_ensembl_ENST_to_RefSeq_ids
from annotation_utils.get_MANE_table import get_MANE_ensembl_transcript_table
from annotation_utils.gtf_utils import parse_gtf


# to get the latest database name, run:
#    mysql -h useastdb.ensembl.org -u anonymous -e "show databases;" | grep homo_sapiens_core
DEFAULT_ENSEMBL_DATABASE = "homo_sapiens_core_115_38"

# this is used to get the list of MANE select and MANE plus clinical ENST transcript ids.
DEFAULT_MANE_URL_BASE = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4"
DEFAULT_MANE_SUMMARY_TABLE_FILENAME = "MANE.GRCh38.v1.4.summary.txt.gz"

def main():
    p = argparse.ArgumentParser(description="""This script takes a Gencode .gtf.gz file
        and outputs an annotation file which can be passed to SpliceAI instead of 
        the default SpliceAI annotations which are still on Gencode v24. 
    """)
    p.add_argument("--mane-url-base", default=DEFAULT_MANE_URL_BASE)
    p.add_argument("--mane-summary-table-filename", default=DEFAULT_MANE_SUMMARY_TABLE_FILENAME)
    p.add_argument("-e", "--ensembl-database", default=DEFAULT_ENSEMBL_DATABASE)
    p.add_argument("gtf_gz_path", help="Path of gene annotations file in GTF format")
    args = p.parse_args()

    mane_summary_table_url = os.path.join(args.mane_url_base, args.mane_summary_table_filename)
    MANE_df = get_MANE_ensembl_transcript_table(mane_summary_table_url=mane_summary_table_url)

    print(f"Initalizing transcript priority annotation function")
    compute_transcript_priority = get_transcript_priority_annotation_function(
        ensembl_database=args.ensembl_database, MANE_df=MANE_df)

    esnembl_ENST_to_RefSeq_ids = get_ensembl_ENST_to_RefSeq_ids(database=args.ensembl_database)
    print(f"Downloaded {len(esnembl_ENST_to_RefSeq_ids):,d} ENST to RefSeq mappings")
    for key, refseq_ids in esnembl_ENST_to_RefSeq_ids.items():
        esnembl_ENST_to_RefSeq_ids[key] = list(sorted(refseq_ids))

    MANE_df["ensembl_ENST_without_version"] = MANE_df["Ensembl_nuc"].apply(lambda s: s.split(".")[0])
    MANE_ensembl_ENST_to_RefSeq_id = dict(MANE_df[["ensembl_ENST_without_version", "RefSeq_nuc"]].itertuples(index=False))
    print(f"Got {len(MANE_ensembl_ENST_to_RefSeq_id):,d} MANE ENST to RefSeq mappings, of which "
          f"{len(set(MANE_ensembl_ENST_to_RefSeq_id) - set(esnembl_ENST_to_RefSeq_ids)):,d} are unique.")
    MANE_ensembl_ENST_to_RefSeq_id = {k: [v] for k, v in MANE_ensembl_ENST_to_RefSeq_id.items()}
    esnembl_ENST_to_RefSeq_ids.update(MANE_ensembl_ENST_to_RefSeq_id)

    gene_coordinates_lookup = {}  # used for checking coordinate consistency between gene and transcript records (eg. GPR143)
    max_transcript_coordinates_lookup = {}
    for record in parse_gtf(os.path.expanduser(args.gtf_gz_path), feature_type="gene"):
        gene_id_without_version = record["gene_id"].split(".")[0]
        gene_coordinates_lookup[(record["gene_name"], gene_id_without_version)] = (record["chrom"], record["start"], record["end"])

    print(f"Parsing {args.gtf_gz_path}")
    output_json = {}
    for record in parse_gtf(os.path.expanduser(args.gtf_gz_path), feature_type="transcript"):
        transcript_id_without_version = record["transcript_id"].split(".")[0]
        transcript_priority = compute_transcript_priority(transcript_id=transcript_id_without_version)
        refseq_transcript_ids = esnembl_ENST_to_RefSeq_ids.get(transcript_id_without_version)
        output_json[transcript_id_without_version] = {
            "g_name": record["gene_name"],
            "g_id": record["gene_id"],
            "t_id": record["transcript_id"],
            "t_type": record["transcript_type"],
            "t_strand": record["strand"],
            "t_priority": transcript_priority,
            "t_refseq_ids": refseq_transcript_ids,
        }

        gene_id_without_version = record["gene_id"].split(".")[0]
        if gene_id_without_version in max_transcript_coordinates_lookup:
            transcript_chrom, transcript_start, transcript_end = max_transcript_coordinates_lookup[gene_id_without_version]
            max_transcript_coordinates_lookup[(record["gene_name"], gene_id_without_version)] = (
                transcript_chrom, min(transcript_start, record["start"]), max(transcript_end, record["end"]))
        else:
            max_transcript_coordinates_lookup[(record["gene_name"], gene_id_without_version)] = (
                record["chrom"], record["start"], record["end"])

    # check consistency of gene vs. transcript start/end coordinates
    warning1_counter = 0
    warning2_counter = 0
    for gene_name, gene_id_without_version in gene_coordinates_lookup:
        key = gene_name, gene_id_without_version
        gene_chrom, gene_start, gene_end = gene_coordinates_lookup[key]
        transcript_chrom, transcript_start, transcript_end = max_transcript_coordinates_lookup.get(key, (None, None, None))
        if transcript_chrom is None:
            warning1_counter += 1
            print(f"WARNING: Gene {gene_name} ({gene_id_without_version}) has no transcript records")
        elif transcript_start > gene_start or transcript_end < gene_end:
            warning2_counter += 1
            #start_diff = f". Start diff is {transcript_start - gene_start:,d}bp" if transcript_start - gene_start else ""
            #end_diff = f". End diff is {gene_end - transcript_end:,d}bp" if gene_end - transcript_end else ""
            #print(f"WARNING: Gene {gene_name} ({gene_id_without_version}) has inconsistent coordinates: "
            #      f"gene={gene_chrom}:{gene_start}-{gene_end}, "
            #      f"transcript={transcript_chrom}:{transcript_start}-{transcript_end}"
            #      f"{start_diff}{end_diff}")
    if warning1_counter > 0:
        print(f"WARNING: {warning1_counter:,d} genes don't have any transcript records")
    if warning2_counter > 0:
        print(f"WARNING: {warning2_counter:,d} out of {len(gene_coordinates_lookup):,d} genes have a genomic interval "
              f"that is wider than the interval of any of their transcripts.")

    output_path = re.sub(".gtf.gz$", "", os.path.basename(args.gtf_gz_path)) + ".transcript_annotations.json.gz"
    with gzip.open(output_path, "wt") as f:
        json.dump(output_json, f, indent=4, sort_keys=True)

    print(f"Done writing {len(output_json):,d} transcript annotations to {output_path}")


def get_transcript_priority_annotation_function(ensembl_database, MANE_df):
    """Initializes annotation data and returns the compute_transcript_priority function."""

    MANE_select_transcript_ids = {
        t_id.split(".")[0] for t_id in MANE_df[MANE_df["MANE_status"] == "MANE Select"]["Ensembl_nuc"]}
    print(f"Got {len(MANE_select_transcript_ids):,d} MANE select transcript ids")

    MANE_plus_clinical_transcript_ids = {
        t_id.split(".")[0] for t_id in MANE_df[MANE_df["MANE_status"] == "MANE Plus Clinical"]["Ensembl_nuc"]}
    print(f"Got {len(MANE_plus_clinical_transcript_ids):,d} MANE plus clinical transcript ids")

    gene_id_to_canonical_transcript_id = get_gene_id_to_canonical_transcript_id(database=ensembl_database)
    canonical_transcript_ids = {
        t_id.split(".")[0] for t_id in gene_id_to_canonical_transcript_id.values()}
    print(f"Got {len(canonical_transcript_ids):,d} canonical transcript ids")

    def compute_transcript_priority(transcript_id):
        """Returns a string indicating the priority of the given transcript.
        The return value can be (in order from higher to lower priority):
            "MS" (for MANE select)
            "MP" (for MANE plus clinical)
            "C" (for canonical)
            "N" (for none of the above)
        """
        transcript_id = transcript_id.split(".")[0]

        if transcript_id in MANE_select_transcript_ids:
            transcript_priority = "MS"
        elif transcript_id in MANE_plus_clinical_transcript_ids:
            transcript_priority = "MP"
        elif transcript_id in canonical_transcript_ids:
            transcript_priority = "C"
        else:
            transcript_priority = "N"

        return transcript_priority

    return compute_transcript_priority


if __name__ == "__main__":
    main()
