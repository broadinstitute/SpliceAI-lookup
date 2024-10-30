import argparse
import gzip
import os
import pandas as pd
import subprocess
import tqdm

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-g", "--gene-thresholds-csv", default="PrimateAI_3D.per_gene_percentile_thresholds.csv.gz",
						help="Input CSV file")
	parser.add_argument("-p", "--show-progress-bar", action="store_true", help="Show progress bar")

	parser.add_argument("scores_table", help="Table with precomputed PrimateAI_3D scores")
	args = parser.parse_args()

	if not os.path.isfile(args.gene_thresholds_csv):
		parser.error(f"Gene thresholds CSV file not found: {args.gene_thresholds_csv}")

	if not args.scores_table.endswith(".txt.gz"):
		parser.error("Scores table must have a .txt.gz extension")

	if not os.path.isfile(args.scores_table):
		parser.error(f"Scores file not found: {args.scores_table}")

	return args

def main():

	args = parse_args()

	# Load gene thresholds
	print(f"Parsing {args.gene_thresholds_csv}")
	gene_thresholds_df = pd.read_csv(args.gene_thresholds_csv)
	transcript_to_percentile_threshold_map = dict(
		zip(gene_thresholds_df['Transcript'], gene_thresholds_df['PAI3D_Gene_Percentile_Threshold']))

	print(f"Loaded {len(gene_thresholds_df):,d} gene thresholds")
	print(f"Parsing {args.scores_table}")
	with gzip.open(args.scores_table, "rt") as f, open(f"{args.scores_table}.unfinished", "wt") as out_f:
		header = next(f).strip().split("\t")
		transcript_id_index = 4
		percentile_index = 9
		assert header[transcript_id_index] == "gene_name"
		assert header[percentile_index] == "percentile_PAI3D"

		if args.show_progress_bar:
			f = tqdm.tqdm(f, unit=" lines", unit_scale=True)

		output_columns = [
			"chrom", "pos", "ref", "alt", "PAI3D_percentile", "PAI3D_gene_threshold",
		]
		out_f.write("\t".join(output_columns) + "\n")

		for i, line in enumerate(f):
			fields = line.strip().split("\t")
			output_row = fields[:4]
			transcript_id = fields[transcript_id_index]

			if transcript_id not in transcript_to_percentile_threshold_map:
				raise ValueError(f"Transcript ID {transcript_id} from {args.score_table} not found in {args.gene_thresholds_csv}")

			percentile = fields[percentile_index]
			output_row += [percentile, f"{float(transcript_to_percentile_threshold_map[transcript_id]):0.3f}"]
			out_f.write("\t".join(output_row) + "\n")

	output_table_path = args.scores_table.replace(".txt.gz", "") + ".with_gene_thresholds.txt.gz"
	subprocess.check_output(f"bgzip {args.scores_table}.unfinished", shell=True)
	subprocess.check_output(f"mv {args.scores_table}.unfinished.gz {output_table_path}", shell=True)
	subprocess.check_output(f"tabix -f -S 1 -s 1 -b 2 -e 2 {output_table_path}", shell=True)
	#subprocess.check_output(f"gsutil -m cp {output_table_path} {output_table_path}.tbi gs://spliceai-lookup-reference-data/", shell=True)

if __name__ == "__main__":
	main()
