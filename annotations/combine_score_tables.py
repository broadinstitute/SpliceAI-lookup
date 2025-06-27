import argparse
import gzip
import os
import tqdm

chrom_to_index = {f"chr{i}": i for i in range(1, 23)}
chrom_to_index["chrX"] = 23
chrom_to_index["chrY"] = 24
index_to_chrom = {i: v for v, i in chrom_to_index.items()}
base_to_index = {c: i for i, c in enumerate("ACGT")}
index_to_base = {i: c for c, i in base_to_index.items()}

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--genome-version", choices=["hg38", "hg19"], required=True)
args = parser.parse_args()

hg38_or_hg19 = args.genome_version

## Process hg38 tables
all_keys = set()
table1_lookup = {}
table1_path = os.path.expanduser(f"~/code/SpliceAI-lookup/annotations/PrimateAI_3D.{hg38_or_hg19}.with_gene_thresholds.txt.gz")

chrom_column = "chr" if hg38_or_hg19 == "hg38" else "chrom"
ref_column = "non_flipped_ref" if hg38_or_hg19 == "hg38" else "ref"
alt_column = "non_flipped_alt" if hg38_or_hg19 == "hg38" else "alt"
percentile_column = "percentile_PAI3D" if hg38_or_hg19 == "hg38" else "PAI3D_percentile"
gene_threshold_column = "PAI3D_Gene_Percentile_Threshold" if hg38_or_hg19 == "hg38" else "PAI3D_gene_threshold"

print(f"Reading table #1 from {table1_path}")
with gzip.open(table1_path, "rt") as f:
	header = next(f).strip().split("\t")
	header_indices = {c: i for i, c in enumerate(header)}
	counter = 0
	for line in tqdm.tqdm(f, unit=" lines", unit_scale=True, total=70_667_467):
		fields = line.rstrip().split("\t")
		chrom = fields[header_indices[chrom_column]]
		if "_" in chrom:
			# skip supercontigs
			continue

		key = (
			chrom_to_index[chrom],
			int(fields[header_indices["pos"]]),
			base_to_index[fields[header_indices[ref_column]]],
			base_to_index[fields[header_indices[alt_column]]]
		)
		all_keys.add(key)
		table1_lookup[key] = [
			float(fields[header_indices[percentile_column]]), 
			float(fields[header_indices[gene_threshold_column]]), 
		]
		counter += 1

print(f"Parsed {counter:,d} records from table #1")

table2_lookup = {}
if hg38_or_hg19 == "hg38":
	table2_path = os.path.expanduser("~/code/SpliceAI-lookup/annotations/promoterAI_tss500.tsv.gz")
else:
	table2_path = os.path.expanduser("~/code/SpliceAI-lookup/annotations/promoterAI_tss500_hg19.tsv.gz")

print(f"Reading table #2 from {table2_path}")
with gzip.open(table2_path, "rt") as f:
	header = next(f).strip().split("\t")
	header_indices = {c: i for i, c in enumerate(header)}
	counter = 0
	for line in tqdm.tqdm(f, unit=" lines", unit_scale=True, total=261_666_406):
		fields = line.rstrip().split("\t")
		chrom = fields[header_indices[chrom_column]]
		if "_" in chrom:
			# skip supercontigs
			continue

		key = (
			chrom_to_index[chrom],
			int(fields[header_indices["pos"]]),
			base_to_index[fields[header_indices["ref"]]],
			base_to_index[fields[header_indices["alt"]]]
		)
		all_keys.add(key)
		table2_lookup[key] = float(fields[header_indices["promoterAI"]])
		counter += 1

print(f"Parsed {counter:,d} records from table #2")
output_path = f"PrimateAI_and_PromoterAI_scores.{hg38_or_hg19}.tsv"
print(f"Writing output to {output_path}")
with open(output_path, "wt") as f:
	f.write("\t".join([
		"chrom",
		"pos",
		"ref",
		"alt",
		"PAI3D_percentile",
		"PAI3D_gene_threshold",
		"PromoterAI_score",
	]) + "\n")

	for key in tqdm.tqdm(sorted(all_keys), unit=" records", unit_scale=True, total=len(all_keys)):
		chrom_index, pos, ref_index, alt_index = key
		percentile_PAI3D, PAI3D_gene_threshold = table1_lookup.get(key, (None, None))
		promoterAI_score = table2_lookup.get(key, None)

		f.write("\t".join([
			index_to_chrom[chrom_index],
			str(pos),
			index_to_base[ref_index],
			index_to_base[alt_index],
			f"{percentile_PAI3D:.3f}" if percentile_PAI3D is not None else "",
			f"{PAI3D_gene_threshold:.2f}" if PAI3D_gene_threshold is not None else "",
			f"{promoterAI_score:.3f}" if promoterAI_score is not None else "",
		]) + "\n")

os.system(f"bgzip {output_path}")
os.system(f"tabix -f -S 1 -s 1 -b 2 -e 2 {output_path}.gz")

#%%
