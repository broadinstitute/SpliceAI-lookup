set -ex

python3 -c 'pd.read_table("PrimateAI_3D.hg19.txt.gz").sort_values(["chr", "pos"], ascending=[True, True]).to_csv("PrimateAI_3D.hg19.sorted.txt.gz", sep="\t", header=True, index=False)'
gunzip -c PrimateAI_3D.hg19.sorted.txt.gz | bgzip > PrimateAI_3D.hg19.txt.gz
tabix -S 1 -s 1 -b 2 -e 2 PrimateAI_3D.hg19.txt.gz
tabix -S 1 -s 1 -b 2 -e 2 PrimateAI_3D.hg38.txt.gz

