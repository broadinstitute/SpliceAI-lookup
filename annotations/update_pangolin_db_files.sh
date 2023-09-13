gzcat gencode.v44lift37.annotation.gtf.gz | sed 's/chr//g' | bgzip > gencode.v44lift37.annotation.without_chr_prefix.gtf.gz

for p in gencode.v44.annotation.gtf.gz  gencode.v44lift37.annotation.without_chr_prefix.gtf.gz; do
    set -x
    python3 ~/code/Pangolin/scripts/create_db.py $p
    set +x
done

wait
