gzcat gencode.v43lift37.annotation.gtf.gz | sed 's/chr//g' | bgzip > gencode.v43lift37.annotation.without_chr_prefix.gtf.gz

for i in gencode.v43.annotation.gtf.gz  gencode.v43lift37.annotation.without_chr_prefix.gtf.gz; do
    set -x
    python3 ~/code/Pangolin/scripts/create_db.py $i 
    set +x
done

wait
