set -ex

gencode_version=v44
gzcat gencode.${gencode_version}lift37.basic.annotation.gtf.gz | sed 's/chr//g' | bgzip > gencode.${gencode_version}lift37.basic.annotation.without_chr_prefix.gtf.gz
for p in gencode.${gencode_version}.basic.annotation.gtf.gz  gencode.${gencode_version}lift37.basic.annotation.without_chr_prefix.gtf.gz; do
    set -x
    python3 ~/code/Pangolin/scripts/create_db.py $p &
    set +x
done

wait
