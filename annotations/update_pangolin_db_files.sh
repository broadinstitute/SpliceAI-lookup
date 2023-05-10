for i in gencode.v43.annotation.gtf.gz  gencode.v43lift37.annotation.gtf.gz; do
    set -x
    python3 ~/code/Pangolin/scripts/create_db.py $i 
    set +x
done

wait
