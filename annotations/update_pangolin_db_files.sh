for i in *v42*.annotation.gtf.gz; do
    set -x
    python3 ~/code/Pangolin/scripts/create_db.py $i &
    set +x
done

wait
