export PYTHONPATH=~/code/methods/:$PYTHONPATH
set -x
for p in gencode.v43.annotation.gtf.gz  gencode.v43lift37.annotation.gtf.gz; do
    time python3 convert_gencode_gtf_to_spliceai_annotation_input_file.py $p  | tee process_$(echo $p | sed s/.gtf.gz//).log
done
