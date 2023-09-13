set -ex

for p in gencode.v44.annotation.gtf.gz  gencode.v44lift37.annotation.gtf.gz; do
    time python3 convert_gtf_to_SpliceAI_annotation_input_format.py $p  | tee process_$(echo $p | sed s/.gtf.gz//).log
done
