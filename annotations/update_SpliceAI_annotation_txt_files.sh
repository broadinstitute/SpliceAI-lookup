set -ex

# make sure annotation-utils is installed since it's a dependency of convert_gtf_to_SpliceAI_annotation_input_format.py
python3 -m pip install git+https://github.com/bw2/annotation-utils

gencode_version=v44
for p in gencode.${gencode_version}.basic.annotation.gtf.gz  gencode.${gencode_version}lift37.basic.annotation.gtf.gz; do
    log_path=process_$(echo ${p} | sed s/.gtf.gz//).log
    time python3 generate_transcript_annotation_json.py ${p} | tee -a ${log_path}
    json_path=$(echo ${p} | sed 's/.gtf.gz/.transcript_annotations.json/')
    time python3 convert_gtf_to_SpliceAI_annotation_input_format.py -a ${json_path} ${p}  | tee ${log_path}
done
