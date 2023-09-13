To generate or update the transcript annotation files needed for running SpliceAI and Pangolin:

1. Download the latest "basic" gene annotations in GTF format from Gencode for both [GRCh38](https://www.gencodegenes.org/human/) and GRCh37.
2. Update the Gencode version string at the top of these bash scripts, and then run them:
    - update_json_annotation_files.sh
    - update_SpliceAI_annotation_txt_files.sh
    - update_Pangolin_db_files.sh
3. Update the GENCODE_VERSION string in ../server.py