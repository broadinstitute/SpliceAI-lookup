set -ex

spliceai -R ~/p1/ref/GRCh38/hg38.fa  -I test.vcf -O results.vcf -A ../annotations/gencode.v43.annotation.txt.gz

