set -ex
mysql -h useastdb.ensembl.org -u anonymous -e "show databases;" | grep -i homo_sapiens_core
