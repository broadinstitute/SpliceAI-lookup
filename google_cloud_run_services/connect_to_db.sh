#set -ex

echo DB password: $(cat db_password.txt)
psql -h 34.173.33.168 -d spliceai-lookup-db -U postgres -W -d spliceai-lookup-db
