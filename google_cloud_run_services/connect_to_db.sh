#set -ex

PGPASSWORD=$(cat .pgpass) psql -h 34.173.33.168 -d spliceai-lookup-db -U postgres -d spliceai-lookup-db


# useful queries:

# count variant consequences (counted once per variant)
# select variant_consequence, count(*) as c from (select variant_consequence, variant from log where length(variant_consequence) > 1 group by variant_consequence, variant) log group by variant_consequence order by c desc;

# count queries per ip per day
# select ip, logtime::timestamp::date, MAX(event_name), count(*) as c from log group by ip, logtime::timestamp::date ORDER BY logtime::timestamp::date desc, c desc
