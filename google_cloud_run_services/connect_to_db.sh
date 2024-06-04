#set -ex

PGPASSWORD=$(cat .pgpass) psql -h 34.173.33.168 -d spliceai-lookup-db -U postgres -d spliceai-lookup-db


# useful queries:

# count variant consequences (counted once per variant)
# select variant_consequence, count(*) as c from (select variant_consequence, variant from log where length(variant_consequence) > 1 group by variant_consequence, variant) log group by variant_consequence order by c desc;

# count queries per ip per day
# select ip, logtime::timestamp::date, MAX(event_name), count(*) as c from log group by ip, logtime::timestamp::date ORDER BY logtime::timestamp::date desc, c desc


# check intergenic variants
# select distinct variant, genome from log where variant_consequence = 'intergenic' and genome='38';


# compute % of queried variants that are splice-region
# select c as splice_region_variants, d as total_variants, c::float/d::float as percent from ( select count(*) as c from ( select variant_consequence, variant from log where length(variant_consequence) > 1 group by variant_consequence, variant ) temp1 where variant_consequence = 'splice_region_variant' or variant_consequence = 'splice_donor_variant' or variant_consequence = 'splice_acceptor_variant' or variant_consequence = 'splice_polypyrimidine_tract_variant' or variant_consequence = 'splice_donor_region_variant' )  temp2 full outer join ( select count(*) as d from (	select variant from log where length(variant_consequence) > 1 group by variant_consequence, variant ) temp3 ) temp4  on 1=1;	
