
gcloud --project spliceai-lookup-412920 sql instances list

# view settings
gcloud --project spliceai-lookup-412920 sql instances describe spliceai-lookup-db

# adjust DB settings
# https://cloud.google.com/sql/docs/postgres/flags#gcloud
gcloud --project spliceai-lookup-412920 sql instances patch spliceai-lookup-db --database-flags=max_connections=50

# restart DB
gcloud --project spliceai-lookup-412920 sql instances restart spliceai-lookup-db
