set -x

redis-cli flushall  #  clear all keys from redis

gunicorn -w 5 -t 1800 -b 0.0.0.0:80  -b 0.0.0.0:443 \
  --keyfile=../spliceailookup-api.broadinstitute.org.key \
  --certfile=../spliceailookup-api_broadinstitute_org.cert \
  server:app
