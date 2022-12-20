set -x

#redis-cli flushall  #  clear all keys from redis

gunicorn -w 8 -t 1800 -b 0.0.0.0:80  -b 0.0.0.0:443 \
  --keyfile=../spliceailookup-api.broadinstitute.org.key \
  --certfile=../spliceailookup-api.broadinstitute.org.crt \
  server:app
