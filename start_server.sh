set -x

gunicorn -w 6 -t 600 -b 0.0.0.0:80  -b 0.0.0.0:443 \
  --keyfile=../spliceailookup-api.broadinstitute.org.key \
  --certfile=../spliceailookup-api_broadinstitute_org.cert \
  server:app
