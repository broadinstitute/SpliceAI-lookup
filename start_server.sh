set -x

gunicorn -t 300 -w 1 -b 0.0.0.0:80 server:app
