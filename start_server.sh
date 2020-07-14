set -x

gunicorn -t 300 -w 2 -b 0.0.0.0:80 server:app
