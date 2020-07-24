set -x

gunicorn -t 600 -w 2 -b 0.0.0.0:80 server:app
