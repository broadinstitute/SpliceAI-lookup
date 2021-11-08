set -x

NUM_THREADS=1
HOST=127.0.0.1
PORT=8080
TIMEOUT=1800
gunicorn -w ${NUM_THREADS} -t ${TIMEOUT} -b ${HOST}:${PORT}  server:app
