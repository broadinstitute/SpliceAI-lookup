set -x

NUM_THREADS=1    # set this to the number of cores on your machine (or a bit less)
HOST=127.0.0.1   # set this to 0.0.0.0 instead of 127.0.0.1 to allow access from other computers
PORT=8080        # set this to a port number that is not already in use
TIMEOUT=1800     # kill the server thread if it takes more than this many seconds to compute a response

# start the gunicorn server
gunicorn -w ${NUM_THREADS} -t ${TIMEOUT} -b ${HOST}:${PORT}  server:app
