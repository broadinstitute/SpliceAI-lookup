redis-cli flushall

kill -HUP $(pgrep gunicorn | head -n 1)
