import redis
import time

r = redis.Redis(host='localhost')

while True:
    for key in r.keys("request *liftover*"):
        print("Deleting key: ", key.decode("UTF-8"))
        r.delete(key)
    time.sleep(1)
