import redis
import time

r = redis.Redis(host='localhost')

while True:
    for key in r.keys("request *"):
        key = key.decode("UTF-8")
        if "liftover" in key:
            print("Deleting key: ", key)
            r.delete(key)
    time.sleep(1)
