import sys

print(sys.path)
print(sys.version)

import time

print(time.time())
print(time.time_ns())
time.sleep(0.01)
print(time.time_ns())

# Wed Nov 18 20:42:49 2020
print(time.ctime(time.time())) 

#time.struct_time(tm_year=2020, tm_mon=11, tm_mday=18, 
# tm_hour=20, tm_min=43, tm_sec=32, tm_wday=2, tm_yday=323, tm_isdst=0)
print(time.localtime(time.time()))
print(time.localtime(time.time()).tm_mday)
print(time.localtime(time.time()).tm_hour)

import random

print(random.random()) # 0.8247688584701691

print(random.randrange(1,10)) # 4

arr = list(range(100,200))
print(random.choice(arr)) # 116

import hashlib
# 加密算法

print(ord("我")) # 25105

print(hashlib.md5("中文".encode('utf-8')).hexdigest()) # a7bac2239fcdcb3a067903d8077c4a07