import time

def f(n):
    pre = cur = 1
    for i in range(n-1):
        temp = cur
        cur = cur + pre
        pre = temp
    return cur

l = time.time()
a = f(40) # 165580141
l2 = time.time()

print(str(l2-l)+ "  " + str(a))