def adder(fac):
    def add(dfac):
        nonlocal fac
        return fac + dfac

    return add


add5 = adder(5)
add10 = adder(10)

print(add5(5))
print(add5(6))
print(add5(7))

print(add10(2))
print(add10(3))
print(add10(4))


print("-------------")


def add(a):
    if a < 100:

        def add0(b):
            return add(a + b)

        return add0
    else:
        return a


print(add(101))
print(add(1)(1)(21)(99))
t = add(0)
while not isinstance(t, int):
    t = t(1)

print(t)


def add2bound(bound):
    def adder(a):
        if a > bound:
            return a
        else:

            def add0(b):
                return adder(a + b)

            return add0

    return adder


import random

random.randint(0, 50)


def just():
    t = add2bound(200)
    while not isinstance(t, int):
        t = t(random.randint(0, 3))
    print(t)


for i in range(10):
    just()


def hello():
    print("hello world")


hello()


def recur(n):
    if n == 1:
        return 1
    else:
        return n + recur(n - 1)


print(recur(400))

print("---------装饰器-----------")


def wrap(core_method):
    times = 1

    def wraped(*args):
        nonlocal times
        core_method(*args)
        print(f"函数{core_method}被调用{times}次")
        times += 1

    return wraped


@wrap
def core():
    print("执行核心业务")


@wrap
def core2(a):
    print(a)


core()
core()
core()
core()
core()

core2(12)
core2(13)
core2(14)
core2(15)

print("------------------ 带参数装饰器 --------")


routers = {}  # 全局变量 注册的路由-方法


def router(path):
    def wrap(method):
        routers[path] = method  # 注册路由

        def core(*args, **kwargs):
            print(f"访问路径{path}，携带参数{args}和{kwargs}")  # 方法增强
            method(*args, **kwargs)  # 核心业务

        return core

    return wrap


@router("/node1")
def n1(name):
    print(f"核心业务，处理节点1，来人是{name}")


@router("/node2")
def n2():
    print(f"核心业务，处理节点2")


# 看看路由注册
print(routers)

# 网关
def gateway(path, *args, **kwargs):
    m = routers.get(path)
    if m is None:
        print(f"错误，没有可以处理{path}的方法")
    else:
        m(*args, **kwargs)


# 外部网络访问
gateway("/node1", "zrx")
gateway("/node1", "mdk")
gateway("/node2")
gateway("/node3")


print("-----------lambda----------")

f = lambda a, b: a + b

print(f(1, 2))
print(f("a", "b"))
print(f([1, 2], [2, 1]))

f0 = lambda: print("空参λ")

f0()


print(list(map(lambda x: x * 2, [1,2])))


from functools import reduce

print(reduce(
    lambda x,y:x+y,
    range(11),0
))

print(list(
    filter(lambda x:x<10,[4,10,11])
))



import numpy as np

np.ndarray.sum = lambda seq : sum(seq)
arr = np.array([1,2,3])

print(arr)

print(arr.sum())