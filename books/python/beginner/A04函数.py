import random


def breakLine():
    print("-----------------")


def generate_random(number=5, max=10):
    for i in range(number):
        print(random.randint(0, max))


# <function generate_random at 0x0000024D95877F70>
print(generate_random)

generate_random()
generate_random(10)
generate_random(5, 100)

breakLine()

l = [1, 2]
l2 = l
l2.append(3)
print(l)


def add(*args):
    print(args)


add()
add(1)
add(1, 2)
add(1, 2, 3)


def add(*args):
    return sum(args)


print(add())
print(add(1))
print(add(1, 2))
print(add(1, 2, 3))


def fun(a, b):
    print(a, b)


fun(1, 2)
fun(b=1, a=2)


def fun(**kwargs):
    print(kwargs)


fun(a=1)
fun(a="1", b="aa0")

d = {"name": "mdk", "age": 12}
fun(**d)
