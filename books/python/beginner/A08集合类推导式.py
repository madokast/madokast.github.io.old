names = ["abc", "qqqqqq", "aa"]

shortNames = [n for n in names if len(n) <= 3]

print(shortNames)

lens = [len(n) for n in names]

print(lens)

print([i for i in range(101) if i % 3 == 0])

print("---------------")

print([(x, y) for x in range(5) if x % 2 == 0 for y in range(5) if y % 2 != 0])

arr1 = ["a", "b", "c"]
arr2 = [1, 2, 3]
cross = [(x, y) for x in arr1 for y in arr2]
print(cross)


m = {10: "aaa", 2: "b"}
print({k: len(v) for k, v in m.items()})

print("------------------生成器------------------")

gene = (x * 3 for x in range(20))

print(gene)

print(gene.__next__())

print(next(gene))

for i in gene:
    print(i, end=" ")

print()

gene = (x * 3 for x in range(20))

print(list(gene))


print("---------------------函数yield生成器-------------")


def addOne(length=5):
    n = 0
    while n < length:
        n += 1
        yield n


print(addOne())  # <generator object addOne at 0x000001DDBBCD94A0>

for i in addOne():
    print(i)

print("---------------------函数yield生成器 返回值-------------")


def sendYield():
    while True:
        sent = yield 5
        print(f"sent = {sent}")
        if sent == 4:
            break

gene = sendYield()

print(gene.send(None))
print(gene.send(1))
print(gene.send(2))
print(gene.send(3))
# print(gene.send(4)) # StopIteration
# print(gene.send(5)) # StopIteration

print('----------协程-----------------')

def task1(n):
    for i in range(n):
        print(f"搬砖{i}")
        yield

def task2(n):
    for i in range(n):
        print(f"听歌{i}")
        yield

g1 = task1(5)
g2 = task2(5)

for i in range(5):
    g1.__next__()
    g2.__next__()