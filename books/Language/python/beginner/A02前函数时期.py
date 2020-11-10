var1 = "12a"

print(id(var1), var1)

var1 = var1 + "1"

print(id(var1), var1)

var2 = "aaa"

var3 = "aaa"

# py 也存在常量池 1721861033200 1721861033200 aaa aaa
print(id(var2), id(var3), var2, var3)

var4 = str(var3)
print(id(var3), id(var4), var3, var4)


def printIdAndSelf(obj):
    print(id(obj), obj)


# 1393037143024 123
# 1393037173424 123a

var5 = "123"

printIdAndSelf(var5)

var5 += "a"

printIdAndSelf(var5)

# 1988038453552 1
# 1988038453584 2

var6 = 1

printIdAndSelf(var6)

var6 += 1

printIdAndSelf(var6)

# -------------------- ==

var7 = "a" + "aa"

var8 = "a"

var8 = var8 + "aa"

printIdAndSelf(var7)
printIdAndSelf(var8)

print(var7 == var8) # T
print(id(var7) == var8) # F
print(var7 is var8) # F

var9 = 12
var10 = 12

printIdAndSelf(var9)
printIdAndSelf(var10)

print(var9 is var10) # T


# -------------------- bin

var11 = 128

print(var11,bin(var11))

var12 = 0.1

# print(var12,bin(var12))

# ----------------- bigInt

print(type(1))
print(type(1**10))
print(type(1**20))
print(type(1**100))


print(5)
print(~5)

print(5 if 3>1 else 1)
print(5 if 3<1 else 1)

# ------------------ for
for i in range(0,10):
    print(i,end=' ')

print()

for i in range(5):
    print(i,end=' ')

print()

for i in range(10):
    if(i==10):
        break
else:
    print("循环正常结束")



for i in range(10):
    if(i==5):
        break
else:
    print("循环正常结束")

print(list(range(10)))
print(list(range(5,10)))
print(list(range(0,10,2)))
print(list(range(0,10,8)))

for s in "abcd":
    print(s,end=' ')

print()

print('aa' in 'baa')

print("abcdfrg"[-3:])

print("abcdfrg"[::-1])
print("abcdfrg"[-1::-1])

print('aaaa'.count('a'))

print(list(range(2,5)))