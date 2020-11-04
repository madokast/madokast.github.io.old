names = []

print(names)

print(bool(names))

names = ["123", "abd"]

print("123" in names)

names.insert(1, "a")

print(names)

names[names.index("a")] = "b"

print(names)

del names[1]

print(names)

# ------------  切片

arr = list(range(5))

print(arr)

print(arr[2:4])

sub = arr[2:4]

arr[2] = 100

print(arr)

print(sub)

print(arr.__class__)

# ------------  add

arr = []

for i in range(5):
    arr.append(i)

arr.extend("abc")

arr.extend([1, 2, 3, 4, 5])

arr = arr * 2

print(arr)

# -------------- 头加

arr = []

for i in range(-5, 6):
    arr = [i] + arr

print(arr)

arr.sort()

print(arr)

print("aAaaBB".swapcase())

# ------------ tuple


def printType(obj):
    print(type(obj))

"""
<class 'tuple'>
<class 'int'>
<class 'tuple'>
<class 'tuple'>
"""

printType(())
printType((1))
printType((1,))
printType((1, 2))

arr = list(range(5))

t = tuple(arr)

print(arr)

print(t[:3])

a,*x,b = list(range(5))

printType(a)
printType(x)
printType(b)

print(a)
print(x)
print(b)
print(*x)

print(*'zrx')

# SyntaxError: can't use starred expression here
# t = *[1,2,3]
# print(t)
# printType(t)

print(*'123')
print(*[1,2,3])
print(*(1,2,3))

m = {}
print(m)

m = {'name':'mdk'}
print(m)

m = dict([('name','mdk'),('age',12)])
print(m)

m['id']=123
print(m)

m['id']=111
print(m)

for k in m:
    print(k,m[k])

print(m.items())
print(m.values())
print(m.keys())

for k,v in m.items():
    print(k,v)

print(m)
print(m.pop('name'))
print(m)

arr = [1,1,1,2,3]
arr = list(set(arr))
print(arr)

s = {1,2,3}
s.remove(1)
print(s)

s.pop()
print(s)

printType(s)
s = frozenset(s)
printType(s)

m = {}
print(type(m))