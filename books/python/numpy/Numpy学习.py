import numpy as np
import time


def duration(fun, *args) -> None:
    l = time.time()
    fun(*args)
    print(f"duration={time.time()-l}")


print("--------验证安装--------")
a1 = np.eye(4, 4)
print(a1)

a2 = np.array([1.0, 2.0, 3.0])
print(a2)

print(a2.dtype.name)


duration(lambda len: np.sin(np.linspace(0, len, 360)), 2 * np.pi)
duration(lambda len: np.array([np.sin(x / len) for x in range(0, 360)]), 2 * np.pi)


def test(len):
    r = []
    for i in range(360):
        r.append(np.sin(i / len))
    return np.array(r)


duration(test, 2 * np.pi)

print("---------------- 广播 -----------------")

a3 = np.array([1, 2, 3, 4], dtype=np.float64)
print(a3)
print(np.sum(a3))

a4 = np.array([[1, 2, 3], [2, 3, 4], [11, 22, 33]])
print(a4)
print(np.sum(a4))
print(np.sum(a4, axis=0))  # 列求和
print(np.sum(a4, axis=1))  # 行求和

print("-------------- 叉乘 ------------")
x = np.array([1, 0, 0], dtype=np.float64)
y = np.array([0, 1, 0], dtype=np.float64)
z = np.cross(x, y)
print(x, " X ", y, " = ", z)


x = np.array([1, 2, 3], dtype=np.float64)
y = np.array([4, 5, 6], dtype=np.float64)
z = np.cross(x, y)
print(x, " X ", y, " = ", z)


x = np.array([2, 2, 3], dtype=np.float64)
y = np.array([4, 5, 6], dtype=np.float64)
z = np.cross(x, y)
print(x, " X ", y, " = ", z)

# [1. 0. 0.]  X  [0. 1. 0.]  =  [0. 0. 1.]
# [1. 2. 3.]  X  [4. 5. 6.]  =  [-3.  6. -3.]
# [2. 2. 3.]  X  [4. 5. 6.]  =  [-3.  0.  2.]

xs = np.array([[1, 0, 0], [1, 2, 3], [2, 2, 3]], dtype=np.float64)
ys = np.array([[0, 1, 0], [4, 5, 6], [4, 5, 6]], dtype=np.float64)
zs = np.cross(xs, ys)
print(xs, ys, zs, sep="\n")

