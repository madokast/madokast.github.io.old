import numpy as np

tb = np.array(
    [
        [1, 2, 3, 4],
        [1.1, 2.1, 3.1, 4.1],
        [1.2, 2.3, 3.3, 4.3],
    ],
    dtype=np.float64,
)

print(tb)
print("------------------------")

x0 = tb[:, [0]]
x1 = tb[:, [1]]
x2 = tb[:, [2]]
x3 = tb[:, [3]]

print(x0, x1, x2, x3, sep="\n\n")
print("------------------------")
f1 = x1 + x2
f2 = x2 + x3

print(f1, f2, sep="\n\n")

print("------------------------")

print(np.hstack([f1, f2]))

print(np.column_stack((f1, f2)))
