import numpy as np

a = np.random.randn(56016*3)

b =  np.pad(a,(0,(2000*120-56016)*3),mode='constant')

print(a)
print(b)

print(a.shape)
print(b.shape)