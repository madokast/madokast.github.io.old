import numpy as np 
x=np.linspace(0,10,11)
y = 3*x**2+2*x+1
a=np.polyfit(x,y,2)#用2次多项式拟合x，y数组

print(a) # [3. 2. 1.]

print(type(a)) # <class 'numpy.ndarray'>