import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

fig = plt.figure()  # 定义新的三维坐标轴
ax = plt.axes()

# 绘制横纵坐标
ax.spines['bottom'].set_position(('data', 0))
ax.spines['left'].set_position(('data',0))
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.set_xlim(-2,2)
ax.set_ylim(-2,2)


plt.show()