from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

fig = plt.figure()  # 定义新的三维坐标轴
ax3 = plt.axes(projection="3d")

# 圆柱半径
r = 1.0
# 圆柱长度
length = 5.0
# 圆柱圆周方向分段
ksi_number = 500
# 圆柱轴向分段 注意 Y 方向是圆柱轴向
z_number = 500

# 定义三维数据
ksi_steps = np.linspace(0, 2 * np.pi, ksi_number)
# 圆柱数据
xx = r * np.cos(ksi_steps)
zz = r * np.sin(ksi_steps)
yy = np.linspace(0, length, z_number)
X, Y = np.meshgrid(xx, yy)
Z, _ = np.meshgrid(zz, yy)

# 圆柱作图
ax3.plot_surface(X, Y, Z, color="blue", alpha=0.1)
ax3.plot_wireframe(X, Y, Z, rstride=100, cstride=125, color="grey",alpha=0.1)

# 扩大图形范围，让圆柱长一点
ax3.plot3D(-2, 0, -2)
ax3.plot3D(2, 0, 2)

# 螺线管定义
ksi_steps = np.linspace(0, 2 * 2 * np.pi, 1000)
x = r * np.cos(ksi_steps)
z = r * np.sin(ksi_steps)
y = 2.0* ksi_steps / (2 * np.pi)
ax3.plot3D(x,y,z,color='r',alpha=0.8,lw=2)

# 坐标
ax3.plot3D([0,2],[0,0],[0,0],color='k')
ax3.plot3D([0,0],[0,0],[0,2],color='k')
ax3.plot3D([0,0],[0,8],[0,0],color='k')

plt.axis("off")
ax3.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
plt.show()