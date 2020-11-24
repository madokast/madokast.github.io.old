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
# 圆柱轴向分段
z_number= 50

# 定义三维数据
ksi_steps = np.linspace(0, 2 * np.pi, ksi_number)

xx = r * np.cos(ksi_steps)
zz = r * np.sin(ksi_steps)
yy = np.linspace(0, length, z_number)
X, Y = np.meshgrid(xx, yy)
Z, _ = np.meshgrid(zz, yy)


# 作图
ax3.plot_surface(X, Y, Z, color="r", alpha=0.8)
ax3.plot_wireframe(X, Y, Z, rstride=10, cstride=125,color='grey')

ax3.plot3D(-2,0,-2)
ax3.plot3D(2,0,2)


plt.axis('off')
ax3.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
plt.show()