"""
2020年11月25日 FIXED 禁止修改
用于 books\cct\CCT几何分析并解决rib宽度问题.md

增加坐标 x y z

gif 压缩：https://gifcompressor.com/zh/
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib.animation as animation

fig = plt.figure(figsize=(11, 4.8))
# 绝对坐标系，确定绘图板位置
ax = plt.axes([0.15, 0.10, 0.35, 0.80])
ax3 = plt.axes([0.5, -0.05, 0.6, 1.3], projection="3d")

# 圆柱半径
r = 1.0
# 圆柱长度
length = 3.0
# 圆柱圆周方向分段
ksi_number = 500
# 圆柱轴向分段 注意 Y 方向是圆柱轴向
z_number = 500

x2s=np.array([])
y2s=np.array([])
x3s=np.array([])
y3s=np.array([])
z3s=np.array([])

# CCT0
if True: 
    # CCT
    number = 10
    theta_steps = np.linspace(0 ,number * 2 * np.pi, number*360)
    z0 = 0.08
    start_z = 0.0
    c0 = 0.0
    c1 = 0.3
    # 二维螺线
    x2 = theta_steps
    y2 = (
        z0 / (2 * np.pi) * theta_steps
        + start_z
        + c1 * np.sin(2 * theta_steps)
    )
    # 三维螺线
    x3 = r * np.cos(theta_steps)
    z3 = r * np.sin(theta_steps)
    y3 = y2

    # ax.plot(x2, y2, "r-")
    # ax3.plot(x3, y3, z3, "r-")

    x2s = np.concatenate((x2s,x2))
    y2s = np.concatenate((y2s,y2))
    x3s = np.concatenate((x3s,x3))
    y3s = np.concatenate((y3s,y3))
    z3s = np.concatenate((z3s,z3))

# CCT1
if True: 
    # CCT
    number = 10
    number2 = 8
    theta_steps = np.linspace(number * 2 * np.pi ,(number-number2) * 2 * np.pi, number2*360)
    z0 = -0.08
    start_z = 0.08*number*2
    c0 = 0.0
    c1 = 0.3
    # 二维螺线
    x2 = theta_steps
    y2 = (
        z0 / (2 * np.pi) * theta_steps
        + start_z
        + c1 * np.sin(2 * theta_steps)
    )
    # 三维螺线
    x3 = r * np.cos(theta_steps)
    z3 = r * np.sin(theta_steps)
    y3 = y2

    # ax.plot(x2, y2, "b-")
    # ax3.plot(x3, y3, z3, "b-")

    x2s = np.concatenate((x2s,x2))
    y2s = np.concatenate((y2s,y2))
    x3s = np.concatenate((x3s,x3))
    y3s = np.concatenate((y3s,y3))
    z3s = np.concatenate((z3s,z3))

# CCT2
if True: 
    # CCT
    number = 10
    number2 = 8
    number3 = 15
    theta_steps = np.linspace((number-number2) * 2 * np.pi ,(number-number2+number3) * 2 * np.pi, number3*360)
    z0 = 0.08
    start_z = 0.08*16
    c0 = 0.0
    c1 = 0.3
    # 二维螺线
    x2 = theta_steps
    y2 = (
        z0 / (2 * np.pi) * theta_steps
        + start_z
        + c1 * np.sin(2 * theta_steps)
    )
    # 三维螺线
    x3 = r * np.cos(theta_steps)
    z3 = r * np.sin(theta_steps)
    y3 = y2

    # ax.plot(x2, y2, "g-")
    # ax3.plot(x3, y3, z3, "g-")

    x2s = np.concatenate((x2s,x2))
    y2s = np.concatenate((y2s,y2))
    x3s = np.concatenate((x3s,x3))
    y3s = np.concatenate((y3s,y3))
    z3s = np.concatenate((z3s,z3))


# CCT01
if True: 
    # CCT
    number = 10
    theta_steps = np.linspace(0 ,-number * 2 * np.pi, number*360)
    z0 = -0.08
    start_z = 0.0
    c0 = 0.0
    c1 = 0.3
    # 二维螺线
    x2 = theta_steps
    y2 = (
        z0 / (2 * np.pi) * theta_steps
        + start_z
        + c1 * np.sin(2 * theta_steps)
    )
    # 三维螺线
    x3 = r * np.cos(theta_steps)
    z3 = r * np.sin(theta_steps)
    y3 = y2

    # ax.plot(x2, y2, "g-")
    # ax3.plot(x3, y3, z3, "g-")

    x2s = np.concatenate((x2s,x2))
    y2s = np.concatenate((y2s,y2))
    x3s = np.concatenate((x3s,x3))
    y3s = np.concatenate((y3s,y3))
    z3s = np.concatenate((z3s,z3))

# CCT11
if True: 
    # CCT
    number = 10
    number2 = 8
    theta_steps = np.linspace(-number * 2 * np.pi ,-(number-number2) * 2 * np.pi, number2*360)
    z0 = 0.08
    start_z = 0.08*number*2
    c0 = 0.0
    c1 = 0.3
    # 二维螺线
    x2 = theta_steps
    y2 = (
        z0 / (2 * np.pi) * theta_steps
        + start_z
        + c1 * np.sin(2 * theta_steps)
    )
    # 三维螺线
    x3 = r * np.cos(theta_steps)
    z3 = r * np.sin(theta_steps)
    y3 = y2

    # ax.plot(x2, y2, "y-")
    # ax3.plot(x3, y3, z3, "y-")

    x2s = np.concatenate((x2s,x2))
    y2s = np.concatenate((y2s,y2))
    x3s = np.concatenate((x3s,x3))
    y3s = np.concatenate((y3s,y3))
    z3s = np.concatenate((z3s,z3))

# CCT21
if True: 
    # CCT
    number = 10
    number2 = 8
    number3 = 15
    theta_steps = np.linspace(-(number-number2) * 2 * np.pi ,-(number-number2+number3) * 2 * np.pi, number3*360)
    z0 = -0.08
    start_z = 0.08*16
    c0 = 0.0
    c1 = 0.3
    # 二维螺线
    x2 = theta_steps
    y2 = (
        z0 / (2 * np.pi) * theta_steps
        + start_z
        + c1 * np.sin(2 * theta_steps)
    )
    # 三维螺线
    x3 = r * np.cos(theta_steps)
    z3 = r * np.sin(theta_steps)
    y3 = y2

    # ax.plot(x2, y2, "r-")
    # ax3.plot(x3, y3, z3, "r-")

    x2s = np.append(x2s,x2)
    y2s = np.append(y2s,y2)
    x3s = np.append(x3s,x3)
    y3s = np.append(y3s,y3)
    z3s = np.append(z3s,z3)

(point_ani_2d,) = ax.plot(x2s[0], y2s[0], "r.",ms=0.5)
(point_ani_3d,) = ax3.plot(x3s[0], y3s[0], z3s[0], "r-")

def update_2d(index):
        point_ani_2d.set_data(x2s[:index], y2s[:index])
        return point_ani_2d

def update_3d(index):
    point_ani_3d.set_data(x3s[:index], y3s[:index])
    point_ani_3d.set_3d_properties(z3s[:index])
    return point_ani_3d

def update(index):
    update_2d(index)
    update_3d(index)
    return point_ani_2d, point_ani_3d

ani = animation.FuncAnimation(
    fig, update, np.arange(0, x2s.shape[0], 180), interval=1, blit=True
)

print(x2s.shape[0])


# 二维面板设置
if True:
    ax.spines["bottom"].set_position(("data", 0))
    ax.spines["left"].set_position(("data", 0))
    ax.spines["right"].set_color("none")
    ax.spines["top"].set_color("none")

    ax.set_xlim(-30 * np.pi, 30 * np.pi)
    ax.set_ylim(-4, 4)

    ax.text(15, -1, "θ", fontsize=15)
    ax.text(1, 4, "z", fontsize=15)

# 绘制三维圆柱 （2020年11月24日 加上视角）
if True:
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
    ax3.plot_wireframe(X, Y, Z, rstride=100, cstride=125, color="grey", alpha=0.1)

    # 扩大图形范围，让圆柱长一点
    ax3.plot3D(-2, 0, -2)
    ax3.plot3D(2, 0, 2)

    # 坐标
    ax3.plot3D([0, r * 2], [0, 0], [0, 0], color="k", lw=0.8)
    ax3.plot3D([0, 0], [0, 0], [0, r * 2], color="k", lw=0.8)
    ax3.plot3D([0, 0], [0, length], [0, 0], color="k", lw=0.8)

    plt.axis("off")
    ax3.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    # 设置摄像机位置（实际上是设置坐标轴范围）
    ax3.set_xlim3d(-1.5, 1.5)
    ax3.set_ylim3d(0.5, 3.5)
    ax3.set_zlim3d(-1.2, 1.2)

    # 设置摄像机方位角
    ax3.view_init(elev=43, azim=-26)

    # 文字
    ax3.text(0, 3, 0, "z", (0, 1, 0), fontsize=15)
    ax3.text(2 * r, 0, 0, "y", (1, 0, 0), fontsize=15)
    ax3.text(0, 0, r * 2, "x", (0, 0, 0), fontsize=15)


# ani.save('books\cct\img\F01左二维有圆柱三维_AGCCT设计.gif', writer='imagemagick', fps=100)
plt.show()