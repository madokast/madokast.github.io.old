"""
2020年11月24日 FIXED 禁止修改
用于 books\cct\CCT几何分析并解决rib宽度问题.md
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib.animation as animation

fig = plt.figure(figsize=(7, 4.8))
# 绝对坐标系，确定绘图板位置
# ax = plt.axes([0, 0, 0.35, 0.80])
ax3 = plt.axes([0.1, -0.4, 0.6*1.5, 1.3*1.5], projection="3d")

# 圆柱半径
r = 1.0
# 圆柱长度
length = 3.0
# 圆柱圆周方向分段
ksi_number = 500
# 圆柱轴向分段 注意 Y 方向是圆柱轴向
z_number = 500

# 螺线
number = 720
theta_steps = np.linspace(0, 2 * 2 * np.pi, number)
z0 = 1.0
# 二维螺线
x2 = theta_steps
y2 = z0 / (2 * np.pi) * theta_steps
# 三维螺线
x3 = r * np.cos(theta_steps)
z3 = r * np.sin(theta_steps)
y3 = y2

# ax.plot(x2, y2, "r-")
ax3.plot(x3[0 : 90 + 40], y3[0 : 90 + 40], z3[0 : 90 + 40], "r-")
ax3.plot(x3[90:300], y3[90:300], z3[90:300], "r--")
ax3.plot(x3[300 : 450 + 40], y3[300 : 450 + 40], z3[300 : 450 + 40], "r-")
ax3.plot(x3[450:660], y3[450:660], z3[450:660], "r--")
ax3.plot(x3[660:], y3[660:], z3[660:], "r-")


# 动画绘制
if True:
    # (point_ani_2d,) = ax.plot(np.array([0, x2[50]]), np.array([y2[50], y2[50]]), "b.")
    (point_ani_3d,) = ax3.plot(
        np.array([0, x3[0]]), np.array([y3[0], y3[0]]), np.array([0, z3[0]]), "b."
    )
    # (point_ani_2d_z_projecct,) = ax.plot(
    #     np.array([0, x2[50]]), np.array([y2[50], y2[50]]), "b"
    # )
    (point_ani_3d_z_projecct,) = ax3.plot(
        np.array([0, x3[0]]), np.array([y3[0], y3[0]]), np.array([0, z3[0]]), "b"
    )

    # def update_2d(index):
    #     point_ani_2d.set_data(
    #         np.array([0, x2[index]]), np.array([y2[index], y2[index]])
    #     )
    #     point_ani_2d_z_projecct.set_data(
    #         np.array([0, x2[index]]), np.array([y2[index], y2[index]])
    #     )

    def update_3d(index):
        point_ani_3d.set_data(
            np.array([0, x3[index]]), np.array([[y3[index], y3[index]]])
        )
        point_ani_3d.set_3d_properties(np.array([0, z3[index]]))
        point_ani_3d_z_projecct.set_data(
            np.array([0, x3[index]]), np.array([[y3[index], y3[index]]])
        )
        point_ani_3d_z_projecct.set_3d_properties(np.array([0, z3[index]]))

    def update(index):
        # update_2d(index)
        update_3d(index)
        return (
            # point_ani_2d,
            point_ani_3d,
            # point_ani_2d_z_projecct,
            point_ani_3d_z_projecct,
        )

    ani = animation.FuncAnimation(
        fig, update, np.arange(0, number, 2), interval=10, blit=True
    )


# 二维面板设置
# if True:
#     ax.spines["bottom"].set_position(("data", 0))
#     ax.spines["left"].set_position(("data", 0))
#     ax.spines["right"].set_color("none")
#     ax.spines["top"].set_color("none")

#     ax.set_xlim(-2.5 * 2 * np.pi, 2.5 * 2 * np.pi)
#     ax.set_ylim(-4, 4)

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

    # 文字
    ax3.text(0,3,0,'z',(0,1,0),fontsize=15)
    ax3.text(2*r,0,0,'y',(1,0,0),fontsize=15)
    ax3.text(0,0,r*2,'x',(0,0,0),fontsize=15)


# 保存动画
ani.save("books\cct\img\C02三维_螺线_动画_z轴投影.gif", writer="imagemagick", fps=100)
plt.show()