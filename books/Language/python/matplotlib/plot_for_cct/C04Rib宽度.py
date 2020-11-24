"""
2020年11月24日 FIXED 禁止修改
用于 books\cct\CCT几何分析并解决rib宽度问题.md
"""

from typing import Callable
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib.animation as animation


def diff(fun: Callable, delta=1e-4) -> Callable:
    """
    求导数
    """
    return lambda x: (fun(x + delta) - fun(x)) * (1 / delta)


fig = plt.figure(figsize=(7, 4.8))
# 绝对坐标系，确定绘图板位置
# ax = plt.axes([0, 0, 0.35, 0.80])
ax3 = plt.axes([0.1, -0.4, 0.6 * 1.5, 1.3 * 1.5], projection="3d")

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

# 开槽
width = 0.2
depth = 0.3

# 二维螺线
x2 = theta_steps
y2 = z0 / (2 * np.pi) * theta_steps
# 三维螺线
path3_fun = lambda t: np.array([r * np.cos(t), r * np.sin(t), z0 / (2 * np.pi) * t])
path3_tangential = diff(path3_fun)  # 切向
project_z = lambda t: np.array([0.0, 0.0, z0 / (2 * np.pi) * t])
main_norm = lambda t: (path3_fun(t) - project_z(t)) / np.linalg.norm(
    path3_fun(t) - project_z(t), ord=2
)
second_norm = lambda t: np.cross(path3_tangential(t), main_norm(t)) / np.linalg.norm(
    np.cross(path3_tangential(t), main_norm(t)), ord=2
)

x3 = np.array([path3_fun(t)[0] for t in theta_steps])  # r * np.cos(theta_steps)
z3 = np.array([path3_fun(t)[1] for t in theta_steps])  # r * np.sin(theta_steps)
y3 = np.array([path3_fun(t)[2] for t in theta_steps])

# ax.plot(x2, y2, "r-")
ax3.plot(x3[0 : 90 + 40], y3[0 : 90 + 40], z3[0 : 90 + 40], "r-")
ax3.plot(x3[90:300], y3[90:300], z3[90:300], "r--")
ax3.plot(x3[300 : 450 + 40], y3[300 : 450 + 40], z3[300 : 450 + 40], "r-")
ax3.plot(x3[450:660], y3[450:660], z3[450:660], "r--")
ax3.plot(x3[660:], y3[660:], z3[660:], "r-")


# 动画绘制
if True:
    (point_ani_3d,) = ax3.plot(
        np.array([0, x3[0]]), np.array([y3[0], y3[0]]), np.array([0, z3[0]]), "b."
    )
    (point_ani_3d_z_projecct,) = ax3.plot(
        np.array([0, x3[0]]), np.array([y3[0], y3[0]]), np.array([0, z3[0]]), "b"
    )
    (point_ani_3d_box,) = ax3.plot(
        np.array([
            x3[0] + main_norm(theta_steps[0])[0] * depth / 2 + second_norm(theta_steps[0])[0] * width / 2,
            x3[0] - main_norm(theta_steps[0])[0] * depth / 2 + second_norm(theta_steps[0])[0] * width / 2,
            x3[0] - main_norm(theta_steps[0])[0] * depth / 2 - second_norm(theta_steps[0])[0] * width / 2,
            x3[0] + main_norm(theta_steps[0])[0] * depth / 2 - second_norm(theta_steps[0])[0] * width / 2,
            x3[0] + main_norm(theta_steps[0])[0] * depth / 2 + second_norm(theta_steps[0])[0] * width / 2,
            ]),
        np.array([
            y3[0] + main_norm(theta_steps[0])[2] * depth / 2 + second_norm(theta_steps[0])[2] * width / 2,
            y3[0] - main_norm(theta_steps[0])[2] * depth / 2 + second_norm(theta_steps[0])[2] * width / 2,
            y3[0] - main_norm(theta_steps[0])[2] * depth / 2 - second_norm(theta_steps[0])[2] * width / 2,
            y3[0] + main_norm(theta_steps[0])[2] * depth / 2 - second_norm(theta_steps[0])[2] * width / 2,
            y3[0] + main_norm(theta_steps[0])[2] * depth / 2 + second_norm(theta_steps[0])[2] * width / 2,
        ]),
        np.array([
            z3[0] + main_norm(theta_steps[0])[1] * depth / 2 + second_norm(theta_steps[0])[1] * width / 2,
            z3[0] - main_norm(theta_steps[0])[1] * depth / 2 + second_norm(theta_steps[0])[1] * width / 2,
            z3[0] - main_norm(theta_steps[0])[1] * depth / 2 - second_norm(theta_steps[0])[1] * width / 2,
            z3[0] + main_norm(theta_steps[0])[1] * depth / 2 - second_norm(theta_steps[0])[1] * width / 2,
            z3[0] + main_norm(theta_steps[0])[1] * depth / 2 + second_norm(theta_steps[0])[1] * width / 2,
        ]),
        "b-",
    )

    (point_ani_3d_box_1,) = ax3.plot(
        np.array([
            x3[0] + main_norm(theta_steps[358])[0] * depth / 2 + second_norm(theta_steps[358])[0] * width / 2,
            x3[0] - main_norm(theta_steps[358])[0] * depth / 2 + second_norm(theta_steps[358])[0] * width / 2,
            x3[0] - main_norm(theta_steps[358])[0] * depth / 2 - second_norm(theta_steps[358])[0] * width / 2,
            x3[0] + main_norm(theta_steps[358])[0] * depth / 2 - second_norm(theta_steps[358])[0] * width / 2,
            x3[0] + main_norm(theta_steps[358])[0] * depth / 2 + second_norm(theta_steps[358])[0] * width / 2,
            ]),
        np.array([
            y3[0] + main_norm(theta_steps[358])[2] * depth / 2 + second_norm(theta_steps[358])[2] * width / 2,
            y3[0] - main_norm(theta_steps[358])[2] * depth / 2 + second_norm(theta_steps[358])[2] * width / 2,
            y3[0] - main_norm(theta_steps[358])[2] * depth / 2 - second_norm(theta_steps[358])[2] * width / 2,
            y3[0] + main_norm(theta_steps[358])[2] * depth / 2 - second_norm(theta_steps[358])[2] * width / 2,
            y3[0] + main_norm(theta_steps[358])[2] * depth / 2 + second_norm(theta_steps[358])[2] * width / 2,
        ]),
        np.array([
            z3[0] + main_norm(theta_steps[358])[1] * depth / 2 + second_norm(theta_steps[358])[1] * width / 2,
            z3[0] - main_norm(theta_steps[358])[1] * depth / 2 + second_norm(theta_steps[358])[1] * width / 2,
            z3[0] - main_norm(theta_steps[358])[1] * depth / 2 - second_norm(theta_steps[358])[1] * width / 2,
            z3[0] + main_norm(theta_steps[358])[1] * depth / 2 - second_norm(theta_steps[358])[1] * width / 2,
            z3[0] + main_norm(theta_steps[358])[1] * depth / 2 + second_norm(theta_steps[358])[1] * width / 2,
        ]),
        "b-",
    )

    (point_ani_3d_box_2,) = ax3.plot(
        np.array([
            x3[0] + main_norm(theta_steps[360])[0] * depth / 2 + second_norm(theta_steps[360])[0] * width / 2,
            x3[0] - main_norm(theta_steps[360])[0] * depth / 2 + second_norm(theta_steps[360])[0] * width / 2,
            x3[0] - main_norm(theta_steps[360])[0] * depth / 2 - second_norm(theta_steps[360])[0] * width / 2,
            x3[0] + main_norm(theta_steps[360])[0] * depth / 2 - second_norm(theta_steps[360])[0] * width / 2,
            x3[0] + main_norm(theta_steps[360])[0] * depth / 2 + second_norm(theta_steps[360])[0] * width / 2,
            ]),
        np.array([
            y3[0] + main_norm(theta_steps[360])[2] * depth / 2 + second_norm(theta_steps[360])[2] * width / 2,
            y3[0] - main_norm(theta_steps[360])[2] * depth / 2 + second_norm(theta_steps[360])[2] * width / 2,
            y3[0] - main_norm(theta_steps[360])[2] * depth / 2 - second_norm(theta_steps[360])[2] * width / 2,
            y3[0] + main_norm(theta_steps[360])[2] * depth / 2 - second_norm(theta_steps[360])[2] * width / 2,
            y3[0] + main_norm(theta_steps[360])[2] * depth / 2 + second_norm(theta_steps[360])[2] * width / 2,
        ]),
        np.array([
            z3[0] + main_norm(theta_steps[360])[1] * depth / 2 + second_norm(theta_steps[360])[1] * width / 2,
            z3[0] - main_norm(theta_steps[360])[1] * depth / 2 + second_norm(theta_steps[360])[1] * width / 2,
            z3[0] - main_norm(theta_steps[360])[1] * depth / 2 - second_norm(theta_steps[360])[1] * width / 2,
            z3[0] + main_norm(theta_steps[360])[1] * depth / 2 - second_norm(theta_steps[360])[1] * width / 2,
            z3[0] + main_norm(theta_steps[360])[1] * depth / 2 + second_norm(theta_steps[360])[1] * width / 2,
        ]),
        "b-",
    )

    (point_ani_3d_box_3,) = ax3.plot(
        np.array([
            x3[0] + main_norm(theta_steps[362])[0] * depth / 2 + second_norm(theta_steps[362])[0] * width / 2,
            x3[0] - main_norm(theta_steps[362])[0] * depth / 2 + second_norm(theta_steps[362])[0] * width / 2,
            x3[0] - main_norm(theta_steps[362])[0] * depth / 2 - second_norm(theta_steps[362])[0] * width / 2,
            x3[0] + main_norm(theta_steps[362])[0] * depth / 2 - second_norm(theta_steps[362])[0] * width / 2,
            x3[0] + main_norm(theta_steps[362])[0] * depth / 2 + second_norm(theta_steps[362])[0] * width / 2,
            ]),
        np.array([
            y3[0] + main_norm(theta_steps[362])[2] * depth / 2 + second_norm(theta_steps[362])[2] * width / 2,
            y3[0] - main_norm(theta_steps[362])[2] * depth / 2 + second_norm(theta_steps[362])[2] * width / 2,
            y3[0] - main_norm(theta_steps[362])[2] * depth / 2 - second_norm(theta_steps[362])[2] * width / 2,
            y3[0] + main_norm(theta_steps[362])[2] * depth / 2 - second_norm(theta_steps[362])[2] * width / 2,
            y3[0] + main_norm(theta_steps[362])[2] * depth / 2 + second_norm(theta_steps[362])[2] * width / 2,
        ]),
        np.array([
            z3[0] + main_norm(theta_steps[362])[1] * depth / 2 + second_norm(theta_steps[362])[1] * width / 2,
            z3[0] - main_norm(theta_steps[362])[1] * depth / 2 + second_norm(theta_steps[362])[1] * width / 2,
            z3[0] - main_norm(theta_steps[362])[1] * depth / 2 - second_norm(theta_steps[362])[1] * width / 2,
            z3[0] + main_norm(theta_steps[362])[1] * depth / 2 - second_norm(theta_steps[362])[1] * width / 2,
            z3[0] + main_norm(theta_steps[362])[1] * depth / 2 + second_norm(theta_steps[362])[1] * width / 2,
        ]),
        "b-",
    )

    def update_3d(index):
        point_ani_3d.set_data(
            np.array([0, x3[index]]), np.array([[y3[index], y3[index]]])
        )
        point_ani_3d.set_3d_properties(np.array([0, z3[index]]))
        point_ani_3d_z_projecct.set_data(
            np.array([0, x3[index]]), np.array([[y3[index], y3[index]]])
        )
        point_ani_3d_z_projecct.set_3d_properties(np.array([0, z3[index]]))

        point_ani_3d_box.set_data(
            np.array([
                x3[index] + main_norm(theta_steps[index])[0] * depth / 2+ second_norm(theta_steps[index])[0] * width / 2,
                x3[index] - main_norm(theta_steps[index])[0] * depth / 2+ second_norm(theta_steps[index])[0] * width / 2,
                x3[index] - main_norm(theta_steps[index])[0] * depth / 2- second_norm(theta_steps[index])[0] * width / 2,
                x3[index] + main_norm(theta_steps[index])[0] * depth / 2- second_norm(theta_steps[index])[0] * width / 2,
                x3[index] + main_norm(theta_steps[index])[0] * depth / 2+ second_norm(theta_steps[index])[0] * width / 2,
            ]),
            np.array([
                y3[index] + main_norm(theta_steps[index])[2] * depth / 2+ second_norm(theta_steps[index])[2] * width / 2,
                y3[index] - main_norm(theta_steps[index])[2] * depth / 2+ second_norm(theta_steps[index])[2] * width / 2,
                y3[index] - main_norm(theta_steps[index])[2] * depth / 2- second_norm(theta_steps[index])[2] * width / 2,
                y3[index] + main_norm(theta_steps[index])[2] * depth / 2- second_norm(theta_steps[index])[2] * width / 2,
                y3[index] + main_norm(theta_steps[index])[2] * depth / 2+ second_norm(theta_steps[index])[2] * width / 2,
            ]),
        )
        point_ani_3d_box.set_3d_properties(
            np.array([
                z3[index] + main_norm(theta_steps[index])[1] * depth / 2+ second_norm(theta_steps[index])[1] * width / 2,
                z3[index] - main_norm(theta_steps[index])[1] * depth / 2+ second_norm(theta_steps[index])[1] * width / 2,
                z3[index] - main_norm(theta_steps[index])[1] * depth / 2- second_norm(theta_steps[index])[1] * width / 2,
                z3[index] + main_norm(theta_steps[index])[1] * depth / 2- second_norm(theta_steps[index])[1] * width / 2,
                z3[index] + main_norm(theta_steps[index])[1] * depth / 2+ second_norm(theta_steps[index])[1] * width / 2,
            ]),
        )

        ###########################################
        if index+362<720:
            point_ani_3d_box_1.set_data(
                np.array([
                    x3[358+index] + main_norm(theta_steps[358+index])[0] * depth / 2+ second_norm(theta_steps[358+index])[0] * width / 2,
                    x3[358+index] - main_norm(theta_steps[358+index])[0] * depth / 2+ second_norm(theta_steps[358+index])[0] * width / 2,
                    x3[358+index] - main_norm(theta_steps[358+index])[0] * depth / 2- second_norm(theta_steps[358+index])[0] * width / 2,
                    x3[358+index] + main_norm(theta_steps[358+index])[0] * depth / 2- second_norm(theta_steps[358+index])[0] * width / 2,
                    x3[358+index] + main_norm(theta_steps[358+index])[0] * depth / 2+ second_norm(theta_steps[358+index])[0] * width / 2,
                ]),
                np.array([
                    y3[358+index] + main_norm(theta_steps[358+index])[2] * depth / 2+ second_norm(theta_steps[358+index])[2] * width / 2,
                    y3[358+index] - main_norm(theta_steps[358+index])[2] * depth / 2+ second_norm(theta_steps[358+index])[2] * width / 2,
                    y3[358+index] - main_norm(theta_steps[358+index])[2] * depth / 2- second_norm(theta_steps[358+index])[2] * width / 2,
                    y3[358+index] + main_norm(theta_steps[358+index])[2] * depth / 2- second_norm(theta_steps[358+index])[2] * width / 2,
                    y3[358+index] + main_norm(theta_steps[358+index])[2] * depth / 2+ second_norm(theta_steps[358+index])[2] * width / 2,
                ]),
            )
            point_ani_3d_box_1.set_3d_properties(
                np.array([
                    z3[358+index] + main_norm(theta_steps[358+index])[1] * depth / 2+ second_norm(theta_steps[358+index])[1] * width / 2,
                    z3[358+index] - main_norm(theta_steps[358+index])[1] * depth / 2+ second_norm(theta_steps[358+index])[1] * width / 2,
                    z3[358+index] - main_norm(theta_steps[358+index])[1] * depth / 2- second_norm(theta_steps[358+index])[1] * width / 2,
                    z3[358+index] + main_norm(theta_steps[358+index])[1] * depth / 2- second_norm(theta_steps[358+index])[1] * width / 2,
                    z3[358+index] + main_norm(theta_steps[358+index])[1] * depth / 2+ second_norm(theta_steps[358+index])[1] * width / 2,
                ]),
            )

            point_ani_3d_box_2.set_data(
                np.array([
                    x3[360+index] + main_norm(theta_steps[360+index])[0] * depth / 2+ second_norm(theta_steps[360+index])[0] * width / 2,
                    x3[360+index] - main_norm(theta_steps[360+index])[0] * depth / 2+ second_norm(theta_steps[360+index])[0] * width / 2,
                    x3[360+index] - main_norm(theta_steps[360+index])[0] * depth / 2- second_norm(theta_steps[360+index])[0] * width / 2,
                    x3[360+index] + main_norm(theta_steps[360+index])[0] * depth / 2- second_norm(theta_steps[360+index])[0] * width / 2,
                    x3[360+index] + main_norm(theta_steps[360+index])[0] * depth / 2+ second_norm(theta_steps[360+index])[0] * width / 2,
                ]),
                np.array([
                    y3[360+index] + main_norm(theta_steps[360+index])[2] * depth / 2+ second_norm(theta_steps[360+index])[2] * width / 2,
                    y3[360+index] - main_norm(theta_steps[360+index])[2] * depth / 2+ second_norm(theta_steps[360+index])[2] * width / 2,
                    y3[360+index] - main_norm(theta_steps[360+index])[2] * depth / 2- second_norm(theta_steps[360+index])[2] * width / 2,
                    y3[360+index] + main_norm(theta_steps[360+index])[2] * depth / 2- second_norm(theta_steps[360+index])[2] * width / 2,
                    y3[360+index] + main_norm(theta_steps[360+index])[2] * depth / 2+ second_norm(theta_steps[360+index])[2] * width / 2,
                ]),
            )
            point_ani_3d_box_2.set_3d_properties(
                np.array([
                    z3[360+index] + main_norm(theta_steps[360+index])[1] * depth / 2+ second_norm(theta_steps[360+index])[1] * width / 2,
                    z3[360+index] - main_norm(theta_steps[360+index])[1] * depth / 2+ second_norm(theta_steps[360+index])[1] * width / 2,
                    z3[360+index] - main_norm(theta_steps[360+index])[1] * depth / 2- second_norm(theta_steps[360+index])[1] * width / 2,
                    z3[360+index] + main_norm(theta_steps[360+index])[1] * depth / 2- second_norm(theta_steps[360+index])[1] * width / 2,
                    z3[360+index] + main_norm(theta_steps[360+index])[1] * depth / 2+ second_norm(theta_steps[360+index])[1] * width / 2,
                ]),
            )

            point_ani_3d_box_3.set_data(
                np.array([
                    x3[362+index] + main_norm(theta_steps[362+index])[0] * depth / 2+ second_norm(theta_steps[362+index])[0] * width / 2,
                    x3[362+index] - main_norm(theta_steps[362+index])[0] * depth / 2+ second_norm(theta_steps[362+index])[0] * width / 2,
                    x3[362+index] - main_norm(theta_steps[362+index])[0] * depth / 2- second_norm(theta_steps[362+index])[0] * width / 2,
                    x3[362+index] + main_norm(theta_steps[362+index])[0] * depth / 2- second_norm(theta_steps[362+index])[0] * width / 2,
                    x3[362+index] + main_norm(theta_steps[362+index])[0] * depth / 2+ second_norm(theta_steps[362+index])[0] * width / 2,
                ]),
                np.array([
                    y3[362+index] + main_norm(theta_steps[362+index])[2] * depth / 2+ second_norm(theta_steps[362+index])[2] * width / 2,
                    y3[362+index] - main_norm(theta_steps[362+index])[2] * depth / 2+ second_norm(theta_steps[362+index])[2] * width / 2,
                    y3[362+index] - main_norm(theta_steps[362+index])[2] * depth / 2- second_norm(theta_steps[362+index])[2] * width / 2,
                    y3[362+index] + main_norm(theta_steps[362+index])[2] * depth / 2- second_norm(theta_steps[362+index])[2] * width / 2,
                    y3[362+index] + main_norm(theta_steps[362+index])[2] * depth / 2+ second_norm(theta_steps[362+index])[2] * width / 2,
                ]),
            )
            point_ani_3d_box_3.set_3d_properties(
                np.array([
                    z3[362+index] + main_norm(theta_steps[362+index])[1] * depth / 2+ second_norm(theta_steps[362+index])[1] * width / 2,
                    z3[362+index] - main_norm(theta_steps[362+index])[1] * depth / 2+ second_norm(theta_steps[362+index])[1] * width / 2,
                    z3[362+index] - main_norm(theta_steps[362+index])[1] * depth / 2- second_norm(theta_steps[362+index])[1] * width / 2,
                    z3[362+index] + main_norm(theta_steps[362+index])[1] * depth / 2- second_norm(theta_steps[362+index])[1] * width / 2,
                    z3[362+index] + main_norm(theta_steps[362+index])[1] * depth / 2+ second_norm(theta_steps[362+index])[1] * width / 2,
                ]),
            )
        ###########################################


    def update(index):
        update_3d(index)
        return (
            # point_ani_2d,
            point_ani_3d,
            # point_ani_2d_z_projecct,
            point_ani_3d_z_projecct,
            point_ani_3d_box,
            point_ani_3d_box_1,
            point_ani_3d_box_2,
            point_ani_3d_box_3,
        )

    ani = animation.FuncAnimation(
        fig, update, np.arange(0, number//2, 2), interval=10, blit=True
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
    ax3.text(0, 3, 0, "z", (0, 1, 0), fontsize=15)
    ax3.text(2 * r, 0, 0, "y", (1, 0, 0), fontsize=15)
    ax3.text(0, 0, r * 2, "x", (0, 0, 0), fontsize=15)


# 保存动画
ani.save("books\cct\img\C04Rib宽度.gif", writer="imagemagick", fps=100)
plt.show()