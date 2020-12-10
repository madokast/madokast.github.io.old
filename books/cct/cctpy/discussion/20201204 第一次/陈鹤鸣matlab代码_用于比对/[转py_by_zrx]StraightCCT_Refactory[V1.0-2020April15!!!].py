# 将陈鹤鸣师兄的 matlab 代码改写为 python，边写边学习
# 我的标注用 zrx 标识

import numpy as np
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(12, 6), facecolor='w')
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122)

# 弯转 CCT 的设计参数

# N I1 I2 ln wn step miu

dTheta = 1.0  # 步长一度
step = dTheta * np.pi / 180  # 步长一度转为弧度制

omega_i = 0.01  # CCT的匝间距，unit: m
r1 = 0.06  # 内层CCT/Layer1 孔径半径, unit: m
r2 = 0.08  # 外层CCT/Layer2 孔径半径, unit: m
alpha1 = 20.0  # 内层CCT/Layer1 倾斜角 Alpha1, unit: deg.
alpha2 = 180-alpha1  # 外层CCT/Layer1 倾斜角 Alpha2=180-Alpha1, unit: deg.
alpha1 = alpha1*np.pi/180
alpha2 = alpha2*np.pi/180  # 转为弧度制

N = 8  # CCT的匝数
I1 = 12000  # CCT/Layer1, Current I1
I2 = -I1  # CCT/Layer2, Current I2; Default I2=-I1

# -----------多线模型的设计参数-------------
length = 0.006  # 线圈横截面的长度(m)
width = 0.004  # 线圈横截面的宽度(m)
ln = 1  # 线圈横截面长度方向线的数量
wn = 1  # 线圈横截面宽度方向线的数量


delta_w = width/2/wn  # 线圈横截面宽度方向上两线间距离的一半 —— zrx:这个定义很精妙，不饶真看不懂
delta_l = length/2/ln  # 线圈横截面长度方向上两线间距离的一半
I1 = I1/ln/wn  # 多线模型每根线上分摊的电流值
I2 = I2/ln/wn

miu = 4*np.pi*1e-7  # 磁场计算中的 μ
total_steps = int(2.0*np.pi*N/step)  # 总步数

# Calculate Straight CCT Path, With Two Layers!
# Allocate Memory to Store CCT Sections

# zrx：盲猜这应该是路径方程，中心线和每个导线的路径方程
# 猜错了， PathCenterVector_1_b 是副法线

PathCenterVector_1_b = np.zeros((total_steps+1, 3))
# define array to contain multiple lines (ln*wn) in CCT path
MultiLines_1 = np.zeros((total_steps+1, ln, wn, 3))

PathCenterVector_2_b = np.zeros((total_steps+1, 3))
# define array to contain multiple lines (ln*wn) in CCT path
MultiLines_2 = np.zeros((total_steps+1, ln, wn, 3))


# zrx：路径离散点
for theta in np.linspace(0, N*2*np.pi, total_steps+1):  # 这里必须+1
    # zrx 这里原本是 i=int32(theta/step)+1; 我觉得需要加一个四舍五入，防止 0.9999 -> 0
    i = int(np.round(theta/step))+1

    # 副法线反方向，来自博士论文 14 页，式子 (2.3)
    # bθ = p'(θ) * θ^ - r * z^
    # 其中 θ^ = d(r^)/dθ
    # 而 r^ = (cosθ, sinθ, 0)
    # 则 θ^ = (-sinθ, cosθ ,0)
    # 另外 p(θ) = rcot(α)sin(θ) + wθ/2pi
    # 因此 p'(θ) = rcot(α)cos(θ) + w/2pi
    PathCenterVector_1_b[i-1, :] = np.array([
        -r1*(1/np.tan(alpha1))*np.sin(theta) *
        np.cos(theta)-omega_i*np.sin(theta)/(2*np.pi),
        r1*(1/np.tan(alpha1))*np.cos(theta) *
        np.cos(theta)+omega_i*np.cos(theta)/(2*np.pi),
        -r1,
    ])
    # 归一化
    PathCenterVector_1_b[i-1, :] = PathCenterVector_1_b[i-1, :] / \
        np.linalg.norm(PathCenterVector_1_b[i-1, :])
    for m in range(1, ln+1):
        # 多线模型每根线在径向上的坐标——zrx：验证正确
        radius_1 = r1-length/2+(2*m-1)*delta_l
        for k in range(1, wn+1):
            # 多线模型每根线在副法向上的坐标——zrx：验证正确
            delta_p = (2*(k-wn/2)-1)*delta_w

            MultiLines_1[i-1, m-1, k-1, :] = np.array([
                radius_1*np.cos(theta)+delta_p*PathCenterVector_1_b[i-1, 0],
                radius_1*np.sin(theta)+delta_p*PathCenterVector_1_b[i-1, 1],
                r1*(1/np.tan(alpha1))*np.sin(theta) +
                omega_i*theta/2/np.pi+delta_p*PathCenterVector_1_b[i-1, 2],
            ])


for m in range(1, ln+1):
    for k in range(1, wn+1):
        ax1.plot(
            MultiLines_1[:, m-1, k-1, 0],
            MultiLines_1[:, m-1, k-1, 1],
            MultiLines_1[:, m-1, k-1, 2]
        )


# zrx：外层 CCT
for theta in np.linspace(0, N*2*np.pi, total_steps+1):
    i = int(np.round(theta/step))+1

    PathCenterVector_2_b[i-1, :] = np.array([
        -r2*(1/np.tan(alpha2))*np.sin(theta) *
        np.cos(theta)-omega_i*np.sin(theta)/(2*np.pi),
        r2*(1/np.tan(alpha2))*np.cos(theta) *
        np.cos(theta)+omega_i*np.cos(theta)/(2*np.pi),
        -r2,
    ])
    # 归一化
    PathCenterVector_2_b[i-1, :] = PathCenterVector_2_b[i-1, :] / \
        np.linalg.norm(PathCenterVector_2_b[i-1, :])
    for m in range(1, ln+1):
        radius_2 = r2-length/2+(2*m-1)*delta_l
        for k in range(1, wn+1):
            delta_p = (2*(k-wn/2)-1)*delta_w

            MultiLines_2[i-1, m-1, k-1, :] = np.array([
                radius_2*np.cos(theta)+delta_p*PathCenterVector_2_b[i-1, 0],
                radius_2*np.sin(theta)+delta_p*PathCenterVector_2_b[i-1, 1],
                r2*(1/np.tan(alpha2))*np.sin(theta) +
                omega_i*theta/2/np.pi+delta_p*PathCenterVector_2_b[i-1, 2],
            ])

for m in range(1, ln+1):
    for k in range(1, wn+1):
        ax1.plot(MultiLines_2[:, m-1, k-1, 0], MultiLines_2[:,
                                                            m-1, k-1, 1], MultiLines_2[:, m-1, k-1, 2])


# 计算磁场
# 待计算直线的起始和终结点：[x, y, z] unit: m
start_point = np.array([0, 0, -0.2])
end_point = np.array([0, 0, 0.3])
# 节点个数
point_num = 51  # 所需计算直线上的节点数目
CalcuPoint = np.linspace(start_point, end_point, point_num)

dB1 = np.zeros((total_steps+1, 3))  # 一根线产生的的磁场？
sumdB1 = np.zeros((ln, wn, 3))  # 每根线产生的磁场？
B1 = np.zeros((point_num, 3))  # 每个点的磁场，最终值

dB2 = np.zeros((total_steps+1, 3))  # 外层 CCT
sumdB2 = np.zeros((ln, wn, 3))
B2 = np.zeros((point_num, 3))

B_total = np.zeros((point_num, 3))  # 内外和

for t in range(point_num):  # 注意 t 初值 0，和 matlab 不同，索引起始也不同
    print(t)
    for m in range(ln):
        for k in range(wn):
            # 注意和之前的不同，从 step 开始而非 0
            for theta in np.linspace(step, N*2*np.pi, total_steps):
                i = int(np.round(theta/step))

                # dB = (miu0/4pi) Idl×r^/r^2

                # 电流方向 dl
                Current_vector = np.array([
                    # 注意下标改变
                    MultiLines_1[i, m, k, 0] - MultiLines_1[i-1, m, k, 0],
                    MultiLines_1[i, m, k, 1] - MultiLines_1[i-1, m, k, 1],
                    MultiLines_1[i, m, k, 2] - MultiLines_1[i-1, m, k, 2],
                ])

                # 电流元位置（中心） p0
                Current_Center_Point = np.array([
                    MultiLines_1[i, m, k, 0] + MultiLines_1[i-1, m, k, 0],
                    MultiLines_1[i, m, k, 1] + MultiLines_1[i-1, m, k, 1],
                    MultiLines_1[i, m, k, 2] + MultiLines_1[i-1, m, k, 2],
                ])/2.0

                # 电流元到 p，即 r
                CurrentToPoint_vector = CalcuPoint[t, :] - Current_Center_Point
                # r 归一化 er = r^
                radial_vector_norm = np.linalg.norm(CurrentToPoint_vector)
                er = CurrentToPoint_vector / radial_vector_norm
                # dl×r^
                multi = np.cross(Current_vector, er)
                # (miu0/4pi) Idl×r^/r^2
                dB1[i, :] = miu/4/np.pi * \
                    (I1*multi/radial_vector_norm/radial_vector_norm)

                # 下面计算外层 CCT
                Current_vector = np.array([
                    # 注意下标改变
                    MultiLines_2[i, m, k, 0] - MultiLines_2[i-1, m, k, 0],
                    MultiLines_2[i, m, k, 1] - MultiLines_2[i-1, m, k, 1],
                    MultiLines_2[i, m, k, 2] - MultiLines_2[i-1, m, k, 2],
                ])
                Current_Center_Point = np.array([
                    MultiLines_2[i, m, k, 0] + MultiLines_2[i-1, m, k, 0],
                    MultiLines_2[i, m, k, 1] + MultiLines_2[i-1, m, k, 1],
                    MultiLines_2[i, m, k, 2] + MultiLines_2[i-1, m, k, 2],
                ])/2.0
                CurrentToPoint_vector = CalcuPoint[t, :] - Current_Center_Point
                radial_vector_norm = np.linalg.norm(CurrentToPoint_vector)
                er = CurrentToPoint_vector / radial_vector_norm
                multi = np.cross(Current_vector, er)
                dB2[i, :] = miu/4/np.pi * \
                    (I2*multi/radial_vector_norm/radial_vector_norm)

            sumdB1[k, m, :] = np.sum(dB1, 0)
            sumdB2[k, m, :] = np.sum(dB2, 0)  # 每根线的离散点磁场求和

    # 所有线加起来
    B1[t, :] = np.array([
        np.sum(sumdB1[:, :, 0]),
        np.sum(sumdB1[:, :, 1]),
        np.sum(sumdB1[:, :, 2])
    ])

    B2[t, :] = np.array([
        np.sum(sumdB2[:, :, 0]),
        np.sum(sumdB2[:, :, 1]),
        np.sum(sumdB2[:, :, 2])
    ])

    B_total[t, :] = np.array([
        B1[t, 0]+B2[t, 0],
        B1[t, 1]+B2[t, 1],
        B1[t, 2]+B2[t, 2]
    ])


# 画图
# 起点和终点距离
line_length = end_point-start_point
line_length = np.linalg.norm(line_length)
x_length = np.linspace(0, line_length, point_num)

ax2.plot(x_length, B_total[:, 0], 'r-')
ax2.plot(x_length, B_total[:, 1], 'y-')
ax2.plot(x_length, B_total[:, 2], 'g-')
plt.legend(['Bx', 'By', 'Bz'])

plt.show()
