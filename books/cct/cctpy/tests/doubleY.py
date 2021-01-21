# -*- coding: utf-8 -*-
import os
import sys
curPath = os.path.abspath(os.path.dirname(__file__))
rootPath = os.path.split(curPath)[0]
PathProject = os.path.split(rootPath)[0]
sys.path.append(rootPath)
sys.path.append(PathProject)

from cctpy import *

try:
    from books.cct.cctpy.cctpy import *
except Exception as e:
    pass

R = 0.95
bl = (
    Beamline.set_start_point(start_point=P2(R, BaseUtils.angle_to_radian(-20)*R))
    .first_drift(P2.y_direct(), BaseUtils.angle_to_radian(20)*R)
    .append_agcct(
        big_r=R,
        small_rs=[128*MM + 9.5*MM, 113*MM + 9.5 *
                  MM, 98*MM + 9.5*MM, 83*MM + 9.5*MM],
        bending_angles=[17.05, 27.27, 23.18],  # [15.14, 29.02, 23.34]
        tilt_angles=[[30, 87.076, 91.829, 85.857],
                     [101.317, 30, 75.725, 92.044]],
        winding_numbers=[[128], [25, 40, 34]],
        currents=[9536.310, -6259.974],
        disperse_number_per_winding=36
    ).append_drift(BaseUtils.angle_to_radian(20)*R)
)

# ---------------------------------------------------


import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')

# 2.355T
# 17.45T/m

# 二极CCT硬边
x1 = [-20.0, 0.0, 0.0, 67.5, 67.5, 87.5]
y1 = [0.0, 0.0, 2.355, 2.355, 0.0, 0.0]

# 四极CCT硬边
g = 16.93
x2 = [-20.0, 0.0, 0.0, 12.386, 12.386, 39.730000000000004, 39.730000000000004, 67.5, 67.5, 87.5]
y2 = [0.0, 0.0, -g, -g, g, g, -g, -g, 0.0, 0.0]

# 二极CCT
bz = bl.magnetic_field_bz_along(step=100*MM)
x3=numpy.array(BaseUtils.radian_to_angle(BaseUtils.list_multiply(P2.extract_x(bz),1/R)))-20
y3=P2.extract_y(bz)

# 四极
g = bl.graident_field_along(step=100*MM)
x4=numpy.array(BaseUtils.radian_to_angle(BaseUtils.list_multiply(P2.extract_x(bz),1/R)))-20
y4=P2.extract_y(g)

fig = plt.figure()
ax = fig.add_subplot(111)

lns1 = ax.plot(x1, y1, 'g-',linewidth=3, label = 'dipole CCT (hard edge)')
lns2 = ax.plot(x3, y3, 'g-',linewidth=3, label = 'dipole CCT')
ax2 = ax.twinx()
lns3 = ax2.plot(x2, y2, 'b--',linewidth=3, label = 'AG-CCT (hard edge)')
lns4 = ax2.plot(x4, y4, 'b--',linewidth=3, label = 'AG-CCT')

lns = lns1+lns3+lns2+lns4
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, prop={'size': 18},loc='lower right')

ax.grid()

ax.yaxis.label.set_color('g')
ax.tick_params(axis='y', colors='g')

ax2.yaxis.label.set_color('b')
ax2.tick_params(axis='y', colors='b')

ax.set_xlabel("bending degree(°)", fontdict={'size':18})
ax.set_ylabel("dipole field(T)", fontdict={'size':18})
ax2.set_ylabel("quadrupole gradient(T/m)", fontdict={'size':18})


ax.tick_params(labelsize=18)
ax2.tick_params(labelsize=18)

ax.set_ylim(-4,4)
ax2.set_ylim(-20, 20)
ax.set_xlim(-20,67.5+20)

plt.show()