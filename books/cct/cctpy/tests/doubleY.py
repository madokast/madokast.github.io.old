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
    Beamline.set_start_point(start_point=P2(
        R, BaseUtils.angle_to_radian(-20)*R))
    .first_drift(P2.y_direct(), BaseUtils.angle_to_radian(20)*R)
    .append_agcct(
        big_r=R,
        small_rs=[140.5*MM + 0.1*MM, 124.5*MM+ 0.1*MM, 108.5*MM+ 0.1*MM, 92.5*MM+ 0.1*MM],
        bending_angles=[-67.5*25/(25+40+34), -67.5*40/(25+40+34), -67.5*34/(25+40+34)],  # [15.14, 29.02, 23.34]
        tilt_angles=[[30.0, 88.8, 98.1, 91.7],
                        [101.8, 30.0, 62.7, 89.7]],
        winding_numbers=[[128], [25, 40, 34]],
        currents=[9409.261,	-7107.359],
        disperse_number_per_winding=36
    )
    .append_drift(BaseUtils.angle_to_radian(20)*R)
)

if False:

    Plot3.plot_p3(P3(),'w.')
    ps = CCT.as_cct(bl.magnets[0]).dispersed_path3_points
    Plot3.ax.plot([p.x for p in ps], [p.y for p in ps], [p.z for p in ps], color='b', linewidth=0.2)

    ps = CCT.as_cct(bl.magnets[1]).dispersed_path3_points
    Plot3.ax.plot([p.x for p in ps], [p.y for p in ps], [p.z for p in ps], color='b', linewidth=0.2)

    ps = CCT.as_cct(bl.magnets[2]).dispersed_path3_points
    Plot3.ax.plot([p.x for p in ps], [p.y for p in ps], [p.z for p in ps], color='r', linewidth=0.2)

    ps = CCT.as_cct(bl.magnets[3]).dispersed_path3_points
    Plot3.ax.plot([p.x for p in ps], [p.y for p in ps], [p.z for p in ps], color='r', linewidth=0.2)

    ps = CCT.as_cct(bl.magnets[4]).dispersed_path3_points
    Plot3.ax.plot([p.x for p in ps], [p.y for p in ps], [p.z for p in ps], color='g', linewidth=0.2)

    ps = CCT.as_cct(bl.magnets[5]).dispersed_path3_points
    Plot3.ax.plot([p.x for p in ps], [p.y for p in ps], [p.z for p in ps], color='g', linewidth=0.2)

    ps = CCT.as_cct(bl.magnets[6]).dispersed_path3_points
    Plot3.ax.plot([p.x for p in ps], [p.y for p in ps], [p.z for p in ps], color='r', linewidth=0.2)

    ps = CCT.as_cct(bl.magnets[7]).dispersed_path3_points
    Plot3.ax.plot([p.x for p in ps], [p.y for p in ps], [p.z for p in ps], color='r', linewidth=0.2)
    Plot3.remove_background_color()
    Plot3.off_axis()
    Plot3.show()

    os.exit()
# ---------------------------------------------------


import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')

# 2.355T
# 17.45T/m

fontsize=36

# 二极CCT硬边
x1 = [-20.0, 0.0, 0.0, 67.5, 67.5, 87.5]
y1 = [0.0, 0.0, 2.355, 2.355, 0.0, 0.0]

# 四极CCT硬边
g = -14.54
cct1_ang = 3.7164+8
cct2_ang = 19.93897+8
cct3_ang = 19.844626+8
x2 = [-20.0, 0.0, 0.0, cct1_ang, cct1_ang, cct1_ang+cct2_ang, cct1_ang+cct2_ang, 67.5, 67.5, 87.5]
y2 = [0.0, 0.0, -g, -g, g, g, -g, -g, 0.0, 0.0]

# 二极CCT
bz = bl.magnetic_field_bz_along(step=10*MM)
x3=numpy.array(BaseUtils.radian_to_angle(BaseUtils.list_multiply(P2.extract_x(bz),1/R)))-20
y3=P2.extract_y(bz)

# 四极
g = bl.graident_field_along(step=10*MM)
x4=numpy.array(BaseUtils.radian_to_angle(BaseUtils.list_multiply(P2.extract_x(bz),1/R)))-20
y4=P2.extract_y(g)

fig = plt.figure()
ax = fig.add_subplot(111)

lns1 = ax.plot(x1, y1, 'g--',linewidth=3, label = 'dipole CCT (SCOFF)')
lns2 = ax.plot(x3, y3, 'g-',linewidth=3, label = 'dipole CCT')
ax2 = ax.twinx()
lns3 = ax2.plot(x2, y2, 'b--',linewidth=3, label = 'AG-CCT (SCOFF)')
lns4 = ax2.plot(x4, y4, 'b-',linewidth=3, label = 'AG-CCT')

lns = lns1+lns3+lns2+lns4
# lns = lns3+lns4
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, prop={'size': 25},loc='lower right')

ax.grid()

ax.yaxis.label.set_color('g')
ax.tick_params(axis='y', colors='g')

ax2.yaxis.label.set_color('b')
ax2.tick_params(axis='y', colors='b')

ax.set_xlabel("bending degree(°)", fontdict={'size':fontsize})
ax.set_ylabel("dipole field(T)", fontdict={'size':fontsize})
ax2.set_ylabel("quadrupole gradient(T/m)", fontdict={'size':fontsize})


ax.tick_params(labelsize=fontsize)
ax2.tick_params(labelsize=fontsize)

ax.set_ylim(-4,4)
ax2.set_ylim(-20, 20)
ax.set_xlim(-20,67.5+20)

plt.show()