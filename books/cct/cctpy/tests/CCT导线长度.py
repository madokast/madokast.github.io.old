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
        small_rs=[140.5*MM, 124.5*MM, 108.5*MM, 92.5*MM],
        bending_angles=[17.05, 27.27, 23.18],  # [15.14, 29.02, 23.34]
        tilt_angles=[[30, 88.8, 98.1, 91.7],
                     [101.8, 30, 62.7, 89.7]],
        winding_numbers=[[128], [25, 40, 34]],
        currents=[9536.310, -6259.974],
        disperse_number_per_winding=36
    ).append_drift(BaseUtils.angle_to_radian(20)*R)
)

#########################################
for m in bl.magnets:
    if isinstance(m,CCT):
        # print(m.small_r,m.tilt_angles)
        print(m.conductor_length(16))