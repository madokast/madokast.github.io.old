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

#########################################
for m in bl.magnets:
    if isinstance(m,CCT):
        # print(m.small_r,m.tilt_angles)
        print(m.conductor_length(16))