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

# ---------------------------------------------------

DL1 = 0.8001322
GAP1 = 0.1765959
GAP2 = 0.2960518
qs1_length = 0.2997797
qs2_length = 0.2585548
qs3_length = 0.2382791
DL2 = 2.1162209
GAP3 = 0.1978111


traj = (
    Trajectory.set_start_point(P2.origin()).first_line(P2.x_direct(), DL1)
    .add_arc_line(radius=0.95, clockwise=True, angle_deg=45/2).as_aperture_objrct_on_last(20*MM)
    .add_strait_line(length=GAP1)
    .add_strait_line(length=qs1_length).as_aperture_objrct_on_last(20*MM)
    .add_strait_line(length=GAP2)
    .add_strait_line(length=qs2_length).as_aperture_objrct_on_last(20*MM)
    .add_strait_line(length=GAP2)
    .add_strait_line(length=qs1_length).as_aperture_objrct_on_last(20*MM)
    .add_strait_line(length=GAP1)
    .add_arc_line(radius=0.95, clockwise=True, angle_deg=45/2).as_aperture_objrct_on_last(20*MM)

    .add_strait_line(length=DL1)
    .add_strait_line(length=DL2)
    .add_arc_line(radius=0.95, clockwise=False, angle_deg=135/2).as_aperture_objrct_on_last(80*MM)
    .add_strait_line(length=GAP3)
    .add_strait_line(length=qs3_length).as_aperture_objrct_on_last(60*MM)
    .add_strait_line(length=GAP3)
    .add_arc_line(radius=0.95, clockwise=False, angle_deg=135/2).as_aperture_objrct_on_last(80*MM)
    .add_strait_line(length=DL2)
)


print(f"end = {traj.point_at_end()}")

Plot2.plot(traj,describe=['r-','b-','k-'])
Plot2.plot(StraightLine2(length=traj.point_at_end().x, direct=P2.x_direct(),
                         start_point=P2.origin()), describe='k--')

Plot2.equal()


# ---------------------------------------------------

DL1 = 0.8001322
GAP1 = 0.1765959  * 0 +0.4
GAP2 = 0.2960518
qs1_length = 0.2997797
qs2_length = 0.2585548
qs3_length = 0.2382791
DL2 = 1.9
GAP3 = 0.1978111 * 0 +0.5


traj = (
    Trajectory.set_start_point(P2.origin()).first_line(P2.x_direct(), DL1)
    .add_arc_line(radius=0.95, clockwise=False, angle_deg=45/2).as_aperture_objrct_on_last(20*MM)
    .add_strait_line(length=GAP1)
    .add_strait_line(length=qs1_length).as_aperture_objrct_on_last(20*MM)
    .add_strait_line(length=GAP2)
    .add_strait_line(length=qs2_length).as_aperture_objrct_on_last(20*MM)
    .add_strait_line(length=GAP2)
    .add_strait_line(length=qs1_length).as_aperture_objrct_on_last(20*MM)
    .add_strait_line(length=GAP1)
    .add_arc_line(radius=0.95, clockwise=False, angle_deg=45/2).as_aperture_objrct_on_last(20*MM)

    .add_strait_line(length=DL1)
    .add_strait_line(length=DL2)
    .add_arc_line(radius=0.95, clockwise=True, angle_deg=135/2).as_aperture_objrct_on_last(80*MM)
    .add_strait_line(length=GAP3)
    .add_strait_line(length=qs3_length).as_aperture_objrct_on_last(60*MM)
    .add_strait_line(length=GAP3)
    .add_arc_line(radius=0.95, clockwise=True, angle_deg=135/2).as_aperture_objrct_on_last(80*MM)
    .add_strait_line(length=DL2)
)


print(f"end = {traj.point_at_end()}")

Plot2.plot(traj,describe=['r-','b-','k-'])
Plot2.plot(StraightLine2(length=traj.point_at_end().x, direct=P2.x_direct(),
                         start_point=P2.origin()), describe='k--')

Plot2.equal()
Plot2.show()
