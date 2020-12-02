try:
    from books.cct.cctpy.cctpy import *
except ModuleNotFoundError:
    pass

from cctpy import *

t = (
    Trajectory.set_start_point(start_point=P2.origin())
    .first_line2(direct=P2.x_direct(),length= 1*M)
    .add_arc_line(radius=0.95*M,clockwise=False,angle_deg=90)
    .add_strait_line(length=1*M)
)

CCT.create_cct_along(
    trajectory=t,
    s=1*M,
    big_r=0.95*M,
    small_r=80*MM
)