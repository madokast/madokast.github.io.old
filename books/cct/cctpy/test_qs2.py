try:
    from books.cct.cctpy.cctpy import *
except ModuleNotFoundError:
    pass

from cctpy import *

trajectory = Trajectory.set_start_point(P2.origin()).first_line2(
    direct=P2.y_direct(), length=1 * M
)

qs = QS.create_qs_along(
    trajectory,
    s=0.4 * M,
    length=0.27 * M,
    gradient=20,
    second_gradient=0,
    aperture_radius=60 * MM,
)

qs.local_coordinate_system.location -= qs.local_coordinate_system.XI * 30 * MM

Plot3.plot_line2(trajectory, describe="y--")
Plot3.plot_qs(qs, describe="b-")
Plot3.plot_local_coordinate_system(
    qs.local_coordinate_system, axis_lengths=[0.1, 0.1, 0.3], describe="k-"
)
Plot3.set_box(P3(-0.2, 0, -0.2), P3(0.2, 1, 0.2))
Plot3.off_axis()
Plot3.remove_background_color()
Plot3.show()