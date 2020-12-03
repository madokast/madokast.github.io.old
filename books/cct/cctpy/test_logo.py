try:
    from books.cct.cctpy.cctpy import Plot3
    from books.cct.cctpy.cctpy import LocalCoordinateSystem
    from books.cct.cctpy.cctpy import MM
except ModuleNotFoundError:
    pass

import matplotlib
matplotlib.use('TkAgg')
from cctpy import *

LOGO = Trajectory.__cctpy__()
Plot3.plot_line2s(LOGO, [1 * M], ["r-", "r-", "r-", "b-", "b-"])
Plot3.plot_local_coordinate_system(LocalCoordinateSystem(location=P3(z=-0.5e-6)),axis_lengths=[1000,200,1e-6],describe='k-')
Plot3.off_axis()
Plot3.remove_background_color()
Plot3.ax.view_init(elev=20, azim=-79)
Plot3.show()