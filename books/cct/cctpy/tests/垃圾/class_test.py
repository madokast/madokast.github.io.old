# 绘制 x 和 y 方向相椭圆

try:
    from books.cct.cctpy.cctpy import *
except ModuleNotFoundError:
    pass

import sys
sys.path.append(r'C:\Users\madoka_9900\Documents\github\madokast.github.io\books\cct\cctpy')

from cctpy import Beamline,P2,MM,BaseUtils
import matplotlib.pyplot as plt

bl:Beamline = ( # QS 磁铁加前后 1m 漂移段
    Beamline.set_start_point(P2.origin())
      .first_drift(direct=P2.x_direct(),length=1.0)
      .append_qs(
          length=0.27,gradient=0,
          second_gradient=-1000,aperture_radius=60*MM
      ).append_drift(1.0)
)

x,y = bl.track_phase_ellipse(
    x_sigma_mm=3.5,xp_sigma_mrad=7.5,
    y_sigma_mm=3.5,yp_sigma_mrad=7.5,
    delta=0.0,kinetic_MeV=250,
    particle_number=32,footstep=1*MM,
    concurrency_level = 16
)

plt.subplot(121)
plt.plot(*P2.extract(x),'r.')
plt.xlabel(xlabel='x/mm')
plt.ylabel(ylabel='xp/mr')
plt.title(label='x-plane')
plt.axis("equal")

plt.subplot(122)
plt.plot(*P2.extract(y),'r.')
plt.xlabel(xlabel='y/mm')
plt.ylabel(ylabel='yp/mr')
plt.title(label='y-plane')
plt.axis("equal")

plt.show()