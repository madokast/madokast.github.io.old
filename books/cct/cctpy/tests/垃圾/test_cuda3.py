# -*- coding: utf-8 -*-

from cctpy import (
    BaseUtils,
    P2,
    P3,
    StraightLine2,
    Trajectory,
    Plot2,
    Plot3,
    CCT,
    LocalCoordinateSystem,
    MM,
)

import time
import numpy as np

cct = CCT(
    LocalCoordinateSystem.global_coordinate_system(),
    0.95,
    83 * MM + 15 * MM * 2,
    67.5,
    [30.0, 80.0, 90.0, 90.0],
    128,
    -9664,
    P2(0, 0),
    P2(128 * np.pi * 2, 67.5 / 180.0 * np.pi),
)

times = 2
s = time.time()
for x in np.linspace(0, 0.001, times):
    p = P3(x, 0, 0)
    m = cct.magnetic_field_at(p)
    print(m, p)
    m = cct.magnetic_field_at(p)
    print(m, p)
print(f"GPUINNER-d={time.time()-s}")