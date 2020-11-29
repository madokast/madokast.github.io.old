import numpy as np

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

# Plot2.plot_cct(cct)
# Plot3.plot_cct(cct)

# Plot2.show()
# Plot3.show()


m = cct.magnetic_field_at(P3())

print(m)

import time
s = time.time()
for p in BaseUtils.linspace(P3(),P3(y=2),500):
    print(cct.magnetic_field_at(p))
print(time.time()-s)