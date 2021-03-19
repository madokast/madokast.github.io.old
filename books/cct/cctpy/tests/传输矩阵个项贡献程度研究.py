import os
import sys
curPath = os.path.abspath(os.path.dirname(__file__))
rootPath = os.path.split(curPath)[0]
PathProject = os.path.split(rootPath)[0]
sys.path.append(rootPath)
sys.path.append(PathProject)

from cctpy import *
from cosy_utils import *

try:
    from books.cct.cctpy.cctpy import *
    from books.cct.cctpy.cosy_utils import *
except Exception as e:
    pass


if __name__ == "__main__":
    map = CosyMap("""  0.1685013     -1.599878     0.0000000E+00 0.0000000E+00-0.2143961E-02 100000
  0.6473000    -0.2112797     0.0000000E+00 0.0000000E+00-0.1627490E-03 010000
  0.0000000E+00 0.0000000E+00 -1.112507    -0.2916618     0.0000000E+00 001000
  0.0000000E+00 0.0000000E+00 0.7461328    -0.7032598     0.0000000E+00 000100
  0.0000000E+00 0.0000000E+00 0.0000000E+00 0.0000000E+00  1.000000     000010
  0.1360362E-02-0.1925970E-03 0.0000000E+00 0.0000000E+00  2.585801     000001
  -10.45069     -30.54686     0.0000000E+00 0.0000000E+00 -78.51220     200000
  -26.17113     -12.41180     0.0000000E+00 0.0000000E+00 -44.83269     110000
  -4.110182     -1.222564     0.0000000E+00 0.0000000E+00 -10.63639     020000
  0.0000000E+00 0.0000000E+00  146.7463     -198.2287     0.0000000E+00 101000
  0.0000000E+00 0.0000000E+00  43.16272     -50.64796     0.0000000E+00 011000
   79.41930      27.32562     0.0000000E+00 0.0000000E+00 -38.80589     002000
  0.0000000E+00 0.0000000E+00 -11.55940      37.15282     0.0000000E+00 100100
  0.0000000E+00 0.0000000E+00 -1.315524      6.338662     0.0000000E+00 010100
  -27.68415     -2.450790     0.0000000E+00 0.0000000E+00  7.119962     001100
   94.03482      38.53409     0.0000000E+00 0.0000000E+00 -1.721695     100001
   25.42234      24.55839     0.0000000E+00 0.0000000E+00-0.2942295     010001
  0.0000000E+00 0.0000000E+00  50.05634     -56.73419     0.0000000E+00 001001
   6.020328     0.9735688     0.0000000E+00 0.0000000E+00 -33.65391     000200
  0.0000000E+00 0.0000000E+00  69.56705      24.64584     0.0000000E+00 000101
  0.5520378     0.7667639E-01 0.0000000E+00 0.0000000E+00 -2.264157     000002""")

    print(map)


    # p1_end = map.apply(p1,5,True)
    p1 = PhaseSpaceParticle(3.5*MM/2,7.5*MRAD/2,3.5*MM/2,7.5*MRAD/2,0,0.00)
    p2 = PhaseSpaceParticle(3.5*MM/2,7.5*MRAD/2,3.5*MM/2,7.5*MRAD/2,0,0.05)

    p1 = PhaseSpaceParticle.convert_delta_from_momentum_dispersion_to_energy_dispersion(p1,250)
    p2 = PhaseSpaceParticle.convert_delta_from_momentum_dispersion_to_energy_dispersion(p2,250)

    print(p1)
    print(p2)
    # print(p1_end)

    p2_end = map.apply(p2,5,True)
    print(p2_end)



