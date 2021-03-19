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


print(Protons.convert_momentum_dispersion_to_energy_dispersion(0.1,250))


for i in range(-15,15,1):
    print(i/100,Protons.convert_momentum_dispersion_to_energy_dispersion(i/100,250))