try:
    from books.cct.cctpy.cctpy import *
except ModuleNotFoundError:
    pass

from cctpy import *

nm = Magnet.no_magnet()

print(nm.magnetic_field_at(P3.zeros()))