try:
    from books.cct.cctpy.cctpy import *
except ModuleNotFoundError:
    pass

from cctpy import *


Plot2.plot_p2(P2.x_direct(),'r.')
Plot2.info("xx",'yy','tt',12)
Plot2.show()