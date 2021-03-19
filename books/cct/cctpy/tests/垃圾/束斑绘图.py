import os
import sys
curPath = os.path.abspath(os.path.dirname(__file__))
rootPath = os.path.split(curPath)[0]
PathProject = os.path.split(rootPath)[0]
sys.path.append(rootPath)
sys.path.append(PathProject)

from cctpy import *

try:
    from books.cct.cctpy.cctpy import *
except Exception as e:
    pass

data = numpy.loadtxt("dp5.txt")

miu1 = numpy.average(data[:,0])
miu2 = numpy.average(data[:,1])

miu = P2(miu1,miu2)
print("均值向量",miu)

sigma = numpy.cov(data.T)
print("样本协方差",sigma)

u,v=numpy.linalg.eig(sigma)
lambd1 = u[0]
lambd2 = u[1]
a1 = P2.from_numpy_ndarry(v[0])
a2 = P2.from_numpy_ndarry(v[1])
print(u)
print(v)

p2s = P2.from_numpy_ndarry(data)

# Plot2.plot(p2s,describe='r.')
plt.scatter(data[:,0],data[:,1],s=1,c='r')

lim = 15
Plot2.xlim(-lim,lim)
Plot2.ylim(-lim,lim)
Plot2.info("x/mm","y/mm","",25)

Plot2.plot([miu,miu+a1*numpy.sqrt(lambd1)*2],describe='k-')
Plot2.plot([miu,miu-a1*numpy.sqrt(lambd1)*2],describe='k-')
Plot2.plot([miu,miu+a2*numpy.sqrt(lambd2)*2],describe='k-')
Plot2.plot([miu,miu-a2*numpy.sqrt(lambd2)*2],describe='k-')
Plot2.show()