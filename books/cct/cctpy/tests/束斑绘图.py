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

raw_data = numpy.loadtxt("dp-7.txt")

data = []
for i in range(len(raw_data)):
    x = raw_data[i][0]
    y = raw_data[i][1]
    if (abs(x)<12.5 and abs(y)<12.5):
        data.append([x,y])

data = numpy.array(data)

print(len(raw_data))
print(len(data))

# 数据
x = data[:,0]
y = data[:,1]

sigma_x = numpy.std(x)
sigma_y = numpy.std(y)

print(sigma_x,sigma_y)

# 画图
plt.hist2d(x=x, y=y, bins=100, cmap=plt.cm.Spectral_r,)



# 展示
lim = 15
plt.xlim(-lim,lim)
plt.ylim(-lim,lim)

font_family = "Times New Roman"
font_size = 40
font_label = {
    "family": font_family,
    "weight": "normal",
    "size": font_size,
}
plt.xlabel(xlabel="x/mm", fontdict=font_label)
plt.ylabel(ylabel="y/mm", fontdict=font_label)
plt.title(label="", fontdict=font_label)

plt.xticks(fontproperties=font_family, size=font_size)
plt.yticks(fontproperties=font_family, size=font_size)

if False: # 椭圆拟合
    miu1 = numpy.average(data[:,0])
    miu2 = numpy.average(data[:,1])

    miu = P2(miu1,miu2)
    print("均值向量",miu)

    sigma = numpy.cov(data.T)
    print("样本协方差",sigma)

    u,v=numpy.linalg.eig(sigma)
    lambd1 = u[0]
    lambd2 = u[1]
    a1 = P2.from_numpy_ndarry(v[:,0])
    a2 = P2.from_numpy_ndarry(v[:,1])
    print("lambd",u)
    print("array",v)

    print("v[0]",v[:,0])
    print("v[0].T",numpy.array([v[:,0]]).T)
    print(numpy.mat(sigma)*numpy.mat(numpy.array([v[:,0]]).T))
    print(lambd1*numpy.array([v[:,0]]).T)

    p2s = P2.from_numpy_ndarry(data)

    # Plot2.plot(p2s,describe='r.')
    # plt.scatter(data[:,0],data[:,1],s=3,c='r')

    # Plot2.plot([miu,miu+a1*numpy.sqrt(lambd1)*3],describe='r-')
    # Plot2.plot([miu,miu-a1*numpy.sqrt(lambd1)*3],describe='r-')
    # Plot2.plot([miu,miu+a2*numpy.sqrt(lambd2)*3],describe='r-')
    # Plot2.plot([miu,miu-a2*numpy.sqrt(lambd2)*3],describe='r-')

    Plot2.plot([miu,miu+a1*numpy.sqrt(lambd1)*2],describe='k-')
    Plot2.plot([miu,miu-a1*numpy.sqrt(lambd1)*2],describe='k-')
    Plot2.plot([miu,miu+a2*numpy.sqrt(lambd2)*2],describe='k-')
    Plot2.plot([miu,miu-a2*numpy.sqrt(lambd2)*2],describe='k-')

if True:
    miu1 = numpy.average(data[:,0])
    miu2 = numpy.average(data[:,1])

    miu = P2(miu1,miu2)
    print("均值向量",miu)

    # 2倍标准差
    sx = 2*numpy.std(data[:,0])
    sy = 2*numpy.std(data[:,1])
    print("2倍标准差",sx,sy)

    # 椭圆
    e = BaseUtils.Ellipse.create_standard_ellipse(sx,sy)
    es = [ep + miu for ep in e.uniform_distribution_points_along_edge(64)]

    Plot2.plot(es)
    Plot2.equal()


plt.show()