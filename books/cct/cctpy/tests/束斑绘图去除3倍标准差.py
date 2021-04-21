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

raw_data = numpy.loadtxt("dp5.txt")

# 计算均值标准差
avgx = numpy.average(raw_data[:, 0])
avgy = numpy.average(raw_data[:, 1])

miu = P2(avgx, avgy)
print("均值向量", miu)

# 标准差
stdx = numpy.std(raw_data[:, 0])
stdy = numpy.std(raw_data[:, 1])
print("标准差", stdx, stdy)

stdx2 = stdx*2
stdy2 = stdy*2

stdx3 = stdx*3
stdy3 = stdy*3

data = []
for i in range(len(raw_data)):
    x = raw_data[i][0]
    y = raw_data[i][1]
    if ((x-avgx)**2)/(stdx3**2) + ((y - avgy)**2)/(stdy3**2) < 1:
        data.append([x, y])

data = numpy.array(data)

print(len(raw_data))
print(len(data))

# 数据
x = data[:, 0]
y = data[:, 1]

# 画图
plt.hist2d(x=x, y=y, bins=100, cmap=plt.cm.Spectral_r,)


# 展示
lim = 15
plt.xlim(-lim, lim)
plt.ylim(-lim, lim)

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

# 椭圆
e = BaseUtils.Ellipse.create_standard_ellipse(stdx2, stdy2)
es = [ep + miu for ep in e.uniform_distribution_points_along_edge(64)]
es.append(es[0])

Plot2.plot(es)
Plot2.equal()


plt.show()
