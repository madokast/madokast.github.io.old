# 导入模块
import matplotlib.pyplot as plt
import numpy



raw_data = numpy.loadtxt("dp0.txt")

data = []
for i in range(len(raw_data)):
    x = raw_data[i][0]
    y = raw_data[i][1]
    if (abs(x)<15 and abs(y)<15) or True:
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
font_size = 25
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

plt.show()