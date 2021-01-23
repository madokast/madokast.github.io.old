import matplotlib.pyplot as plt
import numpy

# 线槽宽度深度 （无意义）
# width = 3.2
# depth = 11
# 超导线数目
# line_number = 2*7
# 超导线横截面积
# line_cross_section = math.pi * ((1.25/2)**2)
# 导线数据，电流 - 磁场
# product_data_current = [1100, 885, 640, 390] # r = 1.5
# product_data_magnet = [6, 7, 8, 9] # r = 1.5
product_data_current = [795, 620, 445, 275] # r = 1.25
product_data_magnet = [6, 7, 8, 9] # r = 1.25

# 总电流
current_dicct = 9409 #9536.310 #9426.734 #9488.615
current_agcct = 7107 #6259.974 #5626.101 #7334.914

# 最大磁场
max_magnet_dicct = 4.12 #4.05 #3.78 #4.102784
max_magnet_agcct = 4.52 #4.27 #3.82 #4.596925

# 电流密度，注意除总截面积而不是线槽面积
# current_density_dicct = current_dicct / line_cross_section / line_number
# current_density_agcct = current_agcct / line_cross_section / line_number

# print(f"current_density_dicct={current_density_dicct}")
# print(f"current_density_agcct={current_density_agcct}")


def work_line(current_density, max_manget, end_current_density):
    # 工作点
    plt.plot([max_manget], [current_density], 'r.', markersize=10)
    # 工作线
    plt.plot([0, max_manget], [0, current_density], 'r-')
    # 延长线
    plt.plot([0, end_current_density*max_manget/current_density], [0, end_current_density], 'r--')


# work_line(current_agcct / (2*8), max_magnet_agcct, 1000) # 640
work_line(current_agcct / (2*8), max_magnet_agcct, 656) # 640
# work_line(current_agcct / (2*6), max_magnet_agcct, 1000) # 640
# work_line(current_agcct / (2*5), max_magnet_agcct, 1000) # 640

# work_line(current_dicct / (2*8), max_magnet_dicct, 1000) # 802
work_line(current_dicct / (2*8), max_magnet_dicct, 823) # 802
# work_line(current_dicct / (2*6), max_magnet_dicct, 1000) # 802
# work_line(current_dicct / (2*5), max_magnet_dicct, 1000) # 802

def product_line(product_data_magnet, product_data_current):
    # plt.plot(product_data_magnet, product_data_current,"k--")
    plt.plot(product_data_magnet, product_data_current,"k.", markersize=10)
    plt.xlim(0, 10)
    plt.ylim(0, 1200)
    plt.xlabel('B/T',fontsize=20)
    plt.ylabel('I/A',fontsize=20)

    # 拟合
    fit = numpy.polyfit(product_data_magnet, product_data_current,2)
    f = numpy.poly1d(fit)
    xs = numpy.linspace(0, 10,1000)
    ys = f(xs)
    plt.plot(xs,ys,'k--')

product_line(product_data_magnet, product_data_current)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.show()
