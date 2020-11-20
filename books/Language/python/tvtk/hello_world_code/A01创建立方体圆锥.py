"""
画一个立方体
"""

from tvtk.api import tvtk

# 立方体
s = tvtk.CubeSource(x_length=1.0, y_length=2.0, z_length=3.0)

print(s)

"""
  X Length: 1
  Y Length: 2
  Z Length: 3  三个轴的长度
  Center: (0, 0, 0) 立方体中心
  Output Points Precision: 0 精度
"""

###################################################
# 圆锥
s = tvtk.ConeSource(height=3.0,radius=1.0,resolution=36)

print(s)

