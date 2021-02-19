"""
使用方法
1. 使用 cctpy 中的 ParticleFactory.distributed_particles() 生成随机的相空间粒子
2. 使用 cosy_utils 的 SRConvertor.to_cosy_sr() 将相空间粒子输出为 COSY 的 SR 格式

两个方法的详细说明，见代码注释

以下给出3个使用示例
"""


from cctpy import *
from cosy_utils import *

# 束流参数为 x=y=3.5mm，xp,yp=7.5mr，dp=8%。生成粒子数目20

# 1. 生成x/xp相椭圆圆周上，动量分散为0的粒子，同时 y=3.5mm，yp=7.5mr。在 COSY 绘图为蓝色
ps = ParticleFactory.distributed_particles(
    3.5*MM, 7.5*MRAD, 3.5*MM, 7.5*MRAD, 0.0, number=20,
    distribution_area=ParticleFactory.DISTRIBUTION_AREA_EDGE,
    x_distributed=True, xp_distributed=True
)
cosy_code = SRConvertor.to_cosy_sr(ps,color=SRConvertor.COLOR_BLUE)
print('-----------生成x/xp相椭圆圆周上，动量分散为0的粒子，同时 y=3.5mm，yp=7.5mr。在 COSY 绘图为蓝色---------')
print(cosy_code)
print('\n\n')



# 2. 生成y/yp相椭圆内部，动量分散均为0.05的粒子。在 COSY 绘图为红色
ps = ParticleFactory.distributed_particles(
    0, 0, 3.5*MM, 7.5*MRAD, 0.05, number=20,
    distribution_area=ParticleFactory.DISTRIBUTION_AREA_FULL,
    y_distributed=True, yp_distributed=True
)
cosy_code = SRConvertor.to_cosy_sr(ps,color=SRConvertor.COLOR_RED)
print('-----------生成y/yp相椭圆内部，动量分散均为0.05的粒子，同时 x=xp=0。在 COSY 绘图为红色---------')
print(cosy_code)
print('\n\n')



# 3. 生成 x/xp/delta 三维相椭球球面的粒子。在 COSY 绘图为黑色
ps = ParticleFactory.distributed_particles(
    3.5*MM, 7.5*MRAD, 3.5*MM, 7.5*MRAD, 0.08, number=20,
    distribution_area=ParticleFactory.DISTRIBUTION_AREA_EDGE,
    x_distributed=True, xp_distributed=True, delta_distributed=True
)
cosy_code = SRConvertor.to_cosy_sr(ps,color=SRConvertor.COLOR_BLACK)
print('-----------生成 x/xp/delta 三维相椭球球面的粒子。在 COSY 绘图为蓝色---------')
print(cosy_code)
print('\n\n')

