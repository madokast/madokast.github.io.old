"""
显示三维数据
"""

from tvtk.api import tvtk

# 创建长方体数据源
s = tvtk.CubeSource(x_length=1.0, y_length=2.0, z_length=3.0)

# 数据源 转为图形数据
m = tvtk.PolyDataMapper(input_connection=s.output_port)

# 创建 actor 指定 图形数据
a = tvtk.Actor(mapper=m)

# 创建 renderer，将 actor 传入
r = tvtk.Renderer(background=(0,0,0))
r.add_actor(a)

# 创建一个窗口，把 renderer 传入
w = tvtk.RenderWindow(size=(300,300))
w.add_renderer(r)

# 创建一个交互器，把 窗口 传入
i = tvtk.RenderWindowInteractor(render_window=w)

# 开启
i.initialize()
i.start()
