from cctpy import *

# LOGO
width = 40
c_angle = 300
r = 1 * MM
c_r = 100
C1 = (
    Trajectory.set_start_point(start_point=P2(176,88))
    .first_line2(direct=P2.x_direct().rotate(BaseUtils.angle_to_radian(360-c_angle)/2), length=width)
    .add_arc_line(radius=r, clockwise=False, angle_deg=90)
    .add_arc_line(radius=c_r, clockwise=False, angle_deg=c_angle)
    .add_arc_line(radius=r, clockwise=False, angle_deg=90)
    .add_strait_line(length=width)
    .add_arc_line(radius=r, clockwise=False, angle_deg=90)
    .add_arc_line(radius=c_r-width, clockwise=True, angle_deg=c_angle)
)
C2 = C1+P2(200,0)

t_width = 190
t_height = 190

T = (
    Trajectory.set_start_point(start_point=P2(430,155))
    .first_line2(direct=P2.x_direct(), length=t_width)
    .add_arc_line(radius=r, clockwise=True, angle_deg=90)
    .add_strait_line(length=width)
    .add_arc_line(radius=r, clockwise=True, angle_deg=90)
    .add_strait_line(length=(t_width/2-width/2))
    .add_arc_line(radius=r, clockwise=False, angle_deg=90)
    .add_strait_line(length=t_height-width)
    .add_arc_line(radius=r, clockwise=True, angle_deg=90)
    .add_strait_line(length=width)
    .add_arc_line(radius=r, clockwise=True, angle_deg=90)
    .add_strait_line(length=t_height-width)
    .add_arc_line(radius=r, clockwise=False, angle_deg=90)
    .add_strait_line(length=(t_width/2-width/2))
    .add_arc_line(radius=r, clockwise=True, angle_deg=90)
    .add_strait_line(length=width)
) + P2(0,-5)

p_height = t_height
p_r = 50
width = 45

P_out = (
    Trajectory.set_start_point(start_point=P2(655,155))
    .first_line2(direct=P2.x_direct(), length=2*width)
    .add_arc_line(radius=p_r, clockwise=True, angle_deg=180)
    .add_strait_line(length=width)
    .add_arc_line(radius=r, clockwise=False, angle_deg=90)
    .add_strait_line(length=p_height-p_r*2)
    .add_arc_line(radius=r, clockwise=True, angle_deg=90)
    .add_strait_line(length=width)
    .add_arc_line(radius=r, clockwise=True, angle_deg=90)
    .add_strait_line(length=p_height)
) + P2(0,-5)

P_in = (
    Trajectory.set_start_point(start_point=P_out.point_at(width)-P2(0,width*0.6))
    .first_line2(direct=P2.x_direct(), length=width)
    .add_arc_line(radius=p_r-width*0.6, clockwise=True, angle_deg=180)
    .add_strait_line(length=width)
    .add_arc_line(radius=r, clockwise=True, angle_deg=90)
    .add_strait_line(length=(p_r-width*0.6)*2)
)

width = 40
y_width=50
y_heigt = t_height
y_tilt_len = 120

Y = (
    Trajectory.set_start_point(start_point=P2(810,155))
    .first_line2(direct=P2.x_direct(), length=width)
    .add_arc_line(radius=r, clockwise=True, angle_deg=60)
    .add_strait_line(length=y_tilt_len)
    .add_arc_line(radius=r, clockwise=False, angle_deg=120)
    .add_strait_line(length=y_tilt_len)
    .add_arc_line(radius=r, clockwise=True, angle_deg=60)
    .add_strait_line(length=width)
    .add_arc_line(radius=r, clockwise=True, angle_deg=120)
    .add_strait_line(length=y_tilt_len*1.3)
    .add_arc_line(radius=r, clockwise=False, angle_deg=30)
    .add_strait_line(length=t_height*0.4)
    .add_arc_line(radius=r, clockwise=True, angle_deg=90)
    .add_strait_line(length=width*1.1)
    .add_arc_line(radius=r, clockwise=True, angle_deg=90)
    .add_strait_line(length=t_height*0.4)
    .add_arc_line(radius=r, clockwise=False, angle_deg=30)
    .add_strait_line(length=y_tilt_len*1.3)
)

Plot3.plot_line2(C1,step=1)
Plot3.plot_line2(C2,step=1)
Plot3.plot_line2(T,step=1)
Plot3.plot_line2(P_out,step=1)
Plot3.plot_line2(P_in,step=1)
Plot3.plot_line2(Y,step=1)
Plot3.remove_background_color()
Plot3.PLT.axis('off')
Plot3.show()