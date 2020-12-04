try:
    from books.cct.cctpy.cctpy import *
except ModuleNotFoundError:
    pass

from cctpy import *

bending_angle = 45
winding_number = 10

t = (
    Trajectory.set_start_point(start_point=P2.origin())
    .first_line2(direct=P2.x_direct(),length= 1*M)
    .add_arc_line(radius=0.95*M,clockwise=False,angle_deg=bending_angle)
    .add_strait_line(length=1*M)
)



cct =  CCT.create_cct_along(
    trajectory=t,
    s=1*M,
    big_r=0.95*M,
    small_r=80*MM,
    bending_angle=45,
    tilt_angles=[30],
    winding_number=10,
    current=1000,
    starting_point_in_ksi_phi_coordinate=P2.origin(),
    end_point_in_ksi_phi_coordinate=P2(winding_number*2*math.pi,BaseUtils.angle_to_radian(bending_angle)),
    disperse_number_per_winding=120
)

Plot3.plot_line2(t,describe='y--')
Plot3.plot_cct(cct,describe='r')
Plot3.plot_local_coordinate_system(cct.local_coordinate_system,axis_lengths=[2,2,1],describe='k-')
Plot3.plot_local_coordinate_system(LocalCoordinateSystem.global_coordinate_system(),axis_lengths=[0.5]*3,describe='b-')

Plot3.set_center(t.point_at(t.get_length()/2).to_p3(),cube_size=t.get_length()/2)
Plot3.off_axis()
Plot3.remove_background_color()
Plot3.show()

# bz =  cct.magnetic_field_bz_along(t)
# Plot2.plot_p2s(bz)
# Plot2.show()