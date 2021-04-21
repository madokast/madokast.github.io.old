import os
import sys
curPath = os.path.abspath(os.path.dirname(__file__))
rootPath = os.path.split(curPath)[0]
PathProject = os.path.split(rootPath)[0]
sys.path.append(rootPath)
sys.path.append(PathProject)

from cctpy import *
from cosy_utils import *

try:
    from books.cct.cctpy.cctpy import *
    from books.cct.cctpy.cosy_utils import *
except Exception as e:
    pass

if __name__ == "__main__":

    BaseUtils.i_am_sure_my_code_closed_in_if_name_equal_main()

    data = [4.675,	41.126 	, 88.773,	98.139,
            91.748 	, 101.792,	62.677,	89.705,
            9409.261,	-7107.359, 25, 40, 34]  # *99.8/100

    gantry = HUST_SC_GANTRY(
        qs3_gradient=data[0],
        qs3_second_gradient=data[1],
        dicct345_tilt_angles=[30, data[2], data[3], data[4]],
        agcct345_tilt_angles=[data[5], 30, data[6], data[7]],
        dicct345_current=data[8],
        agcct345_current=data[9],
        agcct3_winding_number=data[10],
        agcct4_winding_number=data[11],
        agcct5_winding_number=data[12],
        agcct3_bending_angle=-67.5*(data[10])/(data[10]+data[11]+data[12]),
        agcct4_bending_angle=-67.5*(data[11])/(data[10]+data[11]+data[12]),
        agcct5_bending_angle=-67.5*(data[12])/(data[10]+data[11]+data[12]),

        DL1=0.9007765,
        GAP1=0.4301517,
        GAP2=0.370816,
        qs1_length=0.2340128,
        qs1_aperture_radius=60 * MM,
        qs1_gradient=0.0,
        qs1_second_gradient=0.0,
        qs2_length=0.200139,
        qs2_aperture_radius=60 * MM,
        qs2_gradient=0.0,
        qs2_second_gradient=0.0,

        DL2=2.35011,
        GAP3=0.43188,
        qs3_length=0.24379,

        agcct345_inner_small_r=92.5 * MM + 0.1*MM,  # 92.5
        agcct345_outer_small_r=108.5 * MM + 0.1*MM,  # 83+15
        dicct345_inner_small_r=124.5 * MM + 0.1*MM,  # 83+30+1
        dicct345_outer_small_r=140.5 * MM + 0.1*MM,  # 83+45 +2
    )
    bl_all = gantry.create_beamline()

    f = gantry.first_bending_part_length()

    sp = bl_all.trajectory.point_at(f)
    sd = bl_all.trajectory.direct_at(f)

    start = time.time()
    bl = gantry.create_second_bending_part(sp, sd)

    if True:
        ps = [
            PhaseSpaceParticle(0,0,0,0,0,delta=0),
            PhaseSpaceParticle(0,0,0,0,0,delta=-0.07),
            PhaseSpaceParticle(0,0,0,0,0,delta=0.05),
        ]

        ip_start = ParticleFactory.create_proton_along(
            bl.trajectory, 0.0, 215)

        rps = ParticleFactory.create_from_phase_space_particles(
            ip_start, ip_start.get_natural_coordinate_system(), ps)

        for p in rps:
            traj = ParticleRunner.run_get_trajectory(p,bl,bl.get_length())
            Plot3.plot_p3s(traj,describe='k-')

        

    ms = bl.magnets

    dicct0 = CCT.as_cct(ms[0])
    dicct1 = CCT.as_cct(ms[1])

    agcct0 = CCT.as_cct(ms[2])
    agcct1 = CCT.as_cct(ms[3])
    agcct2 = CCT.as_cct(ms[4])
    agcct3 = CCT.as_cct(ms[5])
    agcct4 = CCT.as_cct(ms[6])
    agcct5 = CCT.as_cct(ms[7])

    # Plot3.plot_cct(dicct0,describe='r-')
    Plot3.plot_cct(dicct1,describe='r-')

    # Plot3.plot_cct(agcct0,describe='b-')
    # Plot3.plot_cct(agcct1,describe='b-')

    # Plot3.plot_cct(agcct2,describe='y-')
    # Plot3.plot_cct(agcct3,describe='y-')

    # Plot3.plot_cct(agcct4,describe='b-')
    # Plot3.plot_cct(agcct5,describe='b-')

    Plot3.set_center(center=dicct0.global_path3()[len(dicct0.global_path3())//2],cube_size=1)

    Plot3.off_axis()
    Plot3.remove_background_color()

    Plot3.ax.view_init(elev=45., azim=-15)

    Plot3.show()