import os
import sys
curPath = os.path.abspath(os.path.dirname(__file__))
rootPath = os.path.split(curPath)[0]
PathProject = os.path.split(rootPath)[0]
sys.path.append(rootPath)
sys.path.append(PathProject)

from cctpy import *
from agcct_connector import *
from cctpy_ext import *
import numpy

try:
    from books.cct.cctpy.cctpy import *
    from books.cct.cctpy.cctpy_ext import *
    from books.cct.cctpy.agcct_connector import *
except:
    pass


if __name__ == "__main__":
    R = 0.95
    bl = (
        Beamline.set_start_point(start_point=P2(
            R, BaseUtils.angle_to_radian(-20)*R))
        .first_drift(P2.y_direct(), BaseUtils.angle_to_radian(20)*R)
        .append_agcct(
            big_r=R,
            small_rs=[140.5*MM + 0.1*MM, 124.5*MM+ 0.1*MM, 108.5*MM+ 0.1*MM, 92.5*MM+ 0.1*MM],
            bending_angles=[-17.05, -27.27, -23.18],  # [15.14, 29.02, 23.34]
            tilt_angles=[[30.0, 88.8, 98.1, 91.7],
                         [101.8, 30.0, 62.7, 89.7]],
            winding_numbers=[[128], [25, 40, 34]],
            currents=[9409.261, -7107.359],
            disperse_number_per_winding=36
        )
        .append_drift(BaseUtils.angle_to_radian(20)*R)
    )

    # depth, width = 11*MM, 3.2*MM
    depth, width = 10.8*MM, 2.8*MM

    ms = bl.magnets
    dicct_out = CCT.as_cct(ms[0])
    dicct_in = CCT.as_cct(ms[1])
    agcct3_in = CCT.as_cct(ms[2])
    agcct3_out = CCT.as_cct(ms[3])

    agcct4_in = CCT.as_cct(ms[4])
    agcct4_out = CCT.as_cct(ms[5])

    agcct5_in = CCT.as_cct(ms[6])
    agcct5_out = CCT.as_cct(ms[7])

    # connector_34_in = AGCCT_CONNECTOR(agcct3_in,agcct4_in)
    # connector_34_out = AGCCT_CONNECTOR(agcct3_out,agcct4_out)

    # connector_45_in = AGCCT_CONNECTOR(agcct4_in,agcct5_in)
    # connector_45_out = AGCCT_CONNECTOR(agcct4_out,agcct5_out)

    if True: # 二极CCT内层延展
        lcs = dicct_in.local_coordinate_system.copy()
        lcs.location = P3.origin()
        print(lcs)
        def p2_func(ksi): return dicct_in.p2_function(ksi)

        def p3_func(ksi): return lcs.point_to_global_coordinate(
            dicct_in.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        ksi0 = dicct_in.starting_point_in_ksi_phi_coordinate.x # 0.0
        ksi1 = dicct_in.end_point_in_ksi_phi_coordinate.x # 128*2*π

        ksi0_pre = -42*2*numpy.pi + ksi0
        ksi1_post = 42*2*numpy.pi + ksi1

        print(ksi0, ksi1)
        print(ksi0_pre, ksi1_post)

        ksi_list = BaseUtils.linspace(ksi0, ksi1, 8+1)

        ksi_list_pre = BaseUtils.linspace(ksi0_pre, ksi0, 3+1)
        ksi_list_post = BaseUtils.linspace(ksi1, ksi1_post, 3+1)

        print(ksi_list)
        print(ksi_list_pre)
        print(ksi_list_post)

        tangential_direct = BaseUtils.derivative(p3_func)

        def main_normal_direct(ksi): return lcs.point_to_global_coordinate(
            dicct_in.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        def second_normal_direction(ksi): return (
            tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        def path0(ksi): return p3_func(ksi)

        def path1(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path2(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path3(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        def path4(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        numpy.savetxt('dicct_in_pre_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0_pre,ksi0,42*360)]))
        numpy.savetxt('dicct_in_post_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi1,ksi1_post,42*360)]))
        for i in range(len(ksi_list_pre)-1):
            print(i)
            ksi0 = ksi_list_pre[i]
            ksi1 = ksi_list_pre[i+1]
            print(ksi0, ksi1)
            numpy.savetxt(f'layer3_dicct_in_pre_1_part_{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path1(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer3_dicct_in_pre_2_part_{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path2(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer3_dicct_in_pre_3_part_{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path3(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer3_dicct_in_pre_4_part_{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path4(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')

        for i in range(len(ksi_list_post)-1):
            print(i)
            ksi0 = ksi_list_post[i]
            ksi1 = ksi_list_post[i+1]
            print(ksi0, ksi1)
            numpy.savetxt(f'layer3_dicct_in_post_1_part_{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path1(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer3_dicct_in_post_2_part_{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path2(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer3_dicct_in_post_3_part_{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path3(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer3_dicct_in_post_4_part_{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path4(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
        # print('down')

    if True:  # 二极CCT外层延展
        lcs = dicct_out.local_coordinate_system.copy()
        lcs.location = P3.origin()
        print(lcs)
        def p2_func(ksi): return dicct_out.p2_function(ksi)

        def p3_func(ksi): return lcs.point_to_global_coordinate(
            dicct_out.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        ksi0 = dicct_out.starting_point_in_ksi_phi_coordinate.x # 0.0
        ksi1 = dicct_out.end_point_in_ksi_phi_coordinate.x # -128*2*π

        ksi0_pre = ksi0 + 42*2*numpy.pi
        ksi1_post = ksi1 - 42*2*numpy.pi


        print(ksi0, ksi1)
        print(ksi0_pre, ksi1_post)


        ksi_list = BaseUtils.linspace(ksi0, ksi1, 8+1)
        ksi_list_pre = BaseUtils.linspace(ksi0_pre, ksi0, 3+1)
        ksi_list_post = BaseUtils.linspace(ksi1, ksi1_post, 3+1)

        print(ksi_list)

        tangential_direct = BaseUtils.derivative(p3_func)

        def main_normal_direct(ksi): return lcs.point_to_global_coordinate(
            dicct_out.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        def second_normal_direction(ksi): return (
            tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        def path0(ksi): return p3_func(ksi)

        def path1(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path2(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path3(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        def path4(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        numpy.savetxt('dicct_out_pre_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0_pre,ksi0,42*360)]))
        numpy.savetxt('dicct_out_post_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi1,ksi1_post,42*360)]))

        for i in range(len(ksi_list_pre)-1):
            print(i)
            ksi0 = ksi_list_pre[i]
            ksi1 = ksi_list_pre[i+1]
            print(ksi0, ksi1)
            numpy.savetxt(f'layer4_dicct_out_pre_2_part{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path2(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer4_dicct_out_pre_3_part{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path3(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer4_dicct_out_pre_4_part{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path4(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer4_dicct_out_pre_1_part{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path1(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')

        for i in range(len(ksi_list_post)-1):
            print(i)
            ksi0 = ksi_list_post[i]
            ksi1 = ksi_list_post[i+1]
            print(ksi0, ksi1)
            numpy.savetxt(f'layer4_dicct_out_post_2_part{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path2(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer4_dicct_out_post_3_part{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path3(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer4_dicct_out_post_4_part{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path4(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer4_dicct_out_post_1_part{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path1(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')

    if True:  # 二极CCT内层
        lcs = dicct_in.local_coordinate_system.copy()
        lcs.location = P3.origin()
        print(lcs)
        def p2_func(ksi): return dicct_in.p2_function(ksi)

        def p3_func(ksi): return lcs.point_to_global_coordinate(
            dicct_in.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        ksi0 = dicct_in.starting_point_in_ksi_phi_coordinate.x
        ksi1 = dicct_in.end_point_in_ksi_phi_coordinate.x
        print(ksi0, ksi1)
        ksi_list = BaseUtils.linspace(ksi0, ksi1, 8+1)
        print(ksi_list)

        tangential_direct = BaseUtils.derivative(p3_func)

        def main_normal_direct(ksi): return lcs.point_to_global_coordinate(
            dicct_in.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        def second_normal_direction(ksi): return (
            tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        def path0(ksi): return p3_func(ksi)

        def path1(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path2(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path3(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        def path4(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        # numpy.savetxt('dicct_in_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,128*360)]))
        for i in range(len(ksi_list)-1):
            print(i)
            ksi0 = ksi_list[i]
            ksi1 = ksi_list[i+1]
            print(ksi0, ksi1)
            numpy.savetxt(f'layer3_dicct_in_1_part_{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path1(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 16*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer3_dicct_in_2_part_{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path2(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 16*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer3_dicct_in_3_part_{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path3(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 16*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer3_dicct_in_4_part_{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path4(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 16*360+1)]), fmt='%.4f')
        # print('down')

    if True:  # 二极CCT外层
        lcs = dicct_out.local_coordinate_system.copy()
        lcs.location = P3.origin()
        print(lcs)
        def p2_func(ksi): return dicct_out.p2_function(ksi)

        def p3_func(ksi): return lcs.point_to_global_coordinate(
            dicct_out.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        ksi0 = dicct_out.starting_point_in_ksi_phi_coordinate.x
        ksi1 = dicct_out.end_point_in_ksi_phi_coordinate.x
        print(ksi0, ksi1)
        ksi_list = BaseUtils.linspace(ksi0, ksi1, 8+1)
        print(ksi_list)

        tangential_direct = BaseUtils.derivative(p3_func)

        def main_normal_direct(ksi): return lcs.point_to_global_coordinate(
            dicct_out.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        def second_normal_direction(ksi): return (
            tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        def path0(ksi): return p3_func(ksi)

        def path1(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path2(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path3(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        def path4(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        numpy.savetxt('dicct_out_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,128*360)]))
        for i in range(len(ksi_list)-1):
            print(i)
            ksi0 = ksi_list[i]
            ksi1 = ksi_list[i+1]
            print(ksi0, ksi1)
            numpy.savetxt(f'layer4_dicct_out_2_part{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path2(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 16*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer4_dicct_out_3_part{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path3(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 16*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer4_dicct_out_4_part{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path4(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 16*360+1)]), fmt='%.4f')
            numpy.savetxt(f'layer4_dicct_out_1_part{i+1}.txt', 1000*numpy.array([P3.as_p3(
                path1(t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 16*360+1)]), fmt='%.4f')

    if True:  # 四极CCT内层 3 匝
        lcs = agcct3_in.local_coordinate_system.copy()
        lcs.location = P3.origin()
        print(lcs)
        def p2_func3(ksi): return agcct3_in.p2_function(ksi)
        # p2_func4 = lambda ksi:agcct4_in.p2_function(ksi)
        # p2_func5 = lambda ksi:agcct5_in.p2_function(ksi)

        ksi0_3 = agcct3_in.starting_point_in_ksi_phi_coordinate.x
        # ksi0_4 = agcct4_in.starting_point_in_ksi_phi_coordinate.x
        # ksi0_5 = agcct5_in.starting_point_in_ksi_phi_coordinate.x

        ksi1_3 = agcct3_in.end_point_in_ksi_phi_coordinate.x
        # ksi1_4 = agcct4_in.end_point_in_ksi_phi_coordinate.x
        # ksi1_5 = agcct5_in.end_point_in_ksi_phi_coordinate.x

        # connector34_p2_fun = connector_34_in.p2_function
        # connector45_p2_fun = connector_45_in.p2_function

        # 3
        fp = Function_Part(p2_func3, start=agcct3_in.starting_point_in_ksi_phi_coordinate.x,
                           end=agcct3_in.end_point_in_ksi_phi_coordinate.x)
        print(fp.length)
        print(ksi0_3, ksi1_3)
        ksi_list = numpy.array([0, 12, 25])*math.pi*2
        print(ksi_list)

        #############################

        def p2_func(ksi): return fp.valve_at(ksi)

        def p3_func(ksi): return lcs.point_to_global_coordinate(
            agcct3_in.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        # ksi0 = 0
        # ksi1 = fp.length

        tangential_direct = BaseUtils.derivative(p3_func)

        def main_normal_direct(ksi): return lcs.point_to_global_coordinate(
            agcct3_in.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        def second_normal_direction(ksi): return (
            tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        def path0(ksi): return p3_func(ksi)

        def path1(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path2(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path3(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        def path4(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        # numpy.savetxt('agcct_in_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,25*360)]))
        ksi0 = ksi_list[0]
        ksi1 = ksi_list[1]
        print(ksi0, ksi1)
        numpy.savetxt('layer1_agcct1_in_1_part1.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 12*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct1_in_2_part1.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 12*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct1_in_3_part1.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 12*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct1_in_4_part1.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 12*360+1)]), fmt='%.4f')

        ksi0 = ksi_list[1]
        ksi1 = ksi_list[2]
        print(ksi0, ksi1)
        numpy.savetxt('layer1_agcct1_in_1_part2.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct1_in_2_part2.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct1_in_3_part2.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct1_in_4_part2.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        # print('down')

    if True:  # 四极CCT内层 4 匝
        lcs = agcct3_in.local_coordinate_system.copy()
        lcs.location = P3.origin()
        print(lcs)
        # p2_func3 = lambda ksi:agcct3_in.p2_function(ksi)
        def p2_func4(ksi): return agcct4_in.p2_function(ksi)
        # p2_func5 = lambda ksi:agcct5_in.p2_function(ksi)

        # ksi0_3 = agcct3_in.starting_point_in_ksi_phi_coordinate.x
        ksi0_4 = agcct4_in.starting_point_in_ksi_phi_coordinate.x
        # ksi0_5 = agcct5_in.starting_point_in_ksi_phi_coordinate.x

        # ksi1_3 = agcct3_in.end_point_in_ksi_phi_coordinate.x
        ksi1_4 = agcct4_in.end_point_in_ksi_phi_coordinate.x
        # ksi1_5 = agcct5_in.end_point_in_ksi_phi_coordinate.x

        # connector34_p2_fun = connector_34_in.p2_function
        # connector45_p2_fun = connector_45_in.p2_function

        # 3
        fp = Function_Part(p2_func4, start=agcct4_in.starting_point_in_ksi_phi_coordinate.x,
                           end=agcct4_in.end_point_in_ksi_phi_coordinate.x)
        print(fp.length)
        print(ksi0_4, ksi1_4)
        ksi_list = numpy.array([0, 13, 26, 40])*math.pi*2
        print(ksi_list)

        #############################

        def p2_func(ksi): return fp.valve_at(ksi)

        def p3_func(ksi): return lcs.point_to_global_coordinate(
            agcct3_in.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        # ksi0 = 0
        # ksi1 = fp.length

        tangential_direct = BaseUtils.derivative(p3_func)

        def main_normal_direct(ksi): return lcs.point_to_global_coordinate(
            agcct3_in.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        def second_normal_direction(ksi): return (
            tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        def path0(ksi): return p3_func(ksi)

        def path1(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path2(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path3(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        def path4(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        # numpy.savetxt('agcct_in_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,25*360)]))
        ksi0 = ksi_list[0]
        ksi1 = ksi_list[1]
        print(ksi0, ksi1)
        numpy.savetxt('layer1_agcct2_in_1_part1.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct2_in_2_part1.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct2_in_3_part1.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct2_in_4_part1.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')

        ksi0 = ksi_list[1]
        ksi1 = ksi_list[2]
        print(ksi0, ksi1)
        numpy.savetxt('layer1_agcct2_in_1_part2.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct2_in_2_part2.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct2_in_3_part2.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct2_in_4_part2.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')

        ksi0 = ksi_list[2]
        ksi1 = ksi_list[3]
        print(ksi0, ksi1)
        numpy.savetxt('layer1_agcct2_in_1_part3.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct2_in_2_part3.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct2_in_3_part3.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct2_in_4_part3.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')

    if True:  # 四极CCT内层 5 匝
        lcs = agcct3_in.local_coordinate_system.copy()
        lcs.location = P3.origin()
        print(lcs)
        # p2_func3 = lambda ksi:agcct3_in.p2_function(ksi)
        # p2_func4 = lambda ksi:agcct4_in.p2_function(ksi)
        def p2_func5(ksi): return agcct5_in.p2_function(ksi)

        # ksi0_3 = agcct3_in.starting_point_in_ksi_phi_coordinate.x
        # ksi0_4 = agcct4_in.starting_point_in_ksi_phi_coordinate.x
        ksi0_5 = agcct5_in.starting_point_in_ksi_phi_coordinate.x

        # ksi1_3 = agcct3_in.end_point_in_ksi_phi_coordinate.x
        # ksi1_4 = agcct4_in.end_point_in_ksi_phi_coordinate.x
        ksi1_5 = agcct5_in.end_point_in_ksi_phi_coordinate.x

        # connector34_p2_fun = connector_34_in.p2_function
        # connector45_p2_fun = connector_45_in.p2_function

        # 3
        fp = Function_Part(p2_func5, start=agcct5_in.starting_point_in_ksi_phi_coordinate.x,
                           end=agcct5_in.end_point_in_ksi_phi_coordinate.x)
        print(fp.length)
        print(ksi0_5, ksi1_5)
        ksi_list = numpy.array([0, 17, 34])*math.pi*2
        print(ksi_list)

        #############################

        def p2_func(ksi): return fp.valve_at(ksi)

        def p3_func(ksi): return lcs.point_to_global_coordinate(
            agcct3_in.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        # ksi0 = 0
        # ksi1 = fp.length

        tangential_direct = BaseUtils.derivative(p3_func)

        def main_normal_direct(ksi): return lcs.point_to_global_coordinate(
            agcct3_in.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        def second_normal_direction(ksi): return (
            tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        def path0(ksi): return p3_func(ksi)

        def path1(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path2(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path3(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        def path4(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        # numpy.savetxt('agcct_in_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,25*360)]))
        ksi0 = ksi_list[0]
        ksi1 = ksi_list[1]
        print(ksi0, ksi1)
        numpy.savetxt('layer1_agcct3_in_1_part1.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct3_in_2_part1.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct3_in_3_part1.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct3_in_4_part1.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')

        ksi0 = ksi_list[1]
        ksi1 = ksi_list[2]
        print(ksi0, ksi1)
        numpy.savetxt('layer1_agcct3_in_1_part2.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct3_in_2_part2.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct3_in_3_part2.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        numpy.savetxt('layer1_agcct3_in_4_part2.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')

    if True:  # 四极CCT外层 3
        lcs = agcct3_out.local_coordinate_system.copy()
        lcs.location = P3.origin()
        print(lcs)
        def p2_func3(ksi): return agcct3_out.p2_function(ksi)
        # p2_func4 = lambda ksi:agcct4_out.p2_function(ksi)
        # p2_func5 = lambda ksi:agcct5_out.p2_function(ksi)

        ksi0_3 = agcct3_out.starting_point_in_ksi_phi_coordinate.x
        # ksi0_4 = agcct4_out.starting_point_in_ksi_phi_coordinate.x
        # ksi0_5 = agcct5_out.starting_point_in_ksi_phi_coordinate.x

        ksi1_3 = agcct3_out.end_point_in_ksi_phi_coordinate.x
        # ksi1_4 = agcct4_out.end_point_in_ksi_phi_coordinate.x
        # ksi1_5 = agcct5_out.end_point_in_ksi_phi_coordinate.x

        # connector34_p2_fun = connector_34_out.p2_function
        # connector45_p2_fun = connector_45_out.p2_function

        # 3
        fp = Function_Part(p2_func3, start=agcct3_out.starting_point_in_ksi_phi_coordinate.x,
                           end=agcct3_out.end_point_in_ksi_phi_coordinate.x)
        # fp = Function_Part(p2_func4,start=agcct4_out.starting_point_in_ksi_phi_coordinate.x,end=agcct4_out.end_point_in_ksi_phi_coordinate.x)
        # fp = Function_Part(p2_func5,start=agcct5_out.starting_point_in_ksi_phi_coordinate.x,end=agcct5_out.end_point_in_ksi_phi_coordinate.x)

        print(fp.length)
        print(ksi0_3, ksi1_3)
        ksi_list = numpy.array([0, 12, 25])*math.pi*2
        print(ksi_list)

        #############################

        def p2_func(ksi): return fp.valve_at(ksi)

        def p3_func(ksi): return lcs.point_to_global_coordinate(
            agcct3_out.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        ksi0 = 0
        ksi1 = fp.length

        tangential_direct = BaseUtils.derivative(p3_func)

        def main_normal_direct(ksi): return lcs.point_to_global_coordinate(
            agcct3_out.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        def second_normal_direction(ksi): return (
            tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        def path0(ksi): return p3_func(ksi)

        def path1(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path2(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path3(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        def path4(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        # numpy.savetxt('agcct_out_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,128*360)]))
        ksi0 = ksi_list[0]
        ksi1 = ksi_list[1]
        print(ksi0, ksi1)  # layer1_agcct1_in_1_part1
        numpy.savetxt('layer2_agcct1_out_1_part1.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 12*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct1_out_2_part1.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 12*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct1_out_3_part1.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 12*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct1_out_4_part1.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 12*360+1)]), fmt='%.4f')
        print('down')

        ksi0 = ksi_list[1]
        ksi1 = ksi_list[2]
        print(ksi0, ksi1)
        numpy.savetxt('layer2_agcct1_out_1_part2.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct1_out_2_part2.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct1_out_3_part2.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct1_out_4_part2.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        print('down')

    if True:  # 四极CCT外层 4
        lcs = agcct3_out.local_coordinate_system.copy()
        lcs.location = P3.origin()
        print(lcs)
        # p2_func3 = lambda ksi:agcct3_out.p2_function(ksi)
        def p2_func4(ksi): return agcct4_out.p2_function(ksi)
        # p2_func5 = lambda ksi:agcct5_out.p2_function(ksi)

        # ksi0_3 = agcct3_out.starting_point_in_ksi_phi_coordinate.x
        ksi0_4 = agcct4_out.starting_point_in_ksi_phi_coordinate.x
        # ksi0_5 = agcct5_out.starting_point_in_ksi_phi_coordinate.x

        # ksi1_3 = agcct3_out.end_point_in_ksi_phi_coordinate.x
        ksi1_4 = agcct4_out.end_point_in_ksi_phi_coordinate.x
        # ksi1_5 = agcct5_out.end_point_in_ksi_phi_coordinate.x

        # connector34_p2_fun = connector_34_out.p2_function
        # connector45_p2_fun = connector_45_out.p2_function

        # 3
        fp = Function_Part(p2_func4, start=agcct4_out.starting_point_in_ksi_phi_coordinate.x,
                           end=agcct4_out.end_point_in_ksi_phi_coordinate.x)
        # fp = Function_Part(p2_func4,start=agcct4_out.starting_point_in_ksi_phi_coordinate.x,end=agcct4_out.end_point_in_ksi_phi_coordinate.x)
        # fp = Function_Part(p2_func5,start=agcct5_out.starting_point_in_ksi_phi_coordinate.x,end=agcct5_out.end_point_in_ksi_phi_coordinate.x)

        print(fp.length)
        print(ksi0_4, ksi1_4)
        ksi_list = numpy.array([0, 13, 26, 40])*math.pi*2
        print(ksi_list)

        #############################

        def p2_func(ksi): return fp.valve_at(ksi)

        def p3_func(ksi): return lcs.point_to_global_coordinate(
            agcct3_out.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        ksi0 = 0
        ksi1 = fp.length

        tangential_direct = BaseUtils.derivative(p3_func)

        def main_normal_direct(ksi): return lcs.point_to_global_coordinate(
            agcct3_out.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        def second_normal_direction(ksi): return (
            tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        def path0(ksi): return p3_func(ksi)

        def path1(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path2(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path3(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        def path4(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        # numpy.savetxt('agcct_out_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,128*360)]))
        ksi0 = ksi_list[0]
        ksi1 = ksi_list[1]
        print(ksi0, ksi1)  # layer1_agcct1_in_1_part1
        numpy.savetxt('layer2_agcct2_out_2_part1.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct2_out_3_part1.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct2_out_4_part1.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct2_out_1_part1.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        print('down')

        ksi0 = ksi_list[1]
        ksi1 = ksi_list[2]
        print(ksi0, ksi1)
        numpy.savetxt('layer2_agcct2_out_1_part2.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct2_out_2_part2.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct2_out_3_part2.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct2_out_4_part2.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 13*360+1)]), fmt='%.4f')
        print('down')

        ksi0 = ksi_list[2]
        ksi1 = ksi_list[3]
        print(ksi0, ksi1)
        numpy.savetxt('layer2_agcct2_out_1_part3.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct2_out_2_part3.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct2_out_3_part3.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct2_out_4_part3.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 14*360+1)]), fmt='%.4f')
        print('down')

    if True:  # 四极CCT外层 5
        lcs = agcct3_out.local_coordinate_system.copy()
        lcs.location = P3.origin()
        print(lcs)
        # p2_func3 = lambda ksi:agcct3_out.p2_function(ksi)
        # p2_func4 = lambda ksi:agcct4_out.p2_function(ksi)
        def p2_func5(ksi): return agcct5_out.p2_function(ksi)

        # ksi0_3 = agcct3_out.starting_point_in_ksi_phi_coordinate.x
        # ksi0_4 = agcct4_out.starting_point_in_ksi_phi_coordinate.x
        ksi0_5 = agcct5_out.starting_point_in_ksi_phi_coordinate.x

        # ksi1_3 = agcct3_out.end_point_in_ksi_phi_coordinate.x
        # ksi1_4 = agcct4_out.end_point_in_ksi_phi_coordinate.x
        ksi1_5 = agcct5_out.end_point_in_ksi_phi_coordinate.x

        # connector34_p2_fun = connector_34_out.p2_function
        # connector45_p2_fun = connector_45_out.p2_function

        # 3
        fp = Function_Part(p2_func5, start=agcct5_out.starting_point_in_ksi_phi_coordinate.x,
                           end=agcct5_out.end_point_in_ksi_phi_coordinate.x)
        # fp = Function_Part(p2_func4,start=agcct4_out.starting_point_in_ksi_phi_coordinate.x,end=agcct4_out.end_point_in_ksi_phi_coordinate.x)
        # fp = Function_Part(p2_func5,start=agcct5_out.starting_point_in_ksi_phi_coordinate.x,end=agcct5_out.end_point_in_ksi_phi_coordinate.x)

        print(fp.length)
        print(ksi0_5, ksi1_5)
        ksi_list = numpy.array([0, 17, 34])*math.pi*2
        print(ksi_list)

        #############################

        def p2_func(ksi): return fp.valve_at(ksi)

        def p3_func(ksi): return lcs.point_to_global_coordinate(
            agcct3_out.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        ksi0 = 0
        ksi1 = fp.length

        tangential_direct = BaseUtils.derivative(p3_func)

        def main_normal_direct(ksi): return lcs.point_to_global_coordinate(
            agcct3_out.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        def second_normal_direction(ksi): return (
            tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        def path0(ksi): return p3_func(ksi)

        def path1(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path2(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)

        def path3(ksi): return p3_func(ksi) - 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        def path4(ksi): return p3_func(ksi) + 0.5*depth * \
            main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        # numpy.savetxt('agcct_out_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,128*360)]))
        ksi0 = ksi_list[0]
        ksi1 = ksi_list[1]
        print(ksi0, ksi1)  # layer1_agcct1_in_1_part1
        numpy.savetxt('layer2_agcct3_out_1_part1.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct3_out_2_part1.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct3_out_3_part1.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct3_out_4_part1.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        print('down')

        ksi0 = ksi_list[1]
        ksi1 = ksi_list[2]
        print(ksi0, ksi1)
        numpy.savetxt('layer2_agcct3_out_1_part2.txt', 1000*numpy.array([P3.as_p3(path1(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct3_out_2_part2.txt', 1000*numpy.array([P3.as_p3(path2(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct3_out_3_part2.txt', 1000*numpy.array([P3.as_p3(path3(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        numpy.savetxt('layer2_agcct3_out_4_part2.txt', 1000*numpy.array([P3.as_p3(path4(
            t)).to_list() for t in BaseUtils.linspace(ksi0, ksi1, 17*360+1)]), fmt='%.4f')
        print('down')
