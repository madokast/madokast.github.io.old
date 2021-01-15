import os
import sys
curPath = os.path.abspath(os.path.dirname(__file__))
rootPath = os.path.split(curPath)[0]
PathProject = os.path.split(rootPath)[0]
sys.path.append(rootPath)
sys.path.append(PathProject)

from cctpy import *
from cctpy_ext import *
from agcct_connector import *



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
            small_rs=[139.5*MM, 123.5*MM, 107.5*MM, 92.5*MM],
            bending_angles=[-17.05, -27.27, -23.18],  # [15.14, 29.02, 23.34]
            tilt_angles=[[30, 87.076, 91.829, 85.857],
                         [101.317, 30, 75.725, 92.044]],
            winding_numbers=[[128], [25, 40, 34]],
            currents=[9536.310, -6259.974],
            disperse_number_per_winding=36
        ).append_drift(BaseUtils.angle_to_radian(20)*R)
    )

    print([139.5*MM, 123.5*MM, 107.5*MM, 92.5*MM])

    depth, width = 11*MM, 3.2*MM

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

    if False: # 二极CCT内层
        lcs = dicct_in.local_coordinate_system.copy()
        lcs.location = P3.origin()
        print(lcs)
        p2_func = lambda ksi:dicct_in.p2_function(ksi)
        p3_func = lambda ksi:lcs.point_to_global_coordinate(dicct_in.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        ksi0 = dicct_in.starting_point_in_ksi_phi_coordinate.x
        ksi1 = dicct_in.end_point_in_ksi_phi_coordinate.x
        print(ksi0,ksi1)
        ksi_list = BaseUtils.linspace(ksi0,ksi1,8+1)
        print(ksi_list)

        tangential_direct = BaseUtils.derivative(p3_func)
        main_normal_direct = lambda ksi:lcs.point_to_global_coordinate(dicct_in.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        second_normal_direction = lambda ksi:(tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        path0 = lambda ksi:p3_func(ksi)
        path1 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path2 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path3 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        path4 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        # numpy.savetxt('dicct_in_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,128*360)]))
        for i in range(len(ksi_list)-1):
            print(i)
            ksi0 = ksi_list[i]
            ksi1 = ksi_list[i+1]
            print(ksi0,ksi1)
            numpy.savetxt(f'layer3_dicct_in_1_part_{i+1}.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,16*360+1)]))
            numpy.savetxt(f'layer3_dicct_in_2_part_{i+1}.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,16*360+1)]))
            numpy.savetxt(f'layer3_dicct_in_3_part_{i+1}.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,16*360+1)]))
            numpy.savetxt(f'layer3_dicct_in_4_part_{i+1}.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,16*360+1)]))
        # print('down')

        # Plot3.plot_p3s([p3_func(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='r-')
        # Plot3.plot_p3s([path1(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='b-')
        # Plot3.plot_p3s([path2(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='b-')
        # Plot3.plot_p3s([path3(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='b-')
        # Plot3.plot_p3s([path4(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='b-')
        # Plot3.set_center()
        # Plot3.show()
        

    if True: # 二极CCT外层
        lcs = dicct_out.local_coordinate_system.copy()
        lcs.location = P3.origin()
        print(lcs)
        p2_func = lambda ksi:dicct_out.p2_function(ksi)
        p3_func = lambda ksi:lcs.point_to_global_coordinate(dicct_out.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        ksi0 = dicct_out.starting_point_in_ksi_phi_coordinate.x
        ksi1 = dicct_out.end_point_in_ksi_phi_coordinate.x
        print(ksi0,ksi1)
        ksi_list = BaseUtils.linspace(ksi0,ksi1,8+1)
        print(ksi_list)

        tangential_direct = BaseUtils.derivative(p3_func)
        main_normal_direct = lambda ksi:lcs.point_to_global_coordinate(dicct_out.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        second_normal_direction = lambda ksi:(tangential_direct(ksi)@main_normal_direct(ksi)).normalize()


        path0 = lambda ksi:p3_func(ksi)
        path1 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path2 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path3 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        path4 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        
        # numpy.savetxt('dicct_out_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,128*360)]))
        for i in range(len(ksi_list)-1):
            print(i)
            ksi0 = ksi_list[i]
            ksi1 = ksi_list[i+1]
            print(ksi0,ksi1)
            numpy.savetxt(f'layer4_dicct_out_1_part{i+1}.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,16*360+1)]),fmt='%.10f')
            numpy.savetxt(f'layer4_dicct_out_2_part{i+1}.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,16*360+1)]),fmt='%.10f')
            numpy.savetxt(f'layer4_dicct_out_3_part{i+1}.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,16*360+1)]),fmt='%.10f')
            numpy.savetxt(f'layer4_dicct_out_4_part{i+1}.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,16*360+1)]),fmt='%.10f')
        # print('down')

        # Plot3.plot_p3s([p3_func(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='r-')
        # Plot3.plot_p3s([path1(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='b-')
        # Plot3.plot_p3s([path2(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='b-')
        # Plot3.plot_p3s([path3(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='b-')
        # Plot3.plot_p3s([path4(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='b-')
        # Plot3.set_center()
        # Plot3.show()
        

    if False: # 四极CCT内层 3 匝
        lcs = agcct3_in.local_coordinate_system.copy()
        lcs.location=P3.origin()
        print(lcs)
        p2_func3 = lambda ksi:agcct3_in.p2_function(ksi)
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
        fp = Function_Part(p2_func3,start=agcct3_in.starting_point_in_ksi_phi_coordinate.x,end=agcct3_in.end_point_in_ksi_phi_coordinate.x)
        print(fp.length)
        print(ksi0_3,ksi1_3)
        ksi_list = numpy.array([0,12,25])*math.pi*2
        print(ksi_list)

        #############################

        p2_func = lambda ksi:fp.valve_at(ksi)
        p3_func = lambda ksi:lcs.point_to_global_coordinate(agcct3_in.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        # ksi0 = 0
        # ksi1 = fp.length

        tangential_direct = BaseUtils.derivative(p3_func)
        main_normal_direct = lambda ksi:lcs.point_to_global_coordinate(agcct3_in.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        second_normal_direction = lambda ksi:(tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        path0 = lambda ksi:p3_func(ksi)
        path1 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path2 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path3 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        path4 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        
        # numpy.savetxt('agcct_in_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,25*360)]))
        ksi0 = ksi_list[0]
        ksi1 = ksi_list[1]
        print(ksi0,ksi1)
        numpy.savetxt('layer1_agcct1_in_1_part1.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,12*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct1_in_2_part1.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,12*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct1_in_3_part1.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,12*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct1_in_4_part1.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,12*360+1)]),fmt='%.10f')

        ksi0 = ksi_list[1]
        ksi1 = ksi_list[2]
        print(ksi0,ksi1)
        numpy.savetxt('layer1_agcct1_in_1_part2.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct1_in_2_part2.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct1_in_3_part2.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct1_in_4_part2.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        # print('down')


    if False: # 四极CCT内层 4 匝
        lcs = agcct3_in.local_coordinate_system.copy()
        lcs.location=P3.origin()
        print(lcs)
        # p2_func3 = lambda ksi:agcct3_in.p2_function(ksi)
        p2_func4 = lambda ksi:agcct4_in.p2_function(ksi)
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
        fp = Function_Part(p2_func4,start=agcct4_in.starting_point_in_ksi_phi_coordinate.x,end=agcct4_in.end_point_in_ksi_phi_coordinate.x)
        print(fp.length)
        print(ksi0_4,ksi1_4)
        ksi_list = numpy.array([0,13,26,40])*math.pi*2
        print(ksi_list)

        #############################

        p2_func = lambda ksi:fp.valve_at(ksi)
        p3_func = lambda ksi:lcs.point_to_global_coordinate(agcct3_in.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        # ksi0 = 0
        # ksi1 = fp.length

        tangential_direct = BaseUtils.derivative(p3_func)
        main_normal_direct = lambda ksi:lcs.point_to_global_coordinate(agcct3_in.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        second_normal_direction = lambda ksi:(tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        path0 = lambda ksi:p3_func(ksi)
        path1 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path2 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path3 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        path4 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        
        # numpy.savetxt('agcct_in_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,25*360)]))
        ksi0 = ksi_list[0]
        ksi1 = ksi_list[1]
        print(ksi0,ksi1)
        numpy.savetxt('layer1_agcct2_in_1_part1.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct2_in_2_part1.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct2_in_3_part1.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct2_in_4_part1.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')

        ksi0 = ksi_list[1]
        ksi1 = ksi_list[2]
        print(ksi0,ksi1)
        numpy.savetxt('layer1_agcct2_in_1_part2.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct2_in_2_part2.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct2_in_3_part2.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct2_in_4_part2.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')

        ksi0 = ksi_list[2]
        ksi1 = ksi_list[3]
        print(ksi0,ksi1)
        numpy.savetxt('layer1_agcct2_in_1_part3.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,14*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct2_in_2_part3.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,14*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct2_in_3_part3.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,14*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct2_in_4_part3.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,14*360+1)]),fmt='%.10f')


    if False: # 四极CCT内层 5 匝
        lcs = agcct3_in.local_coordinate_system.copy()
        lcs.location=P3.origin()
        print(lcs)
        # p2_func3 = lambda ksi:agcct3_in.p2_function(ksi)
        # p2_func4 = lambda ksi:agcct4_in.p2_function(ksi)
        p2_func5 = lambda ksi:agcct5_in.p2_function(ksi)

        # ksi0_3 = agcct3_in.starting_point_in_ksi_phi_coordinate.x
        # ksi0_4 = agcct4_in.starting_point_in_ksi_phi_coordinate.x
        ksi0_5 = agcct5_in.starting_point_in_ksi_phi_coordinate.x

        # ksi1_3 = agcct3_in.end_point_in_ksi_phi_coordinate.x
        # ksi1_4 = agcct4_in.end_point_in_ksi_phi_coordinate.x
        ksi1_5 = agcct5_in.end_point_in_ksi_phi_coordinate.x

        # connector34_p2_fun = connector_34_in.p2_function
        # connector45_p2_fun = connector_45_in.p2_function

        # 3
        fp = Function_Part(p2_func5,start=agcct5_in.starting_point_in_ksi_phi_coordinate.x,end=agcct5_in.end_point_in_ksi_phi_coordinate.x)
        print(fp.length)
        print(ksi0_5,ksi1_5)
        ksi_list = numpy.array([0,17,34])*math.pi*2
        print(ksi_list)

        #############################

        p2_func = lambda ksi:fp.valve_at(ksi)
        p3_func = lambda ksi:lcs.point_to_global_coordinate(agcct3_in.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        # ksi0 = 0
        # ksi1 = fp.length

        tangential_direct = BaseUtils.derivative(p3_func)
        main_normal_direct = lambda ksi:lcs.point_to_global_coordinate(agcct3_in.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        second_normal_direction = lambda ksi:(tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        path0 = lambda ksi:p3_func(ksi)
        path1 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path2 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path3 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        path4 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        
        # numpy.savetxt('agcct_in_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,25*360)]))
        ksi0 = ksi_list[0]
        ksi1 = ksi_list[1]
        print(ksi0,ksi1)
        numpy.savetxt('layer1_agcct3_in_1_part1.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct3_in_2_part1.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct3_in_3_part1.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct3_in_4_part1.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')

        ksi0 = ksi_list[1]
        ksi1 = ksi_list[2]
        print(ksi0,ksi1)
        numpy.savetxt('layer1_agcct3_in_1_part2.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct3_in_2_part2.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct3_in_3_part2.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        numpy.savetxt('layer1_agcct3_in_4_part2.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')


    if False: # 四极CCT外层 3
        lcs = agcct3_out.local_coordinate_system.copy()
        lcs.location=P3.origin()
        print(lcs)
        p2_func3 = lambda ksi:agcct3_out.p2_function(ksi)
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
        fp = Function_Part(p2_func3,start=agcct3_out.starting_point_in_ksi_phi_coordinate.x,end=agcct3_out.end_point_in_ksi_phi_coordinate.x)
        # fp = Function_Part(p2_func4,start=agcct4_out.starting_point_in_ksi_phi_coordinate.x,end=agcct4_out.end_point_in_ksi_phi_coordinate.x)
        # fp = Function_Part(p2_func5,start=agcct5_out.starting_point_in_ksi_phi_coordinate.x,end=agcct5_out.end_point_in_ksi_phi_coordinate.x)

        print(fp.length)
        print(ksi0_3,ksi1_3)
        ksi_list = numpy.array([0,12,25])*math.pi*2
        print(ksi_list)

        #############################

        p2_func = lambda ksi:fp.valve_at(ksi)
        p3_func = lambda ksi:lcs.point_to_global_coordinate(agcct3_out.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        ksi0 = 0
        ksi1 = fp.length

        tangential_direct = BaseUtils.derivative(p3_func)
        main_normal_direct = lambda ksi:lcs.point_to_global_coordinate(agcct3_out.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        second_normal_direction = lambda ksi:(tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        path0 = lambda ksi:p3_func(ksi)
        path1 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path2 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path3 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        path4 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        
        # numpy.savetxt('agcct_out_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,128*360)]))
        ksi0 = ksi_list[0]
        ksi1 = ksi_list[1]
        print(ksi0,ksi1) # layer1_agcct1_in_1_part1
        numpy.savetxt('layer2_agcct1_out_1_part1.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,12*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct1_out_2_part1.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,12*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct1_out_3_part1.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,12*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct1_out_4_part1.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,12*360+1)]),fmt='%.10f')
        print('down')

        ksi0 = ksi_list[1]
        ksi1 = ksi_list[2]
        print(ksi0,ksi1)
        numpy.savetxt('layer2_agcct1_out_1_part2.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct1_out_2_part2.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct1_out_3_part2.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct1_out_4_part2.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        print('down')

        # Plot3.plot_p3s([p3_func(t) for t in BaseUtils.linspace(ksi0,ksi1,128*3600)],describe='r-')
        # Plot3.plot_p3s([path1(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='b-')
        # Plot3.plot_p3s([path2(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='b-')
        # Plot3.plot_p3s([path3(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='b-')
        # Plot3.plot_p3s([path4(t) for t in BaseUtils.linspace(ksi0,ksi1,10000)],describe='b-')
        # Plot3.set_center()
        # Plot3.show()


    if False: # 四极CCT外层 4
        lcs = agcct3_out.local_coordinate_system.copy()
        lcs.location=P3.origin()
        print(lcs)
        # p2_func3 = lambda ksi:agcct3_out.p2_function(ksi)
        p2_func4 = lambda ksi:agcct4_out.p2_function(ksi)
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
        fp = Function_Part(p2_func4,start=agcct4_out.starting_point_in_ksi_phi_coordinate.x,end=agcct4_out.end_point_in_ksi_phi_coordinate.x)
        # fp = Function_Part(p2_func4,start=agcct4_out.starting_point_in_ksi_phi_coordinate.x,end=agcct4_out.end_point_in_ksi_phi_coordinate.x)
        # fp = Function_Part(p2_func5,start=agcct5_out.starting_point_in_ksi_phi_coordinate.x,end=agcct5_out.end_point_in_ksi_phi_coordinate.x)

        print(fp.length)
        print(ksi0_4,ksi1_4)
        ksi_list = numpy.array([0,13,26,40])*math.pi*2
        print(ksi_list)

        #############################

        p2_func = lambda ksi:fp.valve_at(ksi)
        p3_func = lambda ksi:lcs.point_to_global_coordinate(agcct3_out.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        ksi0 = 0
        ksi1 = fp.length

        tangential_direct = BaseUtils.derivative(p3_func)
        main_normal_direct = lambda ksi:lcs.point_to_global_coordinate(agcct3_out.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        second_normal_direction = lambda ksi:(tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        path0 = lambda ksi:p3_func(ksi)
        path1 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path2 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path3 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        path4 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        
        # numpy.savetxt('agcct_out_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,128*360)]))
        ksi0 = ksi_list[0]
        ksi1 = ksi_list[1]
        print(ksi0,ksi1) # layer1_agcct1_in_1_part1
        numpy.savetxt('layer2_agcct2_out_1_part1.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct2_out_2_part1.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct2_out_3_part1.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct2_out_4_part1.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        print('down')

        ksi0 = ksi_list[1]
        ksi1 = ksi_list[2]
        print(ksi0,ksi1)
        numpy.savetxt('layer2_agcct2_out_1_part2.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct2_out_2_part2.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct2_out_3_part2.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct2_out_4_part2.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,13*360+1)]),fmt='%.10f')
        print('down')

        ksi0 = ksi_list[2]
        ksi1 = ksi_list[3]
        print(ksi0,ksi1)
        numpy.savetxt('layer2_agcct2_out_1_part3.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,14*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct2_out_2_part3.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,14*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct2_out_3_part3.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,14*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct2_out_4_part3.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,14*360+1)]),fmt='%.10f')
        print('down')


    if False: # 四极CCT外层 5
        lcs = agcct3_out.local_coordinate_system.copy()
        lcs.location=P3.origin()
        print(lcs)
        # p2_func3 = lambda ksi:agcct3_out.p2_function(ksi)
        # p2_func4 = lambda ksi:agcct4_out.p2_function(ksi)
        p2_func5 = lambda ksi:agcct5_out.p2_function(ksi)

        # ksi0_3 = agcct3_out.starting_point_in_ksi_phi_coordinate.x
        # ksi0_4 = agcct4_out.starting_point_in_ksi_phi_coordinate.x
        ksi0_5 = agcct5_out.starting_point_in_ksi_phi_coordinate.x

        # ksi1_3 = agcct3_out.end_point_in_ksi_phi_coordinate.x
        # ksi1_4 = agcct4_out.end_point_in_ksi_phi_coordinate.x
        ksi1_5 = agcct5_out.end_point_in_ksi_phi_coordinate.x

        # connector34_p2_fun = connector_34_out.p2_function
        # connector45_p2_fun = connector_45_out.p2_function

        # 3
        fp = Function_Part(p2_func5,start=agcct5_out.starting_point_in_ksi_phi_coordinate.x,end=agcct5_out.end_point_in_ksi_phi_coordinate.x)
        # fp = Function_Part(p2_func4,start=agcct4_out.starting_point_in_ksi_phi_coordinate.x,end=agcct4_out.end_point_in_ksi_phi_coordinate.x)
        # fp = Function_Part(p2_func5,start=agcct5_out.starting_point_in_ksi_phi_coordinate.x,end=agcct5_out.end_point_in_ksi_phi_coordinate.x)

        print(fp.length)
        print(ksi0_5,ksi1_5)
        ksi_list = numpy.array([0,17,34])*math.pi*2
        print(ksi_list)

        #############################

        p2_func = lambda ksi:fp.valve_at(ksi)
        p3_func = lambda ksi:lcs.point_to_global_coordinate(agcct3_out.bipolar_toroidal_coordinate_system.convert(p2_func(ksi)))
        ksi0 = 0
        ksi1 = fp.length

        tangential_direct = BaseUtils.derivative(p3_func)
        main_normal_direct = lambda ksi:lcs.point_to_global_coordinate(agcct3_out.bipolar_toroidal_coordinate_system.main_normal_direction_at(p2_func(ksi))).normalize()
        second_normal_direction = lambda ksi:(tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        path0 = lambda ksi:p3_func(ksi)
        path1 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path2 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) + 0.5*width*second_normal_direction(ksi)
        path3 = lambda ksi:p3_func(ksi) - 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)
        path4 = lambda ksi:p3_func(ksi) + 0.5*depth*main_normal_direct(ksi) - 0.5*width*second_normal_direction(ksi)

        
        # numpy.savetxt('agcct_out_center.txt',1000*numpy.array([P3.as_p3(path0(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,128*360)]))
        ksi0 = ksi_list[0]
        ksi1 = ksi_list[1]
        print(ksi0,ksi1) # layer1_agcct1_in_1_part1
        numpy.savetxt('layer2_agcct3_out_1_part1.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct3_out_2_part1.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct3_out_3_part1.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct3_out_4_part1.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        print('down')

        ksi0 = ksi_list[1]
        ksi1 = ksi_list[2]
        print(ksi0,ksi1)
        numpy.savetxt('layer2_agcct3_out_1_part2.txt',1000*numpy.array([P3.as_p3(path1(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct3_out_2_part2.txt',1000*numpy.array([P3.as_p3(path2(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct3_out_3_part2.txt',1000*numpy.array([P3.as_p3(path3(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        numpy.savetxt('layer2_agcct3_out_4_part2.txt',1000*numpy.array([P3.as_p3(path4(t)).to_list() for t in BaseUtils.linspace(ksi0,ksi1,17*360+1)]),fmt='%.10f')
        print('down')
