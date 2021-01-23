"""
opera 相关工具类
1. 将 CCT 转为 opera 中的 Brick8 导体块，并导出为 cond 文件
    关键方法 Brick8s.create_by_cct()、OperaConductor.to_opera_cond_script()

2. 读取 opera 磁场表格文件，生成 cctpy 中的 Magnet 对象


2021年1月4日 新增 opera 导出的磁场文件读入功能，读入磁场表格文件，即返回一个 Magnet
2021年1月4日 新增读入 opera track 轨迹 log 文件功能，返回一个 Line3

@Author 赵润晓
"""

try:
    from books.cct.cctpy.cctpy import *
    from books.cct.cctpy.cctpy_ext import *
except:
    pass

from typing import List

from cctpy import *
from cctpy_ext import *
import numpy

OPERA_CONDUCTOR_SCRIPT_HEAD: str = "CONDUCTOR\n"
OPERA_CONDUCTOR_SCRIPT_TAIL: str = "QUIT\nEOF\n"


class Brick8:
    """
    opera 中 8 点导线立方体
    对应 opera 中脚本：
    --------------------------------------------------
    DEFINE BR8
    0.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0
    0.0 0.0 0.0
    1.054000e+00 5.651710e-04 -8.249738e-04
    1.046000e+00 5.651710e-04 -8.249738e-04
    1.046000e+00 -5.651710e-04 8.249738e-04
    1.054000e+00 -5.651710e-04 8.249738e-04
    1.004041e+00 1.474080e-01 1.026480e-01
    9.981663e-01 1.465494e-01 9.728621e-02
    9.973407e-01 1.451229e-01 9.841917e-02
    1.003216e+00 1.459815e-01 1.037810e-01
    3.4575E8 1  'layer1'
    0 0 0
    1.0
    --------------------------------------------------
    脚本中各参数的意义如下：
    --------------------------------------------------
    XCENTRE, YCENTRE, ZCENTRE, PHI1, THETA1, PSI1           # Local coordinate system 1
    XCEN2, YCEN2, ZCEN2                                     # Local coordinate system 2 (origin)
    THETA2, PHI2, PSI2                                      # Local coordinate system 2 (Euler angles)
    XP1, YP1, ZP1                                           #  Bottom right corner of front face
    XP2, YP2, ZP2                                           #  Top right corner of front face
    XP3, YP3, ZP3                                           #  Top left corner of front face
    XP4, YP4, ZP4                                           #  Bottom left corner of front face
    XP5, YP5, ZP5                                           #  Bottom right corner of back face
    XP6, YP6, ZP6                                           #  Top right corner of back face
    XP7, YP7, ZP7                                           #  Top left corner of back face
    XP8, YP8, ZP8                                           #  Bottom left corner of back face
    CURD, SYMMETRY, DRIVELABEL                              #  Current density, symmetry and drive label
    IRXY, IRYZ,IRZX                                         #  Reflections in local coordinate system 1 coordinate planes
    TOLERANCE Flux                                          #  density tolerance
    --------------------------------------------------
    """

    HEAD = 'DEFINE BR8\n0.0 0.0 0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0\n'
    TAIL = '0 0 0\n1.0\n'

    def __init__(self,
                 front_face_point1: P3,
                 front_face_point2: P3,
                 front_face_point3: P3,
                 front_face_point4: P3,
                 back_face_point1: P3,
                 back_face_point2: P3,
                 back_face_point3: P3,
                 back_face_point4: P3,
                 current_density: float,
                 label: str
                 ) -> None:
        """
        front_face_point1234 立方体前面的四个点
        back_face_point1234 立方体后面的四个点

        所谓的前面/后面，指的是按照电流方向（电流从前面流入，后面流出，参考 opera ref-3d 手册）
        前面  -> 电流 -> 后面
        """
        self.front_face_point1 = front_face_point1
        self.front_face_point2 = front_face_point2
        self.front_face_point3 = front_face_point3
        self.front_face_point4 = front_face_point4
        self.back_face_point1 = back_face_point1
        self.back_face_point2 = back_face_point2
        self.back_face_point3 = back_face_point3
        self.back_face_point4 = back_face_point4
        self.current_density = current_density
        self.label = label

    def to_opera_cond(self) -> str:
        def p3_str(p: P3) -> str:
            return f'{p.x} {p.y} {p.z}\n'
        front_face_point1_str = p3_str(self.front_face_point1)
        front_face_point2_str = p3_str(self.front_face_point2)
        front_face_point3_str = p3_str(self.front_face_point3)
        front_face_point4_str = p3_str(self.front_face_point4)

        back_face_point1_str = p3_str(self.back_face_point1)
        back_face_point2_str = p3_str(self.back_face_point2)
        back_face_point3_str = p3_str(self.back_face_point3)
        back_face_point4_str = p3_str(self.back_face_point4)

        current_label_str = f"{self.current_density} 1 '{self.label}'\n"

        return "".join((
            Brick8.HEAD,
            front_face_point1_str,
            front_face_point2_str,
            front_face_point3_str,
            front_face_point4_str,

            back_face_point1_str,
            back_face_point2_str,
            back_face_point3_str,
            back_face_point4_str,

            current_label_str,
            Brick8.TAIL
        ))


class Brick8s:
    """
    opera 中连续的 8 点导体立方体
    所谓连续，指的是前后两个立方体，各有一个面重合
    """

    def __init__(self,
                 line1: List[P3],
                 line2: List[P3],
                 line3: List[P3],
                 line4: List[P3],
                 current_density: float,
                 label: str
                 ) -> None:
        self.line1 = line1
        self.line2 = line2
        self.line3 = line3
        self.line4 = line4
        self.current_density = current_density
        self.label = label

    def to_brick8(self) -> List[Brick8]:
        bricks_list = []
        size = len(self.line1)
        for i in range(size-1):
            bricks_list.append(Brick8(
                self.line1[i],
                self.line2[i],
                self.line3[i],
                self.line4[i],
                self.line1[i+1],
                self.line2[i+1],
                self.line3[i+1],
                self.line4[i+1],
                self.current_density,
                self.label
            ))

        return bricks_list

    def to_opera_cond(self) -> str:
        bricks_list = self.to_brick8()
        return "\n".join([e.to_opera_cond() for e in bricks_list])

    @staticmethod
    def create_by_cct(cct: CCT, channel_width: float, channel_depth: float,
                      label: str, disperse_number_per_winding: int) -> 'Brick8s':
        """
        从 CCT 创建 Brick8s
        channel_width channel_depth 槽的宽度和深度
        label 标签
        disperse_number_per_winding 每匝分段数目

        注意：转为 Brick8s 时，没有进行坐标转换，即在 CCT 的局部坐标系中建模
        """
        delta = 1e-6

        # 路径方程
        def path3(ksi):
            return cct.p3_function(ksi)

        # 切向 正则归一化
        def tangential_direct(ksi):
            return ((path3(ksi+delta)-path3(ksi))/delta).normalize()

        # 主法线方向 注意：已正则归一化
        def main_normal_direct(ksi):
            return cct.bipolar_toroidal_coordinate_system.main_normal_direction_at(cct.p2_function(ksi))

        # 副法线方向
        def second_normal_direc(ksi):
            return (tangential_direct(ksi)@main_normal_direct(ksi)).normalize()

        def channel_path1(ksi):
            return (path3(ksi)
                    + (channel_depth/2) * main_normal_direct(ksi)
                    + (channel_width/2) * second_normal_direc(ksi)
                    )

        def channel_path2(ksi):
            return (path3(ksi)
                    - (channel_depth/2) * main_normal_direct(ksi)
                    + (channel_width/2) * second_normal_direc(ksi)
                    )

        def channel_path3(ksi):
            return (path3(ksi)
                    - (channel_depth/2) * main_normal_direct(ksi)
                    - (channel_width/2) * second_normal_direc(ksi)
                    )

        def channel_path4(ksi):
            return (path3(ksi)
                    + (channel_depth/2) * main_normal_direct(ksi)
                    - (channel_width/2) * second_normal_direc(ksi)
                    )

        start_ksi = cct.starting_point_in_ksi_phi_coordinate.x
        end_ksi = cct.end_point_in_ksi_phi_coordinate.x
        # +1 为了 linspace 获得正确分段结果
        total_disperse_number = disperse_number_per_winding * cct.winding_number + 1

        ksi_list = BaseUtils.linspace(
            start_ksi, end_ksi, total_disperse_number)

        return Brick8s(
            [channel_path1(ksi) for ksi in ksi_list],
            [channel_path2(ksi) for ksi in ksi_list],
            [channel_path3(ksi) for ksi in ksi_list],
            [channel_path4(ksi) for ksi in ksi_list],
            current_density=cct.current / (channel_width*channel_depth),
            label=label
        )


class OperaConductor:
    @staticmethod
    def to_opera_cond_script(brick8s_list: List[Brick8s]) -> str:
        ret = [OPERA_CONDUCTOR_SCRIPT_HEAD]
        for b in brick8s_list:
            ret.append(b.to_opera_cond())

        ret.append(OPERA_CONDUCTOR_SCRIPT_TAIL)

        return '\n'.join(ret)


class OperaFieldTableMagnet(Magnet):
    """
    利用 opera 磁场表格文件生成的磁场
    文件主体应为 6 列，分别是 x y z bx by bz
    """

    def __init__(self, file_name: str,
                 first_corner_x: float, first_corner_y: float, first_corner_z: float,
                 step_between_points_x: float, step_between_points_y: float, step_between_points_z: float,
                 number_of_points_x: int, number_of_points_y: int, number_of_points_z: int,
                 unit_of_length: float = M, unit_of_field: float = 1
                 ) -> None:
        """
        first_corner_x / y / z 即 opera 导出磁场时 First corner 填写的值
        step_between_points_x / y / z 即 opera 导出磁场时 Step between points 填写的值
        number_of_points_x / y / z 即 opera 导出磁场时 Number of points 填写的值
        unit_of_length 是数据中长度单位，默认米。（注意 opera 中默认毫米，如果 opera 未修改单位，需要设为毫米 MM）
        unit_of_field 是磁场单位，默认特斯拉 T。

        以上所有参数，都可以利用磁场文件中前 8 行的元数据获知，但介于 python 性能一般就不分析了。（其实是懒）
        """
        # 打开磁场表格文件
        """
        5 461 311 2
        1 X [METRE]
        2 Y [METRE]
        3 Z [METRE]
        4 BX [TESLA]
        5 BY [TESLA]
        6 BZ [TESLA]
        0
        -0.550000000000      -1.10000000000     -0.100000000000E-01  0.500216216871E-04  0.431668985488E-05 -0.180396818407E-02
        -0.550000000000      -1.10000000000     -0.500000000000E-02  0.422312886145E-04 -0.426162472137E-05 -0.180415005970E-02
        ... ...
        """
        # 总点数
        self.x0 = first_corner_x*unit_of_length
        self.y0 = first_corner_y*unit_of_length
        self.z0 = first_corner_z*unit_of_length
        self.gap_x = step_between_points_x*unit_of_length
        self.gap_y = step_between_points_y*unit_of_length
        self.gap_z = step_between_points_z*unit_of_length
        self.number_of_points_x = int(number_of_points_x)
        self.number_of_points_y = int(number_of_points_y)
        self.number_of_points_z = int(number_of_points_z)
        self.total_point_number = self.number_of_points_x * \
            self.number_of_points_y*self.number_of_points_z
        print(f"table 文件包含 {self.total_point_number} 个节点")

        # 读取数据，使用 numpy
        data = numpy.loadtxt(fname=file_name, skiprows=8)  # 跳过 8 行，因为前 8 行非数据
        # 找出代表bx by bz 的列
        self.mxs = data[:, 3]
        self.mys = data[:, 4]
        self.mzs = data[:, 5]

        if abs(unit_of_field-1) < 1e-6:
            self.mxs = self.mxs * unit_of_field
            self.mys = self.mys * unit_of_field
            self.mzs = self.mzs * unit_of_field

    def magnetic_field_at(self, point: P3) -> P3:
        """
        核心方法
        """
        raise NotImplementedError


if __name__ == "__main__":
    if True: # 导出 cond 文件
        # 2020 年参数
        # data = [-8.085,73.808,80.988,94.383,91.650,106.654,67.901,90.941,9488.615,-7334.914,24,46,37]

        # 2021.1.1 参数
        # data = [4.378,-90.491,87.076,91.829,85.857,101.317,75.725,92.044,9536.310,-6259.974,25,40,34]

        data = [4.675 ,	41.126 	,88.773 ,	98.139 ,
        	91.748 	,101.792 ,	62.677 ,	89.705 ,
            	9409.261 ,	-7107.359 , 25, 40, 34]

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

            agcct345_inner_small_r=92.5 * MM, 
            agcct345_outer_small_r=108.5 * MM, 
            dicct345_inner_small_r=124.5 * MM,
            dicct345_outer_small_r=140.5 * MM,
        )

        bl_all = gantry.create_beamline()

        f = gantry.first_bending_part_length()

        sp = bl_all.trajectory.point_at(f)
        sd = bl_all.trajectory.direct_at(f)

        bl = gantry.create_second_bending_part(sp, sd)

        diccts = bl.magnets[0:2]
        agccts = bl.magnets[2:8]

        # m = Magnets(*ccts)

        # bz = m.magnetic_field_bz_along(bl.trajectory,step=10*MM)

        # Plot2.plot(bz)
        # Plot2.show()

        b8s_list = [Brick8s.create_by_cct(
            c, 3.2*MM, 11*MM, 'dicct', 10) for c in diccts]
        b8s_list.extend([Brick8s.create_by_cct(
            c, 3.2*MM, 11*MM, 'agcct', 10) for c in agccts])

        operafile = open("opera021.cond", "w")
        operafile.write(OperaConductor.to_opera_cond_script(b8s_list))
        operafile.close()


    if False:# opera 对比研究
        data = [4.378,-90.491,87.076,91.829,85.857,101.317,75.725,92.044,9536.310,-6259.974,25,40,34]

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

            agcct345_inner_small_r=83 * MM + 9.5*MM,
            agcct345_outer_small_r=98 * MM + 9.5*MM,  # 83+15
            dicct345_inner_small_r=114 * MM + 9.5*MM,  # 83+30+1
            dicct345_outer_small_r=130 * MM + 9.5*MM,  # 83+45 +2

            part_per_winding=360
        )

        bl_all = gantry.create_beamline()

        f = gantry.first_bending_part_length()

        sp = bl_all.trajectory.point_at(f)
        sd = bl_all.trajectory.direct_at(f)

        bl = gantry.create_second_bending_part(sp, sd)

        diccts = bl.magnets[0:2]
        agccts = bl.magnets[2:8]
        ccts = Magnets(*diccts).add_all(agccts)
        

        # 局部坐标系
        lcs = CCT.as_cct(diccts[0]).local_coordinate_system
        for cct in ccts.to_list():
            print(lcs == CCT.as_cct(cct).local_coordinate_system)
        
        origin_lcs = P3.origin()
        origin_gcs = lcs.point_to_global_coordinate(origin_lcs)
        print(origin_gcs)
        # [6.640302200026097, 2.7739506645753815, 0.0]

        arc = ArcLine2(
            starting_phi=BaseUtils.angle_to_radian(-120),
            center=P2.origin(),
            radius=0.95,
            total_phi=BaseUtils.angle_to_radian(180),
            clockwise=False
        )

        ms = []
        for s in BaseUtils.linspace(0,0.95*BaseUtils.angle_to_radian(180),721):
            p = arc.point_at(s).to_p3()
            p = lcs.point_to_global_coordinate(p) + P3(z=75*MM) # 上移动 75mm
            m = ccts.magnetic_field_at(p)
            print(BaseUtils.radian_to_angle(s/0.95),m.z)
            ms.append(P2(s,m.z))
        
        Plot2.plot(ms)
        Plot2.show()

        # arc_ps = [arc.point_at(s) for s in BaseUtils.linspace(0,BaseUtils.angle_to_radian(180)*0.95,181)]
