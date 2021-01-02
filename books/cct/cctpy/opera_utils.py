"""
opera 相关工具类

@Author 赵润晓
"""

from typing import List
from cctpy import *
from cctpy_ext import *

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
            return tangential_direct(ksi)@main_normal_direct(ksi)

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


if __name__ == "__main__":
    # 2020 年参数
    # data = [-8.085,73.808,80.988,94.383,91.650,106.654,67.901,90.941,9488.615,-7334.914,24,46,37]

    # 2021.1.1 参数
    data = [5.573, 	-44.622 ,	87.453 ,	92.142, 	90.667, 	94.344,
     	73.471 ,	82.190 	,9426.734 ,	-5626.101 ,	25.000, 	40.000 ,	34.000]

    gantry = HUST_SC_GANTRY(
        qs3_gradient=data[0],
        qs3_second_gradient=data[1],
        dicct345_tilt_angles=[30, data[2], data[3], data[4]],
        agcct345_tilt_angles=[data[5] , 30, data[6], data[7]],
        dicct345_current=data[8],
        agcct345_current=data[9],
        agcct3_winding_number=data[10],
        agcct4_winding_number=data[11],
        agcct5_winding_number=data[12],
        agcct3_bending_angle = -67.5*(data[10])/(data[10]+data[11]+data[12]),
        agcct4_bending_angle = -67.5*(data[11])/(data[10]+data[11]+data[12]),
        agcct5_bending_angle = -67.5*(data[12])/(data[10]+data[11]+data[12]),
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

    b8s_list = [Brick8s.create_by_cct(c,3.2*MM,11*MM,'dicct',10) for c in diccts]
    b8s_list.extend([Brick8s.create_by_cct(c,3.2*MM,11*MM,'agcct',10) for c in agccts])

    operafile = open("opera0101.txt", "w")
    operafile.write(OperaConductor.to_opera_cond_script(b8s_list))
    operafile.close()
