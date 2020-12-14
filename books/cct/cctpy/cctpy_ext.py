"""
CCTPY 辅助功能，一般不必使用
"""

import math
from cctpy import *

try:
    from books.cct.cctpy.cctpy import *
except ModuleNotFoundError:
    pass


class Magnets(Magnet):
    """
    多个磁铁
    """

    def __init__(self, *magnet) -> None:
        self.magnets: List[Magnet] = list(magnet)

    def magnetic_field_at(self, point: P3) -> P3:
        b = P3()
        for m in self.magnets:
            b += m.magnetic_field_at(point)

        return b


class Line3:
    """
    三维空间曲线，带有起点、长度
    一般用于任意线圈的构建、CCT 洛伦兹力的分析等
    """

    def __init__(self, p3_function: Callable[[float], P3], start: float, end: float,
                 direct_function: Optional[Callable[[float], P3]] = None,
                 delta_for_compute_direct_function: float = 0.1*MM) -> None:
        """
        p3_function 曲线方程  p = p(s)
        start 曲线起点对应的自变量 s 值
        end 曲线终点对应的自变量 s 值
        direct_function 曲线方向方程， d = p'(s)，可以为空，若为空则 d = (p(s+Δ) - p(s))/Δ 计算而得
        delta_for_compute_direct_function 计算曲线方向方程 d(s) 时 Δ 取值，默认 0.1 毫米
        """
        self.p3_function = p3_function
        self.start = start
        self.end = end

        if direct_function is None:
            def direct_function(s) -> P3:
                return (
                    self.p3_function(s+delta_for_compute_direct_function) -
                    self.p3_function(s)
                )/delta_for_compute_direct_function

        self.direct_function = direct_function

    def point_at(self, s: float) -> P3:
        return self.p3_function(s)

    def direct_at(self, s: float) -> P3:
        return self.direct_function(s)

    def plot3(self, describe: str = 'r-', number: int = 1000) -> None:
        Plot3.plot_p3s([
            self.point_at(s) for s in BaseUtils.linspace(self.start, self.end, number)
        ], describe=describe)


class Wire(Magnet):
    """
    任意空间三维导线
    """

    def __init__(self, line3: Line3, current: float, delta_length: float = 1*MM) -> None:
        """
        line3 导线路径
        current 电流
        delta_length 导线分割成电流元时，每个电流元的长度
        """
        self.line3 = line3
        self.current = current
        self.start = line3.start
        self.end = line3.end

        # 分段数目
        part_number = int(math.ceil(
            (self.end-self.start)/delta_length
        ))

        self.dispersed_s = numpy.array(
            BaseUtils.linspace(self.start, self.end, part_number+1))

        self.dispersed_path3 = numpy.array([
            self.line3.point_at(s).to_numpy_ndarry3() for s in self.dispersed_s
        ])

        # 电流元 (miu0/4pi) * current * (p[i+1] - p[i])
        # refactor v0.1.1
        # 语法分析：示例
        # a = array([1, 2, 3, 4])
        # a[1:] = array([2, 3, 4])
        # a[:-1] = array([1, 2, 3])
        self.elementary_current = 1e-7 * self.current * (
            self.dispersed_path3[1:] - self.dispersed_path3[:-1]
        )

        # 电流元的位置 (p[i+1]+p[i])/2
        self.elementary_current_position = 0.5 * (
            self.dispersed_path3[1:] + self.dispersed_path3[:-1]
        )

    def magnetic_field_at(self, point: P3) -> P3:
        """
        计算磁场，全局坐标
        """
        if BaseUtils.equal(self.current,0,err=1e-6):
            return P3.zeros()
        
        if BaseUtils.equal(self.line3.start,self.line3.end,err=1e-6):
            return P3.zeros()

        p = point.to_numpy_ndarry3()

        # 点 p 到电流元中点
        r = p - self.elementary_current_position

        # 点 p 到电流元中点的距离的三次方
        rr = (numpy.linalg.norm(r, ord=2, axis=1)
              ** (-3)).reshape((r.shape[0], 1))

        # 计算每个电流元在 p 点产生的磁场 (此时还没有乘系数 μ0/4π )
        dB = numpy.cross(self.elementary_current, r) * rr

        # 求和，即得到磁场，
        # (不用乘乘以系数 μ0/4π = 1e-7)
        # refactor v0.1.1
        B = numpy.sum(dB, axis=0)

        # 转回 P3
        B_P3: P3 = P3.from_numpy_ndarry(B)

        return B_P3

    def magnetic_field_on_wire(self, s: float, delta_length: float,
                               local_coordinate_point: LocalCoordinateSystem,
                               other_magnet: Magnet = Magnet.no_magnet()
                               ) -> Tuple[P3, P3]:
        """
        计算线圈上磁场，
        s 线圈位置，即要计算磁场的位置
        delta_length 分段长度
        local_coordinate_point 局部坐标系
        """
        p = self.line3.point_at(s)
        pre_line3 = Line3(
            p3_function=self.line3.p3_function,
            start=self.line3.start,
            end=s - delta_length/2
        )

        post_line3 = Line3(
            p3_function=self.line3.p3_function,
            start=s + delta_length/2,
            end=self.line3.end
        )

        pre_wire = Wire(
            line3=pre_line3,
            current=self.current,
            delta_length=delta_length
        )

        post_wire = Wire(
            line3=post_line3,
            current=self.current,
            delta_length=delta_length
        )

        B = pre_wire.magnetic_field_at(
            p) + post_wire.magnetic_field_at(p) + other_magnet.magnetic_field_at(p)

        return p, local_coordinate_point.vector_to_local_coordinate(B)

    def lorentz_force_on_wire(self, s: float, delta_length: float,
                              local_coordinate_point: LocalCoordinateSystem,
                              other_magnet: Magnet = Magnet.no_magnet()
                              ) -> Tuple[P3, P3]:
        p, b = self.magnetic_field_on_wire(
            s=s,
            delta_length=delta_length,
            local_coordinate_point=LocalCoordinateSystem.global_coordinate_system(),
            other_magnet=other_magnet
        )

        direct = self.line3.direct_at(s)

        F = self.current * delta_length * (direct@b)

        return p, local_coordinate_point.vector_to_local_coordinate(F)

    @staticmethod
    def create_by_cct(cct: CCT) -> 'Wire':
        """
        由 CCT 创建 wire
        """
        return Wire(
            line3=Line3(
                p3_function=cct.p3_function,
                start=cct.starting_point_in_ksi_phi_coordinate.x,
                end=cct.end_point_in_ksi_phi_coordinate.x
            ),
            current=cct.current,
            delta_length=cct.small_r *
            BaseUtils.angle_to_radian(360/cct.winding_number)
        )


#########################################

R = 0.95
rin = 113*MM
rout = 128*MM
tain = [30]
taout = [150]
wn = 128
I = 9206
benda = 67.5
bendr = BaseUtils.angle_to_radian(benda)
delta_angle = 10

cct_in = CCT(
    local_coordinate_system=LocalCoordinateSystem.global_coordinate_system(),
    big_r=R, small_r=rin, bending_angle=benda,
    tilt_angles=tain, winding_number=wn,
    current=I, starting_point_in_ksi_phi_coordinate=P2.origin(),
    end_point_in_ksi_phi_coordinate=P2(
        wn*2*math.pi, bendr),
    disperse_number_per_winding=36
)

cct_out = CCT(
    local_coordinate_system=LocalCoordinateSystem.global_coordinate_system(),
    big_r=R, small_r=rout, bending_angle=benda,
    tilt_angles=taout, winding_number=wn,
    current=I, starting_point_in_ksi_phi_coordinate=P2.origin(),
    end_point_in_ksi_phi_coordinate=P2(
        -wn*2*math.pi, bendr),
    disperse_number_per_winding=36
)

wcct_in = Wire.create_by_cct(cct_in)
wcct_out = Wire.create_by_cct(cct_out)

# 查看内层洛伦兹力


def task(s):
    fon = wcct_in.lorentz_force_on_wire(
        s=BaseUtils.angle_to_radian(s),
        delta_length=rin * BaseUtils.angle_to_radian(delta_angle),
        # local_coordinate_point=LocalCoordinateSystem.global_coordinate_system(),
        local_coordinate_point=LocalCoordinateSystem(
            location=wcct_in.line3.point_at(BaseUtils.angle_to_radian(s)),
            x_direction=wcct_in.line3.direct_at(
                BaseUtils.angle_to_radian(s)),
            z_direction=cct_in.bipolar_toroidal_coordinate_system.main_normal_direction_at(
                cct_in.p2_function(BaseUtils.angle_to_radian(s))
            )
        ),
        other_magnet=cct_out
    )
    print(fon)
    return fon


if __name__ == "__main__":
    BaseUtils.i_am_sure_my_code_closed_in_if_name_equal_main()
    # 测试
    import unittest

    class Test(unittest.TestCase):
        def test_line3(self):
            line = Line3(lambda s: P3(s, 2*s, 0), 0, 1)
            print(line.point_at(0.5))
            print(line.direct_at(0.5))

            self.assertEqual(line.point_at(0.5), P3(0.5, 1, 0))
            self.assertEqual(line.direct_at(0.5), P3(1, 2, 0))

        def test_cct(self):
            r = 137*MM
            R = 0.95
            wn = 128
            cct = CCT(
                local_coordinate_system=LocalCoordinateSystem.global_coordinate_system(),
                big_r=R, small_r=r, bending_angle=67.5, tilt_angles=[30],
                winding_number=wn, current=9664, starting_point_in_ksi_phi_coordinate=P2(),
                end_point_in_ksi_phi_coordinate=P2(
                    wn*2*numpy.pi, BaseUtils.angle_to_radian(67.5)),
                disperse_number_per_winding=120
            )

            traj = (
                Trajectory.set_start_point(
                    P2(R, -1)).first_line(P2.y_direct(), 1.0)
                .add_arc_line(R, False, 67.5).add_strait_line(1.0)
            )

            # [0.743381541944373, -1.1125490995331133, 1.2364386785183252]
            B = cct.magnetic_field_at(
                traj.point_at(traj.get_length()/2).to_p3())
            print('cct', B)
            ksi = 10
            print(cct.p3_function(ksi))

        def test_wire(self):
            small_r = 137*MM
            big_r = 0.95
            wn = 128
            a = math.sqrt(big_r ** 2 - small_r ** 2)
            eta = 0.5 * math.log((big_r + a) / (big_r - a))
            tilt_angle = 30
            phi0 = BaseUtils.angle_to_radian(67.5)/wn

            def ksi_func(ksi):
                cota = 1/math.tan(BaseUtils.angle_to_radian(tilt_angle))
                return cota / math.sinh(eta) * math.sin(ksi) + phi0*ksi/(2*math.pi)

            def p2_2_p3(p: P2):
                nonlocal a, eta
                ksi = p.x
                phi = p.y
                temp = a / (math.cosh(eta) - math.cos(ksi))
                return P3(
                    temp * math.sinh(eta) * math.cos(phi),
                    temp * math.sinh(eta) * math.sin(phi),
                    temp * math.sin(ksi),
                )

            wire = Wire(
                line3=Line3(
                    lambda ksi: p2_2_p3(P2(ksi, ksi_func(ksi))),
                    0,
                    wn*2*math.pi
                ),
                current=9664,
                delta_length=BaseUtils.angle_to_radian(3)
            )

            traj = (
                Trajectory.set_start_point(
                    P2(0.95, -1)).first_line(P2.y_direct(), 1.0)
                .add_arc_line(big_r, False, 67.5).add_strait_line(1.0)
            )

            p = traj.point_at(traj.get_length()/2).to_p3()

            B = wire.magnetic_field_at(p)

            print('wire', B)

            ksi = 10
            print(p2_2_p3(P2(ksi, ksi_func(ksi))))

    # unittest.main(verbosity=1)

    # 洛伦兹力计算
    if False:

        # ss = BaseUtils.linspace(5, 128*360-5, 128*36)
        ss = BaseUtils.linspace(5+64*360, 360-5+64*360, 36)

        fons = BaseUtils.submit_process_task(task=task, param_list=[
            [s]
            for s in ss])

        data = []
        for i in range(len(fons)):
            p, f = fons[i]
            data.append([i+1, p.x, p.y, p.z, f.x, f.y, f.z])

        data = numpy.array(data)
        Plot2.plot_ndarry2ds(data[:, (0, 4)], describe='r-')
        Plot2.plot_ndarry2ds(data[:, (0, 5)], describe='b-')
        Plot2.plot_ndarry2ds(data[:, (0, 6)], describe='y-')
        Plot2.show()

        # Plot3.plot_cct(cct,describe='b-')
        # Plot3.plot_ndarry3ds(data[:,(1,2,3)],describe='r.')
        # Plot3.plot_xyz(data[19,1],data[19,2],data[19,3],describe='k.')

        # Plot3.show()

    if True:
        ss = BaseUtils.linspace(delta_angle/2, 360-delta_angle/2, 36)
        # ss = BaseUtils.linspace(delta_angle/2, 360-delta_angle/2, 360)

        BaseUtils.i_am_sure_my_code_closed_in_if_name_equal_main()
        fons = BaseUtils.submit_process_task(
            task=task, param_list=[[s] for s in ss])

        data = []
        for i in range(len(fons)):
            p, f = fons[i]
            data.append([i+1, p.x, p.y, p.z, f.x, f.y, f.z])

        data = numpy.array(data)

        # numpy.savetxt('data.txt',data)

        # data = numpy.loadtxt(r'C:\Users\madoka_9900\Documents\github\madokast.github.io\data.txt')
        
        if True: # 画图
            Plot2.plot_ndarry2ds(data[:, (0, 4)], describe='r-')
            Plot2.plot_ndarry2ds(data[:, (0, 5)], describe='b-')
            Plot2.plot_ndarry2ds(data[:, (0, 6)], describe='y-')

            Plot2.legend('绕线方向', 'rib方向', '径向', font_size=18,
                        font_family="Microsoft YaHei")
            Plot2.info('index', 'lorentz_force/N', '双层二极CCT的内层受力',
                    font_size=18, font_family="Microsoft YaHei")
            Plot2.show()
