"""
CCTPY 辅助功能，一般不必使用

主要用于进行线圈洛伦兹力分析、压强分析
为此建立了 Line3 和 Wire，表示任意连续三维曲线及其构成的导线

本模块也用于 agcct_connector，构建 agcct 的连接段

2021年1月9日 新增 Function_Part ，用于构建连接后的 agcct

@Author 赵润晓
"""

import math
from typing import Any, TypeVar
from cctpy import *

_T = TypeVar("_T")

try:
    from books.cct.cctpy.cctpy import *
except:
    pass


class Magnets(Magnet):
    """
    多个磁铁，仅仅就是放进数组
    """

    def __init__(self, *magnet) -> None:
        self.magnets: List[Magnet] = list(magnet)

    def add(self, magnet: Magnet) -> 'Magnets':
        self.magnets.append(magnet)
        return self

    def add_all(self, magnets: List[Magnet]) -> 'Magnets':
        self.magnets.extend(magnets)
        return self

    def to_list(self)->List[Magnet]:
        return self.magnets

    def remove(self,  magnet: Magnet) -> 'Magnets':
        self.magnets.remove(magnet)
        return self

    def magnetic_field_at(self, point: P3) -> P3:
        b = P3()
        for m in self.magnets:
            b += m.magnetic_field_at(point)

        return b

    def __str__(self) -> str:
        return f"Magnets: 共有 {len(self.magnets)} 个磁铁"

    def __repr__(self) -> str:
        return self.__str__()


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
        if BaseUtils.equal(self.current, 0, err=1e-6):
            return P3.zeros()

        if BaseUtils.equal(self.line3.start, self.line3.end, err=1e-6):
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

        返回值 位置点+磁场
        """
        # 所在点（全局坐标系）
        p = self.line3.point_at(s)

        # 所在点之前的导线
        pre_line3 = Line3(
            p3_function=self.line3.p3_function,
            start=self.line3.start,
            end=s - delta_length/2
        )

        # 所在点之后的导线
        post_line3 = Line3(
            p3_function=self.line3.p3_function,
            start=s + delta_length/2,
            end=self.line3.end
        )

        # 所在点之前的导线
        pre_wire = Wire(
            line3=pre_line3,
            current=self.current,
            delta_length=delta_length
        )

        # 所在点之后的导线
        post_wire = Wire(
            line3=post_line3,
            current=self.current,
            delta_length=delta_length
        )

        # 磁场
        B = pre_wire.magnetic_field_at(
            p) + post_wire.magnetic_field_at(p) + other_magnet.magnetic_field_at(p)

        # 返回点 + 磁场
        return p, local_coordinate_point.vector_to_local_coordinate(B)

    def lorentz_force_on_wire(self, s: float, delta_length: float,
                              local_coordinate_point: LocalCoordinateSystem,
                              other_magnet: Magnet = Magnet.no_magnet()
                              ) -> Tuple[P3, P3]:
        """
        计算线圈上 s 位置处洛伦兹力
        delta_length 线圈分段长度（s位置的这一段线圈不参与计算）
        local_coordinate_point 洛伦兹力坐标系
        """
        p, b = self.magnetic_field_on_wire(
            s=s,
            delta_length=delta_length,
            # 我居然写多了，真厉害
            local_coordinate_point=LocalCoordinateSystem.global_coordinate_system(),
            other_magnet=other_magnet
        )

        direct = self.line3.direct_at(s)

        F = self.current * delta_length * (direct@b)

        return p, local_coordinate_point.vector_to_local_coordinate(F)

    def pressure_on_wire_MPa(self, s: float, delta_length: float,
                             local_coordinate_point: LocalCoordinateSystem,
                             channel_width: float, channel_depth: float,
                             other_magnet: Magnet = Magnet.no_magnet(),
                             ) -> Tuple[P3, P3]:
        """
        计算压强，默认为 CCT

        默认在 local_coordinate_point 坐标系下计算的洛伦兹力 F

        Fx 绕线方向（应该是 0 ）
        Fy rib 方向 / 副法线方向
        Fz 径向

        返回值点 P 和三个方向的压强Pr

        Pr.x 绕线方向压强 Fx / (channel_width * channel_depth)
        Pr.y rib 方向压强 Fy / (绕线方向长度 * channel_depth)
        Pr.z 径向压强 Fz / (绕线方向长度 * channel_width)

        关于《绕线方向长度》 == len(point_at(s + delta_length/2)-point_at(s - delta_length/2))

        返回值为 点P和 三方向压强/兆帕
        """
        p, f = self.lorentz_force_on_wire(
            s, delta_length, local_coordinate_point, other_magnet
        )

        winding_direct_length: float = (self.line3.point_at(s+delta_length/2) -
                                        self.line3.point_at(s-delta_length/2)).length()

        pressure: P3 = P3(
            x=f.x/(channel_width*channel_depth),
            y=f.y/(winding_direct_length*channel_depth),
            z=f.z/(winding_direct_length*channel_width)
        )/1e6

        return p, pressure

    @staticmethod
    def create_by_cct(cct: CCT) -> 'Wire':
        """
        由 CCT 创建 wire
        """
        def p3f(ksi):
            return cct.bipolar_toroidal_coordinate_system.convert(
                P2(ksi, cct.phi_ksi_function(ksi))
            )

        return Wire(
            line3=Line3(
                p3_function=p3f,
                start=cct.starting_point_in_ksi_phi_coordinate.x,
                end=cct.end_point_in_ksi_phi_coordinate.x
            ),
            current=cct.current,
            delta_length=cct.small_r *
            BaseUtils.angle_to_radian(360/cct.winding_number)
        )


class Function_Part:
    """
    函数段
    即包含一个函数 func，和自变量的起始值和终止值


    since 2021年1月9日
    """

    def __init__(self, func: Callable[[float], _T], start: float, end: float, scale: float = 1.0) -> None:
        self.func = func
        self.start = start
        self.end = end
        self.scale = scale

        self.length = abs(start-end) * self.scale

        # forward 为正数，说明自变量向着增大的方向，即 end > start
        self.forward = start < end

    def valve_at(self, x: float, err=1e-6) -> _T:
        """
        注意，此时函数的起点变成了 0
        取值范围为 [0, self.length]
        """
        x = x/self.scale
        if x > self.length+err or x < -err:
            print(f"Function_Part：自变量超出范围, x={x}, length={self.length}")
        return self.func(self.start + (
            x if self.forward else (-x)
        ))

    def append(self, func: Callable[[float], _T], start: float, end: float, scale: float = 1.0) -> 'Function_Part':
        """
        原函数段尾加新的函数

        这代码写得极其精妙
        """
        appended = Function_Part(func, start, end, scale)

        def fun_linked(t):
            if t < self.length:
                return self.valve_at(t)
            else:
                return appended.valve_at(t-self.length)

        return Function_Part(fun_linked, 0, self.length+appended.length)


if __name__ == "__main__":
    if True:  # test Function_Part
        fp = Function_Part(lambda x: x, 5, 2)
        print(fp.valve_at(0))
        print(fp.valve_at(1))

        fp = fp.append(lambda x: x**2, 0, 5)
        fp = fp.append(lambda x: -x**2, 0, 5, 20)
        fp = fp.append(lambda x: x**2, 0, 5)
        fp = fp.append(lambda x: -x**2, 0, 5, 20)
        Plot2.plot([P2(x, fp.valve_at(x))
                    for x in BaseUtils.linspace(0, fp.length, 100)])
        Plot2.show()
