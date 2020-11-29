"""
CCT 建模 all in one
"""
from typing import Callable, List
import matplotlib.pyplot as plt
import math
import sys
import numpy


GPU_ON = False
if GPU_ON:
    import pycuda.autoinit
    import pycuda.driver as drv
    from pycuda.compiler import SourceModule

M: float = 1.0
MM: float = 0.001
LIGHT_SPEED: float = 299792458.0 * M
RAD: float = 1.0
MRAD: float = 0.001 * RAD
J: float = 1.0
eV = 1.6021766208e-19 * J
MeV = 1000 * 1000 * eV
MeV_PER_C = 5.3442857792e-22  # kg m/s 动量单位


class BaseUtils:
    @staticmethod
    def equal(a, b, err=1e-6, msg: str = None):
        """
        判断 a b 是否相等，相等返回 true
        当 a b 不相等时，若 msg 为空，返回 flase，否则抛出异常，异常信息即 msg
        """
        if isinstance(a, float) and isinstance(b, float):
            if a == b or abs(a - b) <= err or 2 * abs((a - b) / (a + b)) <= err:
                return True
            else:
                if msg is None:
                    return False
                else:
                    raise AssertionError(msg)
        elif (isinstance(a, P2) and isinstance(b, P2)) or (
            isinstance(a, P3) and isinstance(b, P3)
        ):
            if a.__eq__(b, err=err, msg=msg):
                return True
            else:
                if msg is None:
                    return False
                else:
                    raise AssertionError(msg)
        else:
            if a == b:
                return True
            else:
                if msg is None:
                    return False
                else:
                    raise AssertionError(msg)

    @staticmethod
    def linspace(start, end, number: int) -> List:
        # 除法改成乘法以适应 P2 P3 对象
        d = (end - start) * (1/(number - 1)) 
        # i 转为浮点以适应 P2 P3 对象
        return [start + d * float(i) for i in range(number)]
        


    @staticmethod
    def angle_to_radian(deg):
        if isinstance(deg, float) or isinstance(deg, int):
            return deg / 180.0 * math.pi
        elif isinstance(deg, List):
            return [BaseUtils.angle_to_radian(d) for d in deg]
        else:
            raise NotImplementedError

    @staticmethod
    def print_traceback():
        """
        打印函数调用栈
        Returns
        -------

        """
        f = sys._getframe()
        while f is not None:
            print(f)
            f = f.f_back

    class Ellipse:
        """
        椭圆类
        Ax^2+Bxy+Cy^2=D
        """

        def __init__(self, A: float, B: float, C: float, D: float):
            self.A = A
            self.B = B
            self.C = C
            self.D = D

        def point_at(self, theta: float):
            """
            原点出发，方向th弧度的射线和椭圆Ax^2+Bxy+Cy^2=D的交点
            Parameters
            ----------
            theta 弧度

            Returns 方向th弧度的射线和椭圆Ax^2+Bxy+Cy^2=D的交点
            -------

            """
            d = P2()

            while theta < 0:
                theta += 2 * math.pi

            while theta > 2 * math.pi:
                theta -= 2 * math.pi

            if BaseUtils.equal(theta, 0) or BaseUtils.equal(theta, 2 * math.pi):
                d.x = math.sqrt(self.D / self.A)
                d.y = 0

            if BaseUtils.equal(theta, math.pi):
                d.x = -math.sqrt(self.D / self.A)
                d.y = 0

            t = 0.0

            if 0 < theta < math.pi:
                t = 1 / math.tan(theta)
                d.y = math.sqrt(self.D / (self.A * t * t + self.B * t + self.C))
                d.x = t * d.y

            if math.pi < theta < 2 * math.pi:
                theta -= math.pi
                t = 1 / math.tan(theta)
                d.y = -math.sqrt(self.D / (self.A * t * t + self.B * t + self.C))
                d.x = t * d.y

            return d

        @property
        def circumference(self) -> float:
            """
            计算椭圆周长
            Returns 计算椭圆周长
            -------

            """
            num: int = 3600 * 4
            c: float = 0.0
            for i in range(num):
                c += (
                    self.point_at(2.0 * math.pi / float(num) * (i + 1))
                    - self.point_at(2.0 * math.pi / float(num) * (i))
                ).length()

            return c

        def point_after(self, length: float):
            """
            在椭圆 Ax^2+Bxy+Cy^2=D 上行走 length，返回此时的点
            规定起点：椭圆与X轴正方向的交点
            规定行走方向：逆时针
            Parameters
            ----------
            length 行走距离

            Returns 椭圆 Ax^2+Bxy+Cy^2=D 上行走 length，返回此时的点
            -------

            """
            step_theta = BaseUtils.angle_to_radian(0.05)
            theta = 0.0
            while length > 0.0:
                length -= (
                    self.point_at(theta + step_theta) - self.point_at(theta)
                ).length()

                theta += step_theta

            return self.point_at(theta)

        def uniform_distribution_points_along_edge(self, num: int) -> List:
            points = []
            c = self.circumference
            for i in range(num):
                points.append(self.point_after(c / num * i))

            return points


class P2:
    """
    二维点 / 二维向量
    """

    def __init__(self, x: float = 0.0, y: float = 0.0):
        self.x = float(x)
        self.y = float(y)

    def length(self) -> float:
        """
        求矢量长度
        """
        return math.sqrt(self.x ** 2 + self.y ** 2)

    def normalize(self):
        """
        正则化，返回新矢量
        """
        return self * (1 / self.length())

    def change_length(self, new_length: float):
        """
        改变长度，返回新矢量
        """
        return self.normalize() * float(new_length)

    def copy(self):
        """
        拷贝
        """
        return P2(self.x, self.y)

    def __add__(self, other):
        """
        矢量加法，返回新矢量
        """
        return P2(self.x + other.x, self.y + other.y)

    def __neg__(self):
        """
        相反方向的矢量
        """
        return P2(-self.x, -self.y)

    def __sub__(self, other):
        """
        矢量减法，返回新矢量
        """
        return self.__add__(other.__neg__())

    def __iadd__(self, other):
        """
        矢量原地相加，self 自身改变
        """
        self.x += other.x
        self.y += other.y
        return self  # 必须显式返回

    def __isub__(self, other):
        """
        矢量原地减法，self 自身改变
        """
        self.x -= other.x
        self.y -= other.y
        return self

    def __matmul(self, m: List[List[float]]):
        """
        2*2矩阵和 self 相乘，仅仅用于矢量旋转
        """
        return P2(
            m[0][0] * self.x + m[0][1] * self.y, m[1][0] * self.x + m[1][1] * self.y
        )

    @staticmethod
    def __rotate_r(phi: float) -> List[List[float]]:
        """
        旋转矩阵
        """
        return [[math.cos(phi), -math.sin(phi)], [math.sin(phi), math.cos(phi)]]

    def rotate(self, phi: float):
        """
        矢量自身旋转
        """
        return self.__matmul(P2.__rotate_r(phi))

    def angle_to_x_axis(self) -> float:
        """
        矢量和 x 轴的夹角，弧度
        """
        a = float(math.atan2(self.y, self.x))
        return a if a >= 0 else math.pi * 2 + a

    def __mul__(self, other):
        """
        矢量乘标量，各元素相乘，返回新矢量
        矢量乘矢量，内积，返回标量
        """
        if isinstance(other, float):
            return P2(self.x * other, self.y * other)
        else:
            return self.x * other.x + self.y * other.y

    def angle_to(self, other) -> float:
        """
        矢量 self 和 另一个矢量 other 的夹角
        """
        theta = (self * other) / (self.length() * other.length())
        return float(math.acos(theta))

    def to_p3(self, transformation: Callable = lambda p2: P3(p2.x, p2.y, 0.0)):
        """
        二维矢量转为三维
        默认情况返回 [x,y,0]
        """
        return transformation(self)

    def __str__(self):
        """
        用于打印矢量值
        """
        return f"P2:x({self.x})y({self.y})"

    def __eq__(self, other, err=1e-6, msg=None):
        """
        矢量相等判断
        """
        return BaseUtils.equal(self.x, other.x, err, msg) and BaseUtils.equal(
            self.y, other.y, err, msg
        )

    @staticmethod
    def x_direct(x: float = 1.0):
        return P2(x=x)

    @staticmethod
    def y_direct(y: float = 1.0):
        return P2(y=y)

    @staticmethod
    def origin():
        return P2()

    @staticmethod
    def zeros():
        return P2()

    def to_list(self) -> List[float]:
        return [self.x, self.y]


class P3:
    """
    三维点 / 三维矢量
    """

    def __init__(self, x: float = 0.0, y: float = 0.0, z: float = 0.0):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def length(self) -> float:
        """
        矢量长度
        """
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def normalize(self):
        """
        正则化，返回新矢量
        """
        return self * (1 / self.length())

    def change_length(self, new_length: float):
        """
        改变长度，返回新矢量
        """
        return self.normalize() * new_length

    def copy(self):
        """
        拷贝
        """
        return P3(self.x, self.y, self.z)

    def __add__(self, other):
        """
        矢量相加
        """
        return P3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __neg__(self):
        """
        相反矢量
        """
        return P3(-self.x, -self.y, -self.z)

    def __sub__(self, other):
        """
        矢量相减
        """
        return self.__add__(other.__neg__())

    def __iadd__(self, other):
        """
        矢量原地相加
        """
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __isub__(self, other):
        """
        矢量原地减法
        """
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        return self

    def __mul__(self, other):
        """
        矢量乘标量，各元素相乘，返回新矢量
        矢量乘矢量，内积，返回标量
        """
        if isinstance(other, float):
            return P3(self.x * other, self.y * other, self.z * other)
        else:
            return self.x * other.x + self.y * other.y + self.z * other.z

    def __matmul__(self, other):
        """
        矢量叉乘 / 外积，返回新矢量
        """
        return P3(
            self.y * other.z - self.z * other.y,
            -self.x * other.z + self.z * other.x,
            self.x * other.y - self.y * other.x,
        )

    def __str__(self):
        """
        矢量信息
        """
        return f"P3:x({self.x})y({self.y})z({self.z})"

    def __eq__(self, other, err=1e-6, msg=None):
        """
        矢量相等判断
        """
        return (
            BaseUtils.equal(self.x, other.x, err, msg)
            and BaseUtils.equal(self.y, other.y, err, msg)
            and BaseUtils.equal(self.z, other.z, err, msg)
        )

    @staticmethod
    def x_direct(x: float = 1.0):
        return P3(x=x)

    @staticmethod
    def y_direct(y: float = 1.0):
        return P3(y=y)

    @staticmethod
    def z_direct(z: float = 1.0):
        return P3(z=z)

    @staticmethod
    def origin():
        return P3()

    @staticmethod
    def zeros():
        return P3()

    def to_list(self) -> List[float]:
        return [self.x, self.y, self.z]

    @staticmethod
    def from_numpy_ndarry3(array3: numpy.ndarray):
        return P3(array3[0], array3[1], array3[2])


class Magnet:
    """
    表示一个可以求磁场的对象，如 CCT 、 QS 磁铁
    所有实现此接口的类，可以计算出它在某一点的磁场

    本类（接口）只有一个接口方法 magnetic_field_at
    """

    def magnetic_field_at(self, point: P3) -> P3:
        """
        获得本对象 self 在点 point 处产生的磁场
        这个方法需要在子类中实现/重写
        ----------
        point 三维笛卡尔坐标系中的点，即一个三维矢量，如 [0,0,0]

        Returns 本对象 self 在点 point 处的磁场，用三维矢量表示
        -------
        """
        raise NotImplementedError

    def magnetic_field_along(self, line) -> List[P3]:
        """
        计算本对象在三维曲线 line 上的磁场分布
        ----------
        line 由离散点组成的三维曲线

        Returns 本对象在三维曲线 line 上的磁场分布，用三维矢量的数组表示
        -------
        """
        if isinstance(line, List):
            return [self.magnetic_field_at(p) for p in line]
        elif isinstance(line, Line2):
            return self.magnetic_field_along(line.disperse3d())
        else:
            raise NotImplementedError


class LocalCoordinateSystem:
    """
    局部坐标系，各种磁铁都放置在局部坐标系中，而粒子在全局坐标系中运动，为了求磁铁在粒子位置产生的磁场，需要引入坐标变换
    """

    def __init__(
        self,
        location: P3 = P3.origin(),
        x_direction: P3 = P3.x_direct(),
        z_direction: P3 = P3.z_direct(),
    ):
        """
        Parameters
        ----------
        location 全局坐标系中实体位置，默认全局坐标系的远点
        x_direction 局部坐标系 x 方向，默认全局坐标系 x 方向
        z_direction 局部坐标系 z 方向，默认全局坐标系 z 方向
        """
        BaseUtils.equal(
            x_direction * z_direction,
            0.0,
            msg=f"创建 LocalCoordinateSystem 对象异常，x_direction{x_direction}和z_direction{z_direction}不正交",
        )

        # 局部坐标系，原心
        self.location = location.copy()

        # 局部坐标系的 x y z 三方向
        self.ZI = z_direction.copy().normalize()
        self.XI = x_direction.copy().normalize()
        self.YI = self.ZI @ self.XI

    def point_to_local_coordinate(self, global_coordinate_point: P3) -> P3:
        """
        全局坐标系 -> 局部坐标系
        Parameters
        ----------
        global_coordinate_point 全局坐标系中的点

        Returns 这一点在局部坐标系中的坐标
        -------
        """
        location_to_global_coordinate = global_coordinate_point - self.location
        x = self.XI * location_to_global_coordinate
        y = self.YI * location_to_global_coordinate
        z = self.ZI * location_to_global_coordinate
        return P3(x, y, z)

    def point_to_global_coordinate(self, local_coordinate_point: P3) -> P3:
        """
        局部坐标系 -> 全局坐标系
        Parameters
        ----------
        local_coordinate_point 局部坐标系

        Returns 全局坐标系
        -------

        """

        return self.location + (
            local_coordinate_point.x * self.XI
            + local_coordinate_point.y * self.YI
            + local_coordinate_point.z * self.ZI
        )

    def line_to_local_coordinate(self, global_coordinate_line: List[P3]) -> List[P3]:
        """
        全局坐标系中的 线/点集 坐标转为局部坐标系
        线/点集，为 n*3 的矩阵，矩阵每一行代表一个点
        Parameters
        ----------
        global_coordinate_line 全局坐标系中的 线/点集

        Returns 转为局部坐标系
        -------

        """
        return [self.point_to_local_coordinate(p) for p in global_coordinate_line]

    def line_to_global_coordinate(self, local_coordinate_line: List[P3]) -> List[P3]:
        """
        局部坐标系中的 线/点集 坐标转为全局坐标系
        线/点集，为 n*3 的矩阵，矩阵每一行代表一个点
        Parameters
        ----------
        local_coordinate_line 局部坐标系中的 线/点集

        Returns 转为全局坐标系
        -------

        """
        return [self.point_to_global_coordinate(p) for p in local_coordinate_line]

    def set_direction(self, z_direction: P3, x_direction: P3) -> None:
        BaseUtils.equal(
            z_direction * x_direction,
            0.0,
            msg="LocalCoordinateSystem:set_direction异常，"
            + f"z_direction{z_direction}x_direction{x_direction}不正交",
        )
        self.ZI = z_direction.copy().normalize()
        self.XI = x_direction.copy().normalize()
        self.YI = self.ZI @ self.XI

    def __str__(self) -> str:
        return f"ORIGIN={self.location}, xi={self.XI}, yi={self.YI}, zi={self.ZI}"

    @staticmethod
    def create_by_y_and_z_direction(location: P3, y_direction: P3, z_direction: P3):
        """
        由 原点 location y方向 y_direction 和 z方向 z_direction 创建坐标系
        Parameters
        ----------
        location 原点
        y_direction y方向
        z_direction z方向

        Returns 坐标系
        -------

        """
        BaseUtils.equal(
            y_direction * z_direction,
            0.0,
            msg=f"创建 LocalCoordinateSystem 对象异常，y_direction{y_direction}和z_direction{z_direction}不正交",
        )

        x_direction = y_direction @ z_direction
        return LocalCoordinateSystem(location, z_direction, x_direction)

    @staticmethod
    def global_coordinate_system():
        """
        获取全局坐标系
        Returns 全局坐标系
        -------

        """
        return LocalCoordinateSystem()


class Line2:
    """
    二维 xy 平面的一条有起点和终点的连续曲线段，可以是直线、圆弧
    本类包含 3 个抽象方法，需要实现：
    get_length 获得曲线长度
    point_at 从曲线起点出发，s 位置处的点
    direct_at 从曲线起点出发，s 位置处曲线方向

    说明：这个类主要用于构建 “理想轨道”，理想轨道的用处很多：
    1. 获取理想轨道上的理想粒子；
    2. 研究理想轨道上的磁场分布

    """

    def get_length(self) -> float:
        """
        获得曲线的长度
        Returns 曲线的长度
        -------

        """
        raise NotImplementedError

    def point_at(self, s: float) -> P2:
        """
        获得曲线 s 位置处的点 (x,y)
        即从曲线起点出发，运动 s 长度后的位置
        Parameters
        ----------
        s 长度量，曲线上 s 位置

        Returns 曲线 s 位置处的点 (x,y)
        -------

        """
        raise NotImplementedError

    def direct_at(self, s: float) -> P2:
        """
        获得 s 位置处，曲线的方向
        Parameters
        ----------
        s 长度量，曲线上 s 位置

        Returns s 位置处，曲线的方向
        -------

        """
        raise NotImplementedError

    def right_hand_side_point(self, s: float, d: float) -> P2:
        """
        位于 s 处的点，它右手边 d 处的点

         1    5    10     15
         -----------------@------>
         |2
         |4               *
        如上图，一条直线，s=15，d=4 ,即点 @ 右手边 4 距离处的点 *

        说明：这个方法，主要用于四极场、六极场的计算，因为需要涉及轨道横向位置的磁场

        Parameters
        ----------
        s 长度量，曲线上 s 位置
        d 长度量，d 距离远处

        Returns 位于 s 处的点，它右手边 d 处的点
        -------

        """
        ps = self.point_at(s)

        # 方向
        ds = self.direct_at(s)

        return ps + ds.copy().rotate(-math.pi / 2).change_length(d)

    def left_hand_side_point(self, s: float, d: float) -> P2:
        """
        位于 s 处的点，它左手边 d 处的点
        说明见 right_hand_side_point 方法
        Parameters
        ----------
        s 长度量，曲线上 s 位置
        d 长度量，d 距离远处

        Returns 位于 s 处的点，它左手边 d 处的点
        -------

        """
        return self.right_hand_side_point(s, -d)

    # ------------------------------端点性质-------------------- #
    def point_at_start(self):
        return self.point_at(0.0)

    def point_at_end(self):
        return self.point_at(self.get_length())

    def direct_at_start(self):
        return self.direct_at(0.0)

    def direct_at_end(self):
        return self.direct_at(self.get_length())

    # ------------------------------平移-------------------- #
    def __add__(self, v2: P2):
        """
        Line2 的平移， v2 表示移动的方向和距离
        Parameters
        ----------
        v2 二维向量

        Returns 平移后的 Line2
        -------

        """

        class MovedLine2(Line2):
            def __init__(self, hold):
                self.hold = hold

            def get_length(self) -> float:
                return self.hold.get_length()

            def point_at(self, s: float) -> P2:
                return self.hold.point_at(s) + v2

            def direct_at(self, s: float) -> P2:
                return self.hold.direct_at(s)

        return MovedLine2(self)

    # ------------------------------ 离散 ------------------------#
    def disperse2d(self, step: float = 1.0 * MM) -> List[P2]:
        """
        二维离散轨迹点
        Parameters
        ----------
        step 步长

        Returns 二维离散轨迹点
        -------

        """
        number: int = int(math.ceil(self.get_length() / step))
        return [
            self.point_at(s) for s in BaseUtils.linspace(0, self.get_length(), number)
        ]

    def disperse3d(self, step: float = 1.0 * MM) -> List[P3]:
        """
        三维离散轨迹点，其中第三维 z == 0.0
        Parameters
        ----------
        step 步长

        Returns 二维离散轨迹点
        -------

        """
        return [p.to_p3() for p in self.disperse2d(step)]

    def __str__(self):
        return f"Line2，起点{self.point_at_start()}，长度{self.get_length()}"


class StraightLine2(Line2):
    """
    二维直线段，包含三个参数：长度、方向、起点
    """

    def __init__(self, length: float, direct: P2, start_point: P2):
        self.length = length
        self.direct = direct
        self.start_point = start_point

    def get_length(self) -> float:
        return self.length

    def point_at(self, s: float) -> P2:
        return self.start_point + self.direct.copy().change_length(s)

    def direct_at(self, s: float) -> P2:
        return self.direct

    def __str__(self):
        return f"直线段，起点{self.start_point}，方向{self.direct}，长度{self.length}"


class ArcLine2(Line2):
    """
    二维圆弧段
    借助极坐标的思想来描述圆弧
    基础属性： 圆弧的半径 radius、圆弧的圆心 center
    起点描述：极坐标 phi 值
    弧长：len = radius * totalPhi

    起点start_point、圆心center、半径radius、旋转方向clockwise、角度totalPhi 五个自由度
    起点弧度值 starting_phi、起点处方向、半径radius、旋转方向clockwise、角度totalPhi 五个自由度

    如图： *1 表示起点方向，@ 是圆心，上箭头 ↑ 是起点处方向，旋转方向是顺时针，*5 是终点，因此角度大约是 80 deg
                *5
           *4
       *3
     *2
    *1     ↑       @

    """

    def __init__(
        self,
        starting_phi: float,
        center: P2,
        radius: float,
        total_phi: float,
        clockwise: bool,
    ):
        self.starting_phi = starting_phi
        self.center = center
        self.radius = radius
        self.total_phi = total_phi
        self.clockwise = clockwise
        self.length = radius * total_phi

    def get_length(self) -> float:
        return self.length

    def point_at(self, s: float) -> P2:
        phi = s / self.radius
        current_phi = (
            self.starting_phi - phi if self.clockwise else self.starting_phi + phi
        )

        uc = ArcLine2.unit_circle(current_phi)

        return uc.change_length(self.radius) + self.center

    def direct_at(self, s: float) -> P2:
        phi = s / self.radius
        current_phi = (
            self.starting_phi - phi if self.clockwise else self.starting_phi + phi
        )

        uc = ArcLine2.unit_circle(current_phi)

        return uc.rotate(-math.pi / 2 if self.clockwise else math.pi / 2)

    @staticmethod
    def create(
        start_point: P2,
        start_direct: P2,
        radius: float,
        clockwise: bool,
        total_deg: float,
    ):
        center: P2 = start_point + start_direct.copy().rotate(
            -math.pi / 2 if clockwise else math.pi / 2
        ).change_length(radius)

        starting_phi = (start_point - center).angle_to_x_axis()

        total_phi = BaseUtils.angle_to_radian(total_deg)

        return ArcLine2(starting_phi, center, radius, total_phi, clockwise)

    @staticmethod
    def unit_circle(phi: float) -> P2:
        """
        单位圆（极坐标）
        返回：极坐标(r=1.0,phi=phi)的点的直角坐标(x,y)
        Parameters
        ----------
        phi 极坐标phi

        Returns 单位圆上的一点
        -------

        """
        x = math.cos(phi)
        y = math.sin(phi)

        return P2(x, y)

    def __str__(self):
        return (
            f"弧线段，起点{self.point_at_start()}，"
            + f"方向{self.direct_at_start()}，顺时针{self.clockwise}，半径{self.radius}，角度{self.total_phi}"
        )


class Trajectory(Line2):
    """
    二维轨迹，由直线+圆弧组成
    """

    def __init__(self, first_line2: Line2):
        """
        构造器，传入第一条线 first_line2
        Parameters
        ----------
        first_line2 第一条线
        -------

        """
        self.__trajectoryList = [first_line2]
        self.__length = first_line2.get_length()
        self.__point_at_error_happen = False  # 是否发生 point_at 错误

    def add_strait_line(self, length: float):
        """
        尾接直线
        Parameters
        ----------
        length 直线长度

        Returns self
        -------

        """
        last_line = self.__trajectoryList[-1]
        sp = last_line.point_at_end()
        sd = last_line.direct_at_end()

        sl = StraightLine2(length, sd, sp)

        self.__trajectoryList.append(sl)
        self.__length += length

        return self

    def add_arc_line(self, radius: float, clockwise: bool, angle_deg: float):
        """
        尾接圆弧
        Parameters
        ----------
        radius 半径
        clockwise 顺时针？
        angle_deg 角度

        Returns self
        -------

        """
        last_line = self.__trajectoryList[-1]
        sp = last_line.point_at_end()
        sd = last_line.direct_at_end()

        al = ArcLine2.create(sp, sd, radius, clockwise, angle_deg)

        self.__trajectoryList.append(al)
        self.__length += al.get_length()

        return self

    def get_length(self) -> float:
        return self.__length

    def point_at(self, s: float) -> P2:
        s0 = s

        for line in self.__trajectoryList:
            if s <= line.get_length():
                return line.point_at(s)
            else:
                s -= line.get_length()

        last_line = self.__trajectoryList[-1]

        # 2020年4月2日
        # 解决了一个因为浮点数产生的巨大bug
        if abs(s) <= 1e-8:
            return last_line.point_at_end()

        if not self.__point_at_error_happen:
            self.__point_at_error_happen = True
            print(f"ERROR Trajectory::point_at{s0}")
            BaseUtils.print_traceback()

        return last_line.point_at(s)

    def direct_at(self, s: float) -> P2:
        s0 = s

        for line in self.__trajectoryList:
            if s <= line.get_length():
                return line.direct_at(s)
            else:
                s -= line.get_length()

        last_line = self.__trajectoryList[-1]

        # 2020年4月2日
        # 解决了一个因为浮点数产生的巨大bug
        if abs(s) <= 1e-8:
            return last_line.direct_at_end()

        if not self.__point_at_error_happen:
            self.__point_at_error_happen = True
            print(f"ERROR Trajectory::direct_at{s0}")
            BaseUtils.print_traceback()

        return last_line.direct_at(s)

    def __str__(self):
        details = "\t\n".join(self.__trajectoryList.__str__())
        return f"Trajectory:{details}"


class Protons:
    """
    质子相关常量和计算
    """

    # 静止质量
    STATIC_MASS_KG = 1.672621898e-27

    # 静止能量 = m0 * c ^ 2 单位焦耳
    STATIC_ENERGY_J = STATIC_MASS_KG * LIGHT_SPEED * LIGHT_SPEED

    # 静止能量 eV 为单位
    STATIC_ENERGY_eV = STATIC_ENERGY_J / eV

    # 静止能量 MeV 为单位，应该是 STATIC_ENERGY_J / MeV。但是写成字面量
    STATIC_ENERGY_MeV = 938.2720813

    # 电荷量 库伦
    CHARGE_QUANTITY = 1.6021766208e-19

    @classmethod
    def get_total_energy_MeV(cls, kinetic_energy_MeV: float) -> float:
        """
        计算总能量 MeV = 静止能量 + 动能
        Parameters
        ----------
        kinetic_energy_MeV 动能 MeV 一般为 250 Mev

        Returns 总能量 MeV
        -------

        """
        return cls.STATIC_ENERGY_MeV + kinetic_energy_MeV

    @classmethod
    def get_total_energy_J(cls, kinetic_energy_MeV: float) -> float:
        """
        计算总能量 焦耳
        Parameters
        ----------
        kinetic_energy_MeV 动能 MeV 一般为 250 Mev

        Returns 总能量 焦耳
        -------

        """
        return cls.get_total_energy_MeV(kinetic_energy_MeV) * MeV

    @classmethod
    def get_relativistic_mass(cls, kinetic_energy_MeV: float) -> float:
        """
        计算动质量 kg = 动能 / (c^2)
        Parameters
        ----------
        kinetic_energy_MeV 动能 MeV 一般为 250 Mev

        Returns 动质量 kg
        -------

        """
        return cls.get_total_energy_J(kinetic_energy_MeV) / LIGHT_SPEED / LIGHT_SPEED

    @classmethod
    def get_speed_m_per_s(cls, kinetic_energy_MeV: float) -> float:
        """
        计算速度 m/s = c * sqrt( 1 - (m0/m)^2 )
        Parameters
        ----------
        kinetic_energy_MeV 动能 MeV 一般为 250 Mev

        Returns 速度 m/s
        -------

        """
        return LIGHT_SPEED * math.sqrt(
            1
            - (cls.STATIC_MASS_KG / cls.get_relativistic_mass(kinetic_energy_MeV)) ** 2
        )

    @classmethod
    def get_momentum_kg_m_pre_s(cls, kinetic_energy_MeV: float) -> float:
        """
        动量 kg m/s
        Parameters
        ----------
        kinetic_energy_MeV 动能 MeV 一般为 250 Mev

        Returns 动量 kg m/s
        -------

        """
        return cls.get_relativistic_mass(kinetic_energy_MeV) * cls.get_speed_m_per_s(
            kinetic_energy_MeV
        )

    @classmethod
    def getMomentum_MeV_pre_c(cls, kinetic_energy_MeV: float) -> float:
        """
        动量 MeV/c
        Parameters 动能 MeV 一般为 250 Mev
        ----------
        kinetic_energy_MeV

        Returns 动量 MeV/c
        -------

        """
        return cls.get_momentum_kg_m_pre_s(kinetic_energy_MeV) / MeV_PER_C

    @classmethod
    def get_magnetic_stiffness(cls, kinetic_energy_MeV: float) -> float:
        """
        磁钢度 T/m
        Parameters
        ----------
        kinetic_energy_MeV 动能 MeV 一般为 250 Mev

        Returns 磁钢度 T/m
        -------

        """
        return cls.get_momentum_kg_m_pre_s(kinetic_energy_MeV) / cls.CHARGE_QUANTITY

    # ------------------  动量分散相关  ----------------------
    @classmethod
    def get_kinetic_energy_MeV(cls, momentum_KG_M_PER_S: float) -> float:
        # 速度
        speed = momentum_KG_M_PER_S / math.sqrt(
            cls.STATIC_MASS_KG ** 2 + (momentum_KG_M_PER_S / LIGHT_SPEED) ** 2
        )
        # 动质量
        relativistic_mass = cls.STATIC_MASS_KG / math.sqrt(
            1 - (speed / LIGHT_SPEED) ** 2
        )
        # 总能量 J
        total_energy_J = relativistic_mass * LIGHT_SPEED * LIGHT_SPEED
        # 动能 J
        k = total_energy_J - cls.STATIC_ENERGY_J

        return k / MeV

    @classmethod
    def get动量分散后的动能(cls, 原动能_MeV: float, 动量分散: float):
        """
        英文版见下
        Parameters
        ----------
        原动能_MeV
        动量分散

        Returns 动量分散后的动能 MeV
        -------

        """
        原动量 = cls.get_momentum_kg_m_pre_s(原动能_MeV)

        新动量 = 原动量 * (1 + 动量分散)

        新动能 = cls.get_kinetic_energy_MeV(新动量)

        return 新动能

    @classmethod
    def get_kinetic_energy_MeV_after_momentum_dispersion(
        cls, old_kinetic_energy_MeV: float, momentum_dispersion: float
    ) -> float:
        """
        中文版见上
        Parameters
        ----------
        old_kinetic_energy_MeV 原动能_MeV
        momentum_dispersion 动量分散

        Returns 动量分散后的动能 MeV
        -------

        """
        momentum0 = cls.get_momentum_kg_m_pre_s(old_kinetic_energy_MeV)

        momentum = momentum0 * (1 + momentum_dispersion)

        kinetic_energy = cls.get_kinetic_energy_MeV(momentum)

        return kinetic_energy

    @classmethod
    def convert动量分散_TO_能量分散(cls, 动量分散: float, 动能_MeV: float) -> float:
        """
        下方法的中文版
        Parameters
        ----------
        动量分散
        动能_MeV

        Returns convert动量分散_TO_能量分散
        -------

        """
        k = (动能_MeV + cls.STATIC_ENERGY_MeV) / (动能_MeV + 2 * cls.STATIC_ENERGY_MeV)

        return 动量分散 / k

    @classmethod
    def convert_momentum_dispersion_to_energy_dispersion(
        cls, momentum_dispersion: float, kinetic_energy_MeV: float
    ) -> float:
        """
        上方法的英文版
        Parameters
        ----------
        momentum_dispersion 动量分散
        kinetic_energy_MeV 动能_MeV

        Returns convert动量分散_TO_能量分散
        -------

        """
        k = (kinetic_energy_MeV + cls.STATIC_ENERGY_MeV) / (
            kinetic_energy_MeV + 2 * cls.STATIC_ENERGY_MeV
        )

        return momentum_dispersion / k

    @classmethod
    def convert能量分散_TO_动量分散(cls, 能量分散: float, 动能_MeV: float) -> float:
        k = (动能_MeV + cls.STATIC_ENERGY_MeV) / (动能_MeV + 2 * cls.STATIC_ENERGY_MeV)
        return 能量分散 * k

    @classmethod
    def convert_energy_dispersion_to_momentum_dispersion(
        cls, energyDispersion: float, kineticEnergy_MeV: float
    ) -> float:
        """
        上方法的英文版
        Parameters
        ----------
        energyDispersion 能量分散
        kineticEnergy_MeV 动能，典型值 250

        Returns 动量分散
        -------

        """
        k = (kineticEnergy_MeV + cls.STATIC_ENERGY_MeV) / (
            kineticEnergy_MeV + 2 * cls.STATIC_ENERGY_MeV
        )
        return energyDispersion * k


class RunningParticle:
    """
    在全局坐标系中运动的一个粒子
    position 位置，三维矢量，单位 [m, m, m]
    velocity 速度，三位矢量，单位 [m/s, m/s, m/s]
    relativistic_mass 相对论质量，又称为动质量，单位 kg， M=Mo/√(1-v^2/c^2)
    e 电荷量，单位 C 库伦
    speed 速率，单位 m/s
    distance 运动距离，单位 m
    """

    def __init__(
        self,
        position: P3,
        velocity: P3,
        relativistic_mass: float,
        e: float,
        speed: float,
        distance: float = 0.0,
    ):
        """
        在全局坐标系中运动的一个粒子
        Parameters
        ----------
        position 位置，三维矢量，单位 [m, m, m]
        velocity 速度，三位矢量，单位 [m/s, m/s, m/s]
        relativistic_mass 相对论质量，又称为动质量，单位 kg， M=Mo/√(1-v^2/c^2)
        e 电荷量，单位 C 库伦
        speed 速率，单位 m/s
        distance 运动距离，单位 m
        """
        self.position = position
        self.velocity = velocity
        self.relativistic_mass = relativistic_mass
        self.e = e
        self.speed = speed
        self.distance = distance

    def run_self_in_magnetic_field(
        self, magnetic_field: P3, footstep: float = 1 * MM
    ) -> None:
        """
        粒子在磁场 magnetic_field 中运动 footstep 长度
        Parameters
        ----------
        magnetic_field 磁场，看作恒定场
        footstep 步长，默认 1 MM

        Returns None
        -------
        """
        # 计算受力 qvb
        f = (self.velocity @ magnetic_field) * self.e
        # 计算加速度 a = f/m
        a = f / self.relativistic_mass
        # 计算运动时间
        t = footstep / self.speed
        # 位置变化
        self.position += self.velocity * t
        # 速度变化
        self.velocity += t * a
        # 运动距离
        self.distance += footstep

    def copy(self):
        """
        深拷贝粒子
        Returns 深拷贝粒子
        -------

        """
        return RunningParticle(
            self.position.copy(),
            self.velocity.copy(),
            self.relativistic_mass,
            self.e,
            self.speed,
            self.distance,
        )

    def compute_scalar_momentum(self) -> float:
        """
        获得标量动量
        Returns 标量动量
        -------

        """
        return self.speed * self.relativistic_mass

    def change_scalar_momentum(self, scalar_momentum: float) -> None:
        """
        改变粒子的标量动量。
        注意：真正改变的是粒子的速度和动质量
        这个方法用于生成一组动量分散的粒子

        scalar_momentum 标量动量
        Returns None
        -------

        """
        # 先求 静止质量
        m0 = self.relativistic_mass * math.sqrt(
            1 - (self.speed ** 2) / (LIGHT_SPEED ** 2)
        )
        # 求新的速率
        new_speed = scalar_momentum / math.sqrt(
            m0 ** 2 + (scalar_momentum / LIGHT_SPEED) ** 2
        )
        # 求新的动质量
        new_relativistic_mass = m0 / math.sqrt(1 - (new_speed / LIGHT_SPEED) ** 2)
        # 求新的速度
        new_velocity: P3 = self.velocity.change_length(new_speed)

        # 写入
        self.relativistic_mass = new_relativistic_mass
        self.speed = new_speed
        self.velocity = new_velocity

        # 验证
        BaseUtils.equal(
            scalar_momentum,
            self.compute_scalar_momentum(),
            msg=f"RunningParticle::change_scalar_momentum异常，scalar_momentum{scalar_momentum}!=self.compute_scalar_momentum{self.compute_scalar_momentum}",
            err=1e-6,
        )

        BaseUtils.equal(
            self.speed,
            self.velocity.length(),
            msg=f"RunningParticle::change_scalar_momentum异常,self.speed{self.speed}!=Vectors.length(self.velocity){self.velocity.length()}",
        )

    def get_natural_coordinate_system(
        self, y_direction: P3 = P3.z_direct()
    ) -> LocalCoordinateSystem:
        return LocalCoordinateSystem.create_by_y_and_z_direction(
            self.position, y_direction, self.velocity
        )

    def __str__(self) -> str:
        return f"p={self.position},v={self.velocity},v0={self.speed}"


class ParticleRunner:
    """
    粒子运动工具类
    """

    @staticmethod
    def run_only(
        p: RunningParticle, m: Magnet, length: float, footstep: float = 1 * MM
    ) -> None:
        """
        让粒子 p 在磁场 m 中运动 length 距离，步长 footstep
        Parameters
        ----------
        p 粒子
        m 磁场
        length 运动长度
        footstep 步长

        Returns None
        -------

        """
        distance = 0.0
        while distance < length:
            p.run_self_in_magnetic_field(
                m.magnetic_field_at(p.position), footstep=footstep
            )
            distance += footstep

    @staticmethod
    def run_get_trajectory(
        p: RunningParticle, m: Magnet, length: float, footstep: float = 1 * MM
    ) -> List[P3]:
        """
        让粒子 p 在磁场 m 中运动 length 距离，步长 footstep
        获得粒子的轨迹
        Parameters
        ----------
        p 粒子
        m 磁场
        length 运动长度
        footstep 步长

        Returns 轨迹 np.ndarray，是三维点的数组
        -------

        """
        trajectory: List[P3] = [p.position.copy()]

        i = 1
        distance = 0.0
        while distance < length:
            p.run_self_in_magnetic_field(
                m.magnetic_field_at(p.position), footstep=footstep
            )
            distance += footstep
            trajectory.append(p.position.copy())
            i += 1

        return trajectory

    @staticmethod
    def run_get_all_info(
        p: RunningParticle, m: Magnet, length: float, footstep: float = 1 * MM
    ) -> List[RunningParticle]:
        """
        让粒子 p 在磁场 m 中运动 length 距离，步长 footstep
        获得粒子全部信息
        Parameters
        ----------
        p 粒子
        m 磁场
        length 运动长度
        footstep 步长

        Returns 每一步处的粒子全部信息 List[RunningParticle]
        -------

        """
        all_info: List[RunningParticle] = [p.copy()]
        distance = 0.0
        while distance < length:
            p.run_self_in_magnetic_field(
                m.magnetic_field_at(p.position), footstep=footstep
            )
            distance += footstep
            all_info.append(p.copy())

        return all_info

    @staticmethod
    def run_ps_only_cpu0(
        ps: List[RunningParticle], m: Magnet, length: float, footstep: float = 1 * MM
    ) -> None:
        """
        让粒子群 ps 在磁场 m 中运动 length 距离，步长 footstep
        CPU 计算 单线程
        Parameters
        ----------
        ps 一群粒子
        m 磁场
        length 运动长度
        footstep 步长


        Returns None
        -------

        """
        for p in ps:
            ParticleRunner.run_only(p, m, length, footstep)


class PhaseSpaceParticle:
    XXP_PLANE = 1
    YYP_PLANE = 2

    """
    相空间中的粒子，6个坐标 x xp y yp z delta
    """

    def __init__(
        self, x: float, xp: float, y: float, yp: float, z: float, delta: float
    ):
        self.x = x
        self.xp = xp
        self.y = y
        self.yp = yp
        self.z = z
        self.delta = delta

    def project_to_xxp_plane(self) -> P2:
        """
        投影到 x-xp 平面
        Returns [self.x, self.xp]
        -------

        """
        return P2(self.x, self.xp)

    def project_to_yyp_plane(self) -> P2:
        """
        投影到 y-yp 平面
        Returns [self.y, self.yp]
        -------

        """
        return P2(self.y, self.yp)

    def project_to_plane(self, plane_id: int) -> P2:
        if plane_id == PhaseSpaceParticle.XXP_PLANE:
            return self.project_to_xxp_plane()
        elif plane_id == PhaseSpaceParticle.YYP_PLANE:
            return self.project_to_yyp_plane()
        else:
            raise ValueError(f"没有处理plane_id({plane_id})的方法")

    @staticmethod
    def phase_space_particles_along_positive_ellipse_in_xxp_plane(
        xMax: float, xpMax: float, delta: float, number: int
    ) -> List:
        """
        获取分布于 x xp 平面上 正相椭圆上的 PhaseSpaceParticles
        注意是 正相椭圆
        Parameters
        ----------
        xMax 相椭圆参数 x 最大值
        xpMax 相椭圆参数 xp 最大值
        delta 动量分散
        number 粒子数目

        Returns 分布于 x xp 平面上 正相椭圆上的 PhaseSpaceParticles
        -------

        """
        A: float = 1 / (xMax ** 2)
        B: float = 0
        C: float = 1 / (xpMax ** 2)
        D: float = 1

        return [
            PhaseSpaceParticle(p[0], p[1], 0, 0, 0, delta)
            for p in BaseUtils.Ellipse(
                A, B, C, D
            ).uniform_distribution_points_along_edge(number)
        ]

    @staticmethod
    def phase_space_particles_along_positive_ellipse_in_yyp_plane(
        yMax: float, ypMax: float, delta: float, number: int
    ) -> List:
        """
        获取分布于 y yp 平面上 正相椭圆上的 PhaseSpaceParticles
        注意是 正相椭圆
        Parameters
        ----------
        yMax 相椭圆参数 y 最大值
        ypMax 相椭圆参数 yp 最大值
        delta 动量分散
        number 粒子数目

        Returns 分布于 y yp 平面上 正相椭圆上的 PhaseSpaceParticles
        -------

        """
        A: float = 1 / (yMax ** 2)
        B: float = 0
        C: float = 1 / (ypMax ** 2)
        D: float = 1

        return [
            PhaseSpaceParticle(0, 0, p[0], p[1], 0, delta)
            for p in BaseUtils.Ellipse(
                A, B, C, D
            ).uniform_distribution_points_along_edge(number)
        ]

    @staticmethod
    def phase_space_particles_along_positive_ellipse_in_plane(
        plane_id: int, xMax: float, xpMax: float, delta: float, number: int
    ) -> List:
        """
        获取分布于 x xp 平面上或 y yp 平面上的，正相椭圆上的 PhaseSpaceParticles
        Parameters
        ----------
        xxPlane x 平面或 y 平面，true：x 平面，false:y 平面
        xMax 相椭圆参数 x/y 最大值
        xpMax 相椭圆参数 xp/yp 最大值
        delta 动量分散
        number 粒子数目

        Returns 分布于 x xp 平面上或 y yp 平面上的，正相椭圆上的 PhaseSpaceParticles
        -------

        """
        if plane_id == PhaseSpaceParticle.XXP_PLANE:
            return PhaseSpaceParticle.phase_space_particles_along_positive_ellipse_in_xxp_plane(
                xMax, xpMax, delta, number
            )
        elif plane_id == PhaseSpaceParticle.YYP_PLANE:
            return PhaseSpaceParticle.phase_space_particles_along_positive_ellipse_in_yyp_plane(
                xMax, xpMax, delta, number
            )
        else:
            raise ValueError(f"没有处理plane_id({plane_id})的方法")

    @staticmethod
    def phase_space_particles_project_to_xxp_plane(
        phase_space_particles: List,
    ) -> List[P2]:
        """
        相空间粒子群投影到 x 平面
        Parameters
        ----------
        phase_space_particles 相空间粒子群

        Returns 相空间粒子群投影到 x 平面 [[x1,xp1], [x2,xp2] .. ]
        -------

        """
        return [p.project_to_xxp_plane() for p in phase_space_particles]

    @staticmethod
    def phase_space_particles_project_to_yyp_plane(
        phase_space_particles: List,
    ) -> List[P2]:
        """
        相空间粒子群投影到 y 平面
        Parameters
        ----------
        phase_space_particles 相空间粒子群

        Returns 相空间粒子群投影到 y 平面 [[y1,yp1], [y2,yp2] .. ]
        -------

        """
        return [p.project_to_yyp_plane() for p in phase_space_particles]

    @staticmethod
    def phase_space_particles_project_to_plane(
        phase_space_particles: List, plane_id: int
    ) -> List[P2]:
        """
        相空间粒子群投影到 x/y 平面
        Parameters
        ----------
        phase_space_particles 相空间粒子群
        plane_id 投影到 x 或 y 平面

        Returns 相空间粒子群投影到 x/y 平面
        -------

        """
        if plane_id == PhaseSpaceParticle.XXP_PLANE:
            return PhaseSpaceParticle.phase_space_particles_project_to_xxp_plane(
                phase_space_particles
            )
        elif plane_id == PhaseSpaceParticle.YYP_PLANE:
            return PhaseSpaceParticle.phase_space_particles_project_to_yyp_plane(
                phase_space_particles
            )
        else:
            raise ValueError(f"没有处理plane_id({plane_id})的方法")

    @staticmethod
    def create_from_running_particle(
        ideal_particle: RunningParticle,
        coordinate_system: LocalCoordinateSystem,
        running_particle: RunningParticle,
    ):
        # x y z
        relative_position = running_particle.position - ideal_particle.position
        x = coordinate_system.XI * relative_position
        y = coordinate_system.YI * relative_position
        z = coordinate_system.ZI * relative_position

        # xp yp
        relative_velocity = running_particle.velocity - ideal_particle.velocity
        xp = (coordinate_system.XI * relative_velocity) / ideal_particle.speed
        yp = (coordinate_system.YI * relative_velocity) / ideal_particle.speed

        # delta
        rm = running_particle.compute_scalar_momentum()
        im = ideal_particle.compute_scalar_momentum()
        delta = (rm - im) / im

        return PhaseSpaceParticle(x, xp, y, yp, z, delta)

    @staticmethod
    def create_from_running_particles(
        ideal_particle: RunningParticle,
        coordinate_system: LocalCoordinateSystem,
        running_particles: List[RunningParticle],
    ) -> List:
        return [
            PhaseSpaceParticle.create_from_running_particle(
                ideal_particle, coordinate_system, rp
            )
            for rp in running_particles
        ]

    @staticmethod
    def convert_delta_from_momentum_dispersion_to_energy_dispersion(
        phaseSpaceParticle, centerKineticEnergy_MeV
    ):
        """
        动量分散改动能分散
        Parameters
        ----------
        phaseSpaceParticle 原粒子
        centerKineticEnergy_MeV 中心动能，如 250

        Returns 动量分散改动能分散后的粒子
        -------

        """
        copied: PhaseSpaceParticle = phaseSpaceParticle.copy()
        deltaMomentumDispersion = copied.delta
        deltaEnergyDispersion = (
            Protons.convert_momentum_dispersion_to_energy_dispersion(
                deltaMomentumDispersion, centerKineticEnergy_MeV
            )
        )

        copied.delta = deltaEnergyDispersion

        return copied

    @staticmethod
    def convert_delta_from_momentum_dispersion_to_energy_dispersion_for_list(
        phaseSpaceParticles: List, centerKineticEnergy_MeV
    ):
        """
        动量分散改动能分散，见上方法 convert_delta_from_momentum_dispersion_to_energy_dispersion
        Parameters
        ----------
        phaseSpaceParticles
        centerKineticEnergy_MeV

        Returns
        -------

        """
        return [
            PhaseSpaceParticle.convert_delta_from_momentum_dispersion_to_energy_dispersion(
                pp, centerKineticEnergy_MeV
            )
            for pp in phaseSpaceParticles
        ]

    @staticmethod
    def convert_delta_from_energy_dispersion_to_energy_dispersion_momentum_dispersion(
        phaseSpaceParticle, centerKineticEnergy_MeV: float
    ):
        copied = phaseSpaceParticle.copy()

        EnergyDispersion = copied.getDelta()

        MomentumDispersion = Protons.convert_energy_dispersion_to_momentum_dispersion(
            EnergyDispersion, centerKineticEnergy_MeV
        )

        copied.delta = MomentumDispersion

        return copied

    @staticmethod
    def convert_delta_from_energy_dispersion_to_energy_dispersion_momentum_dispersion_for_list(
        phaseSpaceParticles: List, centerKineticEnergy_MeV: float
    ):
        return [
            PhaseSpaceParticle.convert_delta_from_energy_dispersion_to_energy_dispersion_momentum_dispersion(
                pp, centerKineticEnergy_MeV
            )
            for pp in phaseSpaceParticles
        ]

    def __str__(self) -> str:
        return (
            f"x={self.x},xp={self.xp},y={self.y},yp={self.yp},z={self.z},d={self.delta}"
        )

    def copy(self):
        return PhaseSpaceParticle(self.x, self.xp, self.y, self.yp, self.z, self.delta)

    def getDelta(self):
        return self.delta


class ParticleFactory:
    """
    质子工厂
    """

    @staticmethod
    def create_proton(
        position: P3, direct: P3, kinetic_MeV: float = 250
    ) -> RunningParticle:
        # 速率
        speed = LIGHT_SPEED * math.sqrt(
            1.0
            - (Protons.STATIC_ENERGY_MeV / (Protons.STATIC_ENERGY_MeV + kinetic_MeV))
            ** 2
        )

        # mass kg
        relativistic_mass = Protons.STATIC_MASS_KG / math.sqrt(
            1.0 - (speed ** 2) / (LIGHT_SPEED ** 2)
        )

        return RunningParticle(
            position,
            direct.copy().change_length(speed),
            relativistic_mass,
            Protons.CHARGE_QUANTITY,
            speed,
        )

    @staticmethod
    def create_proton_by_position_and_velocity(
        position: P3, velocity: P3
    ) -> RunningParticle:
        speed = velocity.length()

        relativistic_mass = 0.0

        try:
            relativistic_mass = Protons.STATIC_MASS_KG / math.sqrt(
                1.0 - (speed ** 2) / (LIGHT_SPEED ** 2)
            )
        except RuntimeWarning as e:
            print(
                f"ParticleFactory::create_proton_by_position_and_velocity 莫名其妙的异常 speed={speed} LIGHT_SPEED={LIGHT_SPEED} e={e}"
            )

        return RunningParticle(
            position, velocity, relativistic_mass, Protons.CHARGE_QUANTITY, speed
        )

    @staticmethod
    def create_from_phase_space_particle(
        ideal_particle: RunningParticle,
        coordinate_system: LocalCoordinateSystem,
        phase_space_particle: PhaseSpaceParticle,
    ) -> RunningParticle:
        """
        通过理想粒子，相空间坐标系 和 相空间粒子，来创造粒子
        Parameters
        ----------
        ideal_particle 理想粒子
        coordinate_system 相空间坐标系
        phase_space_particle 相空间粒子

        Returns 通过理想粒子，相空间坐标系 和 相空间粒子，来创造粒子
        -------

        """
        x = phase_space_particle.x
        xp = phase_space_particle.xp
        y = phase_space_particle.y
        yp = phase_space_particle.yp
        z = phase_space_particle.z
        delta = phase_space_particle.delta

        p = ideal_particle.copy()
        # 知道 LocalCoordinateSystem 的用处了吧
        p.position += coordinate_system.XI * x
        p.position += coordinate_system.YI * y
        p.position += coordinate_system.ZI * z

        if delta != 0.0:
            scalar_momentum = p.compute_scalar_momentum() * (1.0 + delta)
            p.change_scalar_momentum(scalar_momentum)  # 这个方法就是为了修改动量而写的

        p.velocity += coordinate_system.XI * (xp * p.speed)
        p.velocity += coordinate_system.YI * (yp * p.speed)

        return p

    @staticmethod
    def create_from_phase_space_particles(
        ideal_particle: RunningParticle,
        coordinate_system: LocalCoordinateSystem,
        phase_space_particles: List[PhaseSpaceParticle],
    ) -> List[RunningParticle]:
        return [
            ParticleFactory.create_from_phase_space_particle(
                ideal_particle, coordinate_system, p
            )
            for p in phase_space_particles
        ]


class CCT(Magnet):
    """
    表示一层弯曲 CCT 线圈
    """

    def __init__(
        self,
        # CCT 局部坐标系
        local_coordinate_system: LocalCoordinateSystem,
        # 大半径：偏转半径
        big_r: float,
        # 小半径（孔径/2）
        small_r: float,
        # 偏转角度，即 phi0*winding_number，典型值 67.5
        bending_angle: float,
        # 各极倾斜角，典型值 [30,90,90,90]
        tilt_angles: List[float],
        # 匝数
        winding_number: int,
        # 电流
        current: float,
        # CCT 路径在二维 ξ-φ 坐标系中的起点
        starting_point_in_ksi_phi_coordinate: P2,
        # CCT 路径在二维 ξ-φ 坐标系中的终点
        end_point_in_ksi_phi_coordinate: P2,
        # 每匝线圈离散电流元数目，数字越大计算精度越高
        disperse_number_per_winding: int = 120,
    ):
        self.local_coordinate_system = local_coordinate_system
        self.big_r = float(big_r)
        self.small_r = float(small_r)
        self.bending_angle = float(bending_angle)
        self.tilt_angles = [float(e) for e in tilt_angles]
        self.winding_number = int(winding_number)
        self.current = float(current)
        self.starting_point_in_ksi_phi_coordinate = starting_point_in_ksi_phi_coordinate
        self.end_point_in_ksi_phi_coordinate = end_point_in_ksi_phi_coordinate
        self.disperse_number_per_winding = int(disperse_number_per_winding)

        # 弯转角度，弧度制
        self.bending_radian = BaseUtils.angle_to_radian(bending_angle)

        # 倾斜角，弧度制
        self.tilt_radians = BaseUtils.angle_to_radian(tilt_angles)

        # 每绕制一匝，φ 方向前进长度
        self.phi0 = self.bending_radian / winding_number

        # 极点 a
        self.a = math.sqrt(big_r ** 2 - small_r ** 2)

        # 双极坐标系另一个常量 η
        self.eta = 0.5 * math.log((big_r + self.a) / (big_r - self.a))

        # 建立 ξ-φ 坐标到三维 xyz 坐标的转换器
        self.bipolar_toroidal_coordinate_system = CCT.BipolarToroidalCoordinateSystem(
            self.a, self.eta, big_r, small_r
        )

        # CCT 路径的在 ξ-φ 坐标的表示 函数 φ(ξ)
        self.phi_ksi_function = lambda ksi: self.__phi_ksi_function(ksi)

        # CCT 路径的在 ξ-φ 坐标的表示 函数 P(ξ)=(ξ,φ(ξ))
        self.p2_function = lambda ksi: P2(ksi, self.phi_ksi_function(ksi))

        # CCT 路径的在 xyz 坐标的表示 函数 P(ξ)=P(x(ξ),y(ξ),z(ξ))
        self.p3_function = lambda ksi: self.bipolar_toroidal_coordinate_system.convert(
            self.p2_function(ksi)
        )

        # 总匝数
        total_disperse_number = winding_number * disperse_number_per_winding

        dispersed_path2: List[List[float]] = [
            self.p2_function(ksi).to_list()
            for ksi in BaseUtils.linspace(
                self.starting_point_in_ksi_phi_coordinate.x,
                self.end_point_in_ksi_phi_coordinate.x,
                total_disperse_number + 1,
            )  # +1 为了满足分段正确性，即匝数 m，需要用 m+1 个点
        ]

        dispersed_path3: List[List[float]] = [
            self.p3_function(ksi).to_list()
            for ksi in BaseUtils.linspace(
                self.starting_point_in_ksi_phi_coordinate.x,
                self.end_point_in_ksi_phi_coordinate.x,
                total_disperse_number + 1,
            )  # +1 为了满足分段正确性，见上
        ]

        # 为了速度，转为 numpy
        self.dispersed_path2: numpy.ndarray = numpy.array(dispersed_path2)
        self.dispersed_path3: numpy.ndarray = numpy.array(dispersed_path3)

        # 电流元 current * (p[i+1] - p[i])
        self.elementary_current = current * (
            self.dispersed_path3[1:] - self.dispersed_path3[:-1]
        )

        # 电流元的位置 (p[i+1]+p[i])/2
        self.elementary_current_position = 0.5 * (
            self.dispersed_path3[1:] + self.dispersed_path3[:-1]
        )

    def __phi_ksi_function(self, ksi: float) -> float:
        x1 = self.starting_point_in_ksi_phi_coordinate.x
        y1 = self.starting_point_in_ksi_phi_coordinate.y
        x2 = self.end_point_in_ksi_phi_coordinate.x
        y2 = self.end_point_in_ksi_phi_coordinate.y

        k = (y2 - y1) / (x2 - x1)
        b = -k * x1 + y1

        phi = k * ksi + b
        for i in range(len(self.tilt_radians)):
            if BaseUtils.equal(self.tilt_angles[i], 90.0):
                continue
            phi += (
                (1 / math.tan(self.tilt_radians[i]))
                / ((i + 1) * math.sinh(self.eta))
                * math.sin((i + 1) * ksi)
            )
        return phi

    class BipolarToroidalCoordinateSystem:
        def __init__(self, a: float, eta: float, big_r: float, small_r: float):
            self.a = a
            self.eta = eta
            self.big_r = big_r
            self.small_r = small_r

            BaseUtils.equal(
                big_r,
                math.sqrt(a * a / (1 - 1 / math.pow(math.cosh(eta), 2))),
                msg=f"BipolarToroidalCoordinateSystem:init 错误1 a({a})eta({eta})R({big_r})r({small_r})",
            )

            BaseUtils.equal(
                small_r,
                big_r / math.cosh(eta),
                msg=f"BipolarToroidalCoordinateSystem:init 错误2 a({a})eta({eta})R({big_r})r({small_r})",
            )

        def convert(self, p: P2) -> P3:
            ksi = p.x
            phi = p.y
            temp = self.a / (math.cosh(self.eta) - math.cos(ksi))
            return P3(
                temp * math.sinh(self.eta) * math.cos(phi),
                temp * math.sinh(self.eta) * math.sin(phi),
                temp * math.sin(ksi),
            )

        def main_normal_direction_at(self, p: P2) -> P3:
            phi = p.y

            center = P3(self.big_r * math.cos(phi), self.big_r * math.sin(phi), 0)

            face_point = self.convert(p)

            return (face_point - center).normalize()

        def __str__(self):
            return f"BipolarToroidalCoordinateSystem a({self.a})eta({self.eta})R({self.big_r})r({self.small_r})"

    def magnetic_field_at(self, point: P3) -> P3:
        if GPU_ON:
            return self.magnetic_field_at_gpu(point)
        else:
            return self.magnetic_field_at_cpu(point)

    def magnetic_field_at_cpu(self, point: P3) -> P3:
        # point 转为局部坐标，并变成 numpy 向量
        p = numpy.array(
            self.local_coordinate_system.point_to_local_coordinate(point).to_list()
        )

        # 点 p 到电流元中点
        r = p - self.elementary_current_position

        # 点 p 到电流元中点的距离的三次方
        rr = (numpy.linalg.norm(r, ord=2, axis=1) ** (-3)).reshape((r.shape[0], 1))

        # 计算每个电流元在 p 点产生的磁场 (此时还没有乘系数 μ0/4π )
        dB = numpy.cross(self.elementary_current, r) * rr

        # 求和，即得到磁场，记得乘以系数 μ0/4π = 1e-7
        B = numpy.sum(dB, axis=0) * 1e-7

        return P3.from_numpy_ndarry3(B)

    def magnetic_field_at_gpu(self, point: P3) -> P3:
        mod = SourceModule(
            """
        #include <stdio.h>
        #include <math.h>
        #include "cuda.h"

        #define MM (0.001f)
        #define DIM (3)
        #define PI (3.1415927f)
        #define X (0)
        #define Y (1)
        #define Z (2)


        __device__ __forceinline__ void vct_cross(float *a, float *b, float *ret) {
            ret[X] = a[Y] * b[Z] - a[Z] * b[Y];
            ret[Y] = -a[X] * b[Z] + a[Z] * b[X];
            ret[Z] = a[X] * b[Y] - a[Y] * b[X];
        }

        __device__ __forceinline__ void vct_add_local(float *a_local, float *b) {
            a_local[X] += b[X];
            a_local[Y] += b[Y];
            a_local[Z] += b[Z];
        }

        __device__ __forceinline__ void vct_add(float *a, float *b, float *ret) {
            ret[X] = a[X] + b[X];
            ret[Y] = a[Y] + b[Y];
            ret[Z] = a[Z] + b[Z];
        }

        __device__ __forceinline__ void vct_dot_a_v(float a, float *v) {
            v[X] *= a;
            v[Y] *= a;
            v[Z] *= a;
        }

        __device__ __forceinline__ void vct_dot_a_v_ret(float a, float *v, float *ret) {
            ret[X] = v[X] * a;
            ret[Y] = v[Y] * a;
            ret[Z] = v[Z] * a;
        }

        __device__ __forceinline__ void vct_copy(float *src, float *des) {
            des[X] = src[X];
            des[Y] = src[Y];
            des[Z] = src[Z];
        }

        __device__ __forceinline__ float vct_len(float *v) {
            return sqrtf(v[X] * v[X] + v[Y] * v[Y] + v[Z] * v[Z]);
        }

        __device__ __forceinline__ void vct_zero(float *v) {
            v[X] = 0.0f;
            v[Y] = 0.0f;
            v[Z] = 0.0f;
        }

        __device__ __forceinline__ void vct_sub(float *a, float *b, float *ret) {
            ret[X] = a[X] - b[X];
            ret[Y] = a[Y] - b[Y];
            ret[Z] = a[Z] - b[Z];
        }

        // 磁场计算 注意，这里计算的不是电流元的磁场，还需要乘以 电流 和 μ0/4π (=1e-7)
        __device__ void dB(float *p0, float *p1, float *p, float *ret) {
            float p01[DIM];
            float r[DIM];
            float rr;

            vct_sub(p1, p0, p01); // p01 = p1 - p0

            vct_add(p0, p1, r); // r = p0 + p1

            vct_dot_a_v(0.5f, r); // r = (p0 + p1)/2

            vct_sub(p, r, r); // r = p - r

            rr = vct_len(r); // rr = len(r)

            vct_cross(p01, r, ret); // ret = p01 x r

            rr = 1.0f / rr / rr / rr; // changed

            vct_dot_a_v(rr, ret); // rr . (p01 x r)
        }

        __global__ void magnet(float *winding, float *p, int *length, float *ret) {
            int tid = threadIdx.x + blockIdx.x * blockDim.x;

            if (tid == 0)vct_zero(ret);

            __syncthreads();

            if (tid < *length - 1) {
                float *p0 = winding + tid * DIM;
                float *p1 = winding + (tid + 1) * DIM;
                float db[3];

                dB(p0, p1, p, db);

                atomicAdd(&ret[X], db[X]);
                atomicAdd(&ret[Y], db[Y]);
                atomicAdd(&ret[Z], db[Z]);
            }
        }"""
        )

        magnet = mod.get_function("magnet")

        # point 转为局部坐标，并变成 numpy 向量
        p = numpy.array(
            self.local_coordinate_system.point_to_local_coordinate(point).to_list()
        )

        length = int(self.dispersed_path3.shape[0])
        winding = self.dispersed_path3.flatten().astype(numpy.float32)
        ret = numpy.empty((3,), dtype=numpy.float32)

        magnet(
            drv.In(winding),
            drv.In(p),
            drv.In(numpy.array([length]).astype(numpy.int32)),
            drv.Out(ret),
            block=(512, 1, 1),
            grid=((length + 511) // 512, 1),
        )

        return P3.from_numpy_ndarry3(ret * self.current * 1e-7)

    def __str__(self):
        return (
            f"CCT: local_coordinate_system({self.local_coordinate_system})big_r({self.big_r})small_r({self.small_r})"
            + f"bending_angle({self.bending_angle})tilt_angles({self.tilt_angles})winding_number({self.winding_number})"
            + f"current({self.current})start_ksi({self.start_ksi})start_phi({self.start_phi})clockwise({self.clockwise})"
            + f"disperse_number_per_winding({self.disperse_number_per_winding})"
        )


class QS(Magnet):
    pass


class Plot3:
    INIT: bool = False
    ax = None

    @staticmethod
    def __init():
        plt.rcParams["font.sans-serif"] = ["SimHei"]  # 用来正常显示中文标签
        plt.rcParams["axes.unicode_minus"] = False  # 用来正常显示负号

        fig = plt.figure()
        Plot3.ax = fig.gca(projection="3d")
        Plot3.ax.grid(False)

        Plot3.INIT = True

    @staticmethod
    def plot_xyz(x, y, z, describe="r") -> None:
        if not Plot3.INIT:
            Plot3.__init()

        Plot3.ax.plot(x, y, z, describe)

    @staticmethod
    def plot_p3(p: P3, describe="r") -> None:
        if not Plot3.INIT:
            Plot3.__init()
        Plot3.ax.plot(p.x, p.y, p.z, describe)

    @staticmethod
    def plot_p3s(ps: List[P3], describe="r") -> None:
        if not Plot3.INIT:
            Plot3.__init()

        Plot3.ax.plot([p.x for p in ps], [p.y for p in ps], [p.z for p in ps], describe)

    @staticmethod
    def plot3d(line: Line2, step: float = 1 * MM, describe="r"):
        Plot3.plot_p3s(line.disperse3d(step), describe)

    @staticmethod
    def plot_ndarry3ds(narray: numpy.ndarray, describe="r-") -> None:
        if not Plot3.INIT:
            Plot3.__init()
        x = narray[:, 0]
        y = narray[:, 1]
        z = narray[:, 2]
        plt.plot(x, y, z, describe)

    @staticmethod
    def plot_cct(cct: CCT, describe="r-") -> None:
        if not Plot3.INIT:
            Plot3.__init()
        cct_path3: numpy.ndarray = cct.dispersed_path3
        Plot3.plot_ndarry3ds(cct_path3)

    @staticmethod
    def show():
        if not Plot3.INIT:
            raise RuntimeError("Plot3::请在show前调用plot")

        plt.show()


class Plot2:
    INIT = False

    @staticmethod
    def __init():
        plt.rcParams["font.sans-serif"] = ["SimHei"]  # 用来正常显示中文标签
        plt.rcParams["axes.unicode_minus"] = False  # 用来正常显示负号

        Plot2.INIT = True

    @staticmethod
    def plot_xy(x: float, y: float, describe="r") -> None:
        if not Plot2.INIT:
            Plot2.__init()

        plt.plot(x, y, describe)

    @staticmethod
    def plot_p2(p: P2, describe="r") -> None:
        if not Plot2.INIT:
            Plot2.__init()

        plt.plot(p.x, p.y, describe)

    @staticmethod
    def plot_p2s(ps: List[P2], describe="r-") -> None:
        if not Plot2.INIT:
            Plot2.__init()

        plt.plot([p.x for p in ps], [p.y for p in ps], describe)

    @staticmethod
    def plot_ndarry2ds(narray: numpy.ndarray, describe="r-") -> None:
        if not Plot2.INIT:
            Plot2.__init()
        x = narray[:, 0]
        y = narray[:, 1]
        plt.plot(x, y, describe)

    @staticmethod
    def plot_cct(cct: CCT, describe="r-") -> None:
        if not Plot2.INIT:
            Plot2.__init()
        cct_path2: numpy.ndarray = cct.dispersed_path2
        Plot2.plot_ndarry2ds(cct_path2)

    @staticmethod
    def show():
        if not Plot2.INIT:
            raise RuntimeError("Plot2::请在show前调用plot")

        plt.show()
