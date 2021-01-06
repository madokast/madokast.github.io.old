"""
agcct 连接段的构建，尽可能的自动化


@Author 赵润晓
"""

try:
    from books.cct.cctpy.cctpy import *
    from books.cct.cctpy.cctpy_ext import *
except:
    pass

from cctpy import *
from cctpy_ext import *
from typing import List


class AGCCT_CONNECTOR(Magnet):
    def __init__(self, agcct1: CCT, agcct2: CCT) -> None:
        # 必要的验证
        if agcct1.local_coordinate_system != agcct2.local_coordinate_system:
            raise ValueError(
                f"需要连接的cct不属于同一坐标系，无法连接。agccct1={agcct1},agcct2={agcct2}")
        BaseUtils.equal(agcct1.big_r, agcct2.big_r,
                        msg=f"需要连接的 cct big_r 不同，无法连接。agccct1={agcct1},agcct2={agcct2}")
        BaseUtils.equal(agcct1.small_r, agcct2.small_r,
                        msg=f"需要连接的 cct small_r 不同，无法连接。agccct1={agcct1},agcct2={agcct2}")
        BaseUtils.equal(agcct1.current, agcct2.current,
                        msg=f"需要连接的 cct 电流不同，无法连接。agccct1={agcct1},agcct2={agcct2}")

        # 坐标系
        self.local_coordinate_system = agcct1.local_coordinate_system

        # 前 cct 的终点
        pre_cct_end_point_in_ksi_phi_coordinate = agcct1.end_point_in_ksi_phi_coordinate
        # 后 cct 的起点
        next_cct_starting_point_in_ksi_phi_coordinate = agcct2.starting_point_in_ksi_phi_coordinate
        BaseUtils.equal(pre_cct_end_point_in_ksi_phi_coordinate.x, next_cct_starting_point_in_ksi_phi_coordinate.x,
                        msg=f"需要连接的前段 cct 终点 ξ 不等于后段 cct 的起点 ξ，无法连接。agccct1={agcct1},agcct2={agcct2}")

        # 前 cct 的终点，在 φR-ξr 坐标系中的点 a
        # 后 cct 的起点  b 
        # 以及切向 va vb
        a = P2(x=pre_cct_end_point_in_ksi_phi_coordinate.y*agcct1.big_r, y=0.0)
        b = P2(x=next_cct_starting_point_in_ksi_phi_coordinate.y*agcct2.big_r, y=0.0)
        va = BaseUtils.derivative(lambda ksi: agcct1.p2_function(ksi))(
            pre_cct_end_point_in_ksi_phi_coordinate.x)
        vb = BaseUtils.derivative(lambda ksi: agcct2.p2_function(ksi))(
            next_cct_starting_point_in_ksi_phi_coordinate.x)
        print(f"a={a}, b={b}")
        print(f"va={va}, vb={vb}")

        # 开始连接
        # 首先是直线
        StraightLine2.intersecting_point()



if __name__ == "__main__":
    print('hello')
