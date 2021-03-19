
import os
import sys
curPath = os.path.abspath(os.path.dirname(__file__))
rootPath = os.path.split(curPath)[0]
PathProject = os.path.split(rootPath)[0]
sys.path.append(rootPath)
sys.path.append(PathProject)


from cosy_utils import *
from cctpy import *

try:
    from books.cct.cctpy.cctpy import *
    from books.cct.cctpy.cctpy_ext import *
    from books.cct.cctpy.cosy_utils import *
except:
    pass


if __name__ == "__main__":
    BaseUtils.i_am_sure_my_code_closed_in_if_name_equal_main()
    if False:  # 圆周均匀分布
        pc = [BaseUtils.Random.uniformly_distributed_along_circumference()
              for i in range(100)]

        Plot2.plot(pc, describe='r.')
        Plot2.equal()
        Plot2.show()

    if False:  # 圆面均匀分布
        pc = [BaseUtils.Random.uniformly_distributed_in_circle()
              for i in range(500)]

        Plot2.plot(pc, describe='r.')
        Plot2.equal()
        Plot2.show()

    if False:  # 球面均匀分布
        pc = [BaseUtils.Random.uniformly_distributed_at_spherical_surface()
              for i in range(500)]

        Plot3.plot_p3s(pc, describe='r.')
        Plot3.set_center(cube_size=2)
        Plot3.show()

    if False:  # 球内均匀分布
        pc = [BaseUtils.Random.uniformly_distributed_in_sphere()
              for i in range(500)]

        Plot3.plot_p3s(pc, describe='r.')
        Plot3.set_center(cube_size=2)
        Plot3.show()

    if False:  # 椭圆圆周均匀分布
        pc = [BaseUtils.Random.uniformly_distributed_along_elliptic_circumference(
            3.5, 7.5) for i in range(100)]

        Plot2.plot(pc, describe='r.')
        Plot2.equal()
        Plot2.show()

    if False:  # 椭圆均匀分布
        pc = [BaseUtils.Random.uniformly_distributed_in_ellipse(
            3.5, 7.5) for i in range(500)]

        Plot2.plot(pc, describe='r.')
        Plot2.equal()
        Plot2.show()

    if False:  # 椭球圆均匀分布
        # pc = [BaseUtils.Random.uniformly_distributed_at_ellipsoidal_surface(3.5,7.5,15) for i in range(500)]

        BaseUtils.i_am_sure_my_code_closed_in_if_name_equal_main()
        pc = BaseUtils.submit_process_task(
            BaseUtils.Random.uniformly_distributed_at_ellipsoidal_surface,
            [(3.5, 7.5, 15) for i in range(500)]
        )

        Plot3.plot_p3s(pc, describe='r.')
        Plot3.set_center(cube_size=30)
        Plot3.show()

    if False:  # 椭球内均匀分布 kirei
        BaseUtils.i_am_sure_my_code_closed_in_if_name_equal_main()
        pc = BaseUtils.submit_process_task(
            BaseUtils.Random.uniformly_distributed_in_ellipsoid,
            [(3.5, 7.5, 15) for i in range(10000)]
        )

        Plot3.plot_p3s(pc, describe='r.')
        Plot3.set_center(cube_size=30)
        Plot3.show()

    if False:
        a = 2.5
        b = 7.5

        # e = BaseUtils.Ellipse(A=1/(a**2), C=1/(b**2), B=0.0, D=1.0)
        # print(e.circumference)

        # c =  math.pi*(a+b)*(1+3*((a-b)/(a+b))**2/(10+math.sqrt(4-3*((a-b)/(a+b))**2))+(4/math.pi-14/11)*((a-b)/(a+b))**(14.233+13.981*((a-b)/(a+b))**6.42))
        # print(c)
        # print(c.real)

        print(BaseUtils.Random.hypersphere_volume(3))
        print(4/3*math.pi)

        print(BaseUtils.Random.hypersphere_volume(4))
        print((math.pi**2)/2)

        print(BaseUtils.Random.hypersphere_area(2))
        print(2*math.pi)

        print(BaseUtils.Random.hypersphere_area(3))
        print(4*math.pi)

        print(BaseUtils.Random.hypersphere_area(4))
        print(2*math.pi**2)

    if False:  # 超球体表面分布 d=2
        pc = [P2.from_numpy_ndarry(numpy.array(
            BaseUtils.Random.uniformly_distributed_at_hyperespherical_surface(2, 1))) for i in range(500)]

        Plot2.plot(pc, describe='r.')
        Plot2.equal()
        Plot2.show()

    if False:  # 超球体表面分布
        pc = [P3.from_numpy_ndarry(numpy.array(
            BaseUtils.Random.uniformly_distributed_at_hyperespherical_surface(3, 1))) for i in range(500)]

        Plot3.plot_p3s(pc, describe='r.')
        Plot3.set_center(cube_size=30)
        Plot3.show()

    if False:  # 超球体表面分布
        pc = [P3.from_numpy_ndarry(numpy.array(
            BaseUtils.Random.uniformly_distributed_at_hyperespherical_surface(4, 1)[:3])) for i in range(5000)]

        Plot3.plot_p3s(pc, describe='r.')
        Plot3.set_center(cube_size=30)
        Plot3.show()

    if False:  # 超椭球球体表面分布
        pc = [P3.from_numpy_ndarry(numpy.array(BaseUtils.Random.uniformly_distributed_at_ellipsoidal_surface(
            (3.5, 7.5, 10)
        )[:3])) for i in range(500)]

        Plot3.plot_p3s(pc, describe='r.')
        Plot3.set_center(cube_size=30)
        Plot3.show()

    if False:  # 超椭球球分布
        pc = [P3.from_numpy_ndarry(numpy.array(BaseUtils.Random.uniformly_distributed_in_hypereellipsoid(
            (3.5, 7.5, 10)
        )[:3])) for i in range(500)]

        Plot3.plot_p3s(pc, describe='r.')
        Plot3.set_center(cube_size=30)
        Plot3.show()

    if False:  # 超椭球球分布
        pc = [P2.from_numpy_ndarry(numpy.array(BaseUtils.Random.uniformly_distributed_in_hypereellipsoid(
            (3.5, 7.5)
        ))) for i in range(500)]

        Plot2.plot(pc, describe='r.')
        Plot2.equal()
        Plot2.show()

    if False:
        ps = ParticleFactory.distributed_particles(
            3.5*MM, 7.5*MRAD, 3.5*MM, 7.5*MM, 0.08, 20,
            ParticleFactory.DISTRIBUTION_AREA_EDGE,
            x_distributed=True, xp_distributed=True, delta_distributed=True
        )

        for s in SRConvertor.to_cosy_sr(ps):
            print(s)

        pc = [P3(p.y, p.yp, p.delta) for p in ps]

    if False:
        p2s = []
        for _ in range(10000):
            p2s.append(P2(
                BaseUtils.Random.gauss(0.0,3.5/2),
                BaseUtils.Random.gauss(0.0,7.5/2),
            ))
        
        e = BaseUtils.Ellipse.create_standard_ellipse(a=3.5,b=7.5)


        Plot2.plot_p2s(p2s,"r.")
        Plot2.plot(e,describe='k-')
        Plot2.equal()
        Plot2.show()

    if False:
        p2s = [P2.from_list(BaseUtils.Random.gauss_multi_dimension([0,0],[3.5/2,7.5/2])) for _ in range(10000)]
        
        e = BaseUtils.Ellipse.create_standard_ellipse(a=3.5,b=7.5)


        Plot2.plot_p2s(p2s,"r.")
        Plot2.plot(e,describe='k-')
        Plot2.equal()
        Plot2.show()

    if True:
        ps = ParticleFactory.distributed_particles(
            3.5*MM/2, 7.5*MRAD/2, 3.5*MM/2, 7.5*MM/2, 0.08, 2000,
            ParticleFactory.DISTRIBUTION_AREA_FULL,
            x_distributed=True, xp_distributed=True,
            distribution_type='gauss'
        )

        pxs = [P3(p.x, p.xp) for p in ps]

        e = BaseUtils.Ellipse.create_standard_ellipse(a=3.5*MM,b=7.5*MRAD)


        Plot2.plot_p2s(pxs,"r.")
        Plot2.plot(e,describe='k-')
        Plot2.equal()
        Plot2.show()