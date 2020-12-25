import time
from cctpy import *
import sys

from numpy.core.defchararray import array
sys.path.append(
    r'C:\Users\madoka_9900\Documents\github\madokast.github.io\books\cct\cctpy')

try:
    from ..cctpy import *
except ImportError:
    pass


beamline = HUST_SC_GANTRY.beamline
beamline_length_part1 = HUST_SC_GANTRY.beamline_length_part1


p0 = ParticleFactory.create_proton_along(
    beamline.trajectory, s=beamline_length_part1, kinetic_MeV=215)

ip = ParticleFactory.create_proton_along(
    beamline.trajectory, s=beamline.get_length(), kinetic_MeV=215)

footsteps = numpy.array([1*MM, 2*MM, 5*MM, 10*MM, 20*MM, 50*MM])

# for s in footsteps:
#     """
#     步长0.001，用时51.842s，x=0.0036250544277036346
#     步长0.002，用时26.128s，x=0.004050980069884211
#     步长0.005，用时10.39s，x=0.005258877946395078
#     步长0.01，用时5.23s，x=0.006977427073071076
#     步长0.02，用时2.635s，x=0.009755044522043907
#     步长0.05，用时1.0442s，x=0.017674440504493513
#     """
#     start = time.time()
#     p = p0.copy()
#     ParticleRunner.run_only(p,beamline,beamline.get_length()-beamline_length_part1,footstep=s)
#     pp = PhaseSpaceParticle.create_from_running_particle(ip,ip.get_natural_coordinate_system(),p)

#     print(f"步长{s:.5}，用时{(time.time()-start):.5}s，x={pp.x}")


def runge_kutta4(t0, y0, y_derived_function, footstep, foot_number):
    for i in range(foot_number):
        k1 = y_derived_function(t0, y0)
        k2 = y_derived_function(t0+footstep/2, y0+footstep/2*k1)
        k3 = y_derived_function(t0+footstep/2, y0+footstep/2*k2)
        k4 = y_derived_function(t0+footstep, y0+footstep*k3)

        t0 += footstep
        y0 += (footstep/6)*(k1+2*k2+2*k3+k4)

    return y0


def func(t, Y):
    pos = Y[0]
    v = Y[1]
    a = (p.e/p.relativistic_mass)*(v@beamline.magnetic_field_at(pos))
    return numpy.array([v, a])


t_end = (beamline.get_length()-beamline_length_part1)/p0.speed
for ds in footsteps[::-1]:
    """
    步长0.05，用时3.8517s，x=0.003614891388300663
    步长0.02，用时9.6791s，x=0.003718388644784866
    步长0.01，用时19.309s，x=0.0033095199856791545
    步长0.005，用时38.551s，x=0.003205233213714436
    步长0.002，用时97.683s，x=0.0032450857436989427
    步长0.001，用时196.94s，x=0.003224365499737452
    """
    p = p0.copy()

    Y0 = numpy.array([p.position, p.velocity])
    t0 = 0
    dt = ds / p.speed

    start = time.time()
    Y1 = runge_kutta4(t0, Y0, func, dt, math.ceil(t_end/dt))
    p_end = ParticleFactory.create_proton_by_position_and_velocity(
        position=Y1[0], velocity=Y1[1])
    pp = PhaseSpaceParticle.create_from_running_particle(
        ip, ip.get_natural_coordinate_system(), p_end)
    print(f"步长{ds:.5}，用时{(time.time()-start):.5}s，x={pp.x}")
