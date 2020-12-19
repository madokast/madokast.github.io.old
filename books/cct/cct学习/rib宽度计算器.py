import math
import numpy as np
from matplotlib import pyplot as plt


MM = 0.001

def dicct_outer():

    R = 0.95
    r = 128*MM
    bend_angle = 67.5
    wind_number = 128
    tilt_angles = [30,80.988,94.383,91.650] # [106.654, 30,67.901,90.941]
    rib_width = 4*MM

    ######################################

    a = math.sqrt(R**2-r**2)
    eta = 0.5 * math.log((R + a) / (R - a))
    sheta = math.sinh(eta)
    cheta = math.cosh(eta)
    bend_r = bend_angle/180*math.pi
    phi0 = bend_r/wind_number
    tilt_rs = [ta/180*math.pi for ta in tilt_angles]


    def phi_ksi_fun(ksi):
        phi = phi0/2/math.pi*ksi
        for i in range(len(tilt_rs)):
            phi += (
                (1 / math.tan(tilt_rs[i]))
                / ((i + 1) * sheta)
                * math.sin((i + 1) * ksi)
            )
        return phi


    delta = 1e-6


    def phi_ksi_fun_diff(ksi):
        return (phi_ksi_fun(ksi+delta)-phi_ksi_fun(ksi))/delta


    # 博士论文定义的 k
    def k(ksi):
        return cheta - math.cos(ksi)


    def wind_distance(ksi):
        return a * phi0 * (1/k(ksi)) * ((
            ((sheta)**(-2)) + ((phi_ksi_fun_diff(ksi))**2)
        )**-0.5)


    # plot
    ksis = np.linspace(0, 2*math.pi, 360)
    wind_dis = [wind_distance(ksi)/MM-rib_width/MM for ksi in ksis]
    plt.plot(ksis/math.pi*180,wind_dis,'k-')
    plt.xticks(np.arange(0,360+60,60))


def dicct_inner():

    R = 0.95
    r = 113*MM
    bend_angle = 67.5
    wind_number = 128
    tilt_angles = [30,80.988,94.383,91.650] # [106.654, 30,67.901,90.941]
    rib_width = 4*MM

    ######################################

    a = math.sqrt(R**2-r**2)
    eta = 0.5 * math.log((R + a) / (R - a))
    sheta = math.sinh(eta)
    cheta = math.cosh(eta)
    bend_r = bend_angle/180*math.pi
    phi0 = bend_r/wind_number
    tilt_rs = [ta/180*math.pi for ta in tilt_angles]


    def phi_ksi_fun(ksi):
        phi = phi0/2/math.pi*ksi
        for i in range(len(tilt_rs)):
            phi += (
                (1 / math.tan(tilt_rs[i]))
                / ((i + 1) * sheta)
                * math.sin((i + 1) * ksi)
            )
        return phi


    delta = 1e-6


    def phi_ksi_fun_diff(ksi):
        return (phi_ksi_fun(ksi+delta)-phi_ksi_fun(ksi))/delta


    # 博士论文定义的 k
    def k(ksi):
        return cheta - math.cos(ksi)


    def wind_distance(ksi):
        return a * phi0 * (1/k(ksi)) * ((
            ((sheta)**(-2)) + ((phi_ksi_fun_diff(ksi))**2)
        )**-0.5)


    # plot
    ksis = np.linspace(0, 2*math.pi, 360)
    wind_dis = [wind_distance(ksi)/MM-rib_width/MM for ksi in ksis]
    plt.plot(ksis/math.pi*180,wind_dis,'k-')
    plt.xticks(np.arange(0,360+60,60))


def agcct_inner():

    R = 0.95
    r = 83*MM
    bend_angle = 67.5
    wind_number = 128
    tilt_angles = [106.654, 30,67.901,90.941]
    rib_width = 4*MM

    ######################################

    a = math.sqrt(R**2-r**2)
    eta = 0.5 * math.log((R + a) / (R - a))
    sheta = math.sinh(eta)
    cheta = math.cosh(eta)
    bend_r = bend_angle/180*math.pi
    phi0 = bend_r/wind_number
    tilt_rs = [ta/180*math.pi for ta in tilt_angles]


    def phi_ksi_fun(ksi):
        phi = phi0/2/math.pi*ksi
        for i in range(len(tilt_rs)):
            phi += (
                (1 / math.tan(tilt_rs[i]))
                / ((i + 1) * sheta)
                * math.sin((i + 1) * ksi)
            )
        return phi


    delta = 1e-6


    def phi_ksi_fun_diff(ksi):
        return (phi_ksi_fun(ksi+delta)-phi_ksi_fun(ksi))/delta


    # 博士论文定义的 k
    def k(ksi):
        return cheta - math.cos(ksi)


    def wind_distance(ksi):
        return a * phi0 * (1/k(ksi)) * ((
            ((sheta)**(-2)) + ((phi_ksi_fun_diff(ksi))**2)
        )**-0.5)


    # plot
    ksis = np.linspace(0, 2*math.pi, 360)
    wind_dis = [wind_distance(ksi)/MM-rib_width/MM for ksi in ksis]
    plt.plot(ksis/math.pi*180,wind_dis,'r-')
    plt.xticks(np.arange(0,360+60,60))

def agcct_outer():

    R = 0.95
    r = 98*MM
    bend_angle = 67.5
    wind_number = 128
    tilt_angles = [106.654, 30,67.901,90.941]
    rib_width = 4*MM

    ######################################

    a = math.sqrt(R**2-r**2)
    eta = 0.5 * math.log((R + a) / (R - a))
    sheta = math.sinh(eta)
    cheta = math.cosh(eta)
    bend_r = bend_angle/180*math.pi
    phi0 = bend_r/wind_number
    tilt_rs = [ta/180*math.pi for ta in tilt_angles]


    def phi_ksi_fun(ksi):
        phi = phi0/2/math.pi*ksi
        for i in range(len(tilt_rs)):
            phi += (
                (1 / math.tan(tilt_rs[i]))
                / ((i + 1) * sheta)
                * math.sin((i + 1) * ksi)
            )
        return phi


    delta = 1e-6


    def phi_ksi_fun_diff(ksi):
        return (phi_ksi_fun(ksi+delta)-phi_ksi_fun(ksi))/delta


    # 博士论文定义的 k
    def k(ksi):
        return cheta - math.cos(ksi)


    def wind_distance(ksi):
        return a * phi0 * (1/k(ksi)) * ((
            ((sheta)**(-2)) + ((phi_ksi_fun_diff(ksi))**2)
        )**-0.5)


    # plot
    ksis = np.linspace(0, 2*math.pi, 360)
    wind_dis = [wind_distance(ksi)/MM-rib_width/MM for ksi in ksis]
    plt.plot(ksis/math.pi*180,wind_dis,'r-')
    plt.xticks(np.arange(0,360+60,60))



dicct_inner()
dicct_outer()
agcct_inner()
agcct_outer()

plt.xlabel('angle/deg')
plt.ylabel('rib width/mm')
plt.legend(['dicct_in','dicct_out','agcct_in','agcct_out'])

plt.show()
