from sympy import *


ksi, phi0, PI = symbols("ksi phi0 PI")
print("定义符号 ksi phi0 PI -- ", ksi, phi0, PI, "\n")

phi = Function("f")(ksi) + phi0 / (2 * PI) * ksi
print("定义函数 phi = f(ksi) + phi0/(2*PI)*ksi -- ", phi, "\n")

print("phi 对 ksi 求导 -- ", phi.diff(ksi), "\n")

a, sinh_eta0, cosh_eta0 = symbols("a sinh_eta0 cosh_eta0")
print("定义符号 a,sinh_eta0,cosh_eta0 -- ", a, sinh_eta0, cosh_eta0, "\n")

x = a / (cosh_eta0 - cos(ksi)) * (sinh_eta0 * cos(phi))
y = a / (cosh_eta0 - cos(ksi)) * (sinh_eta0 * sin(phi))
z = a / (cosh_eta0 - cos(ksi)) * sin(ksi)
print("构造路径方程 (x,y,z) -- ", x, y, z, "\n", sep="\n")

vx = x.diff(ksi)
vy = y.diff(ksi)
vz = z.diff(ksi)
print("求路径切向 (vx,vy,vz) -- ", vx, vy, vz, "\n", sep="\n")

axis_x = a * cos(phi)
axis_y = a * sin(phi)
axis_z = 0
print("构造圆环中轴线 (axis) -- ", axis_x, axis_y, axis_z, "\n", sep="\n")


nm_x = x - axis_x
nm_y = y - axis_y
nm_z = z - axis_z
print("求主法线方向 (nm) -- ", nm_x, nm_y, nm_z, "\n", sep="\n")

ns_x = vy * nm_z - vz * nm_y
ns_y = vz * nm_x - vx * nm_z
ns_z = vx * nm_y - vy * nm_x
print("求副法线方向 (ns) -- ", ns_x, ns_y, ns_z, "\n", sep="\n")

ns_len = sqrt(ns_x ** 2 + ns_y ** 2 + ns_z ** 2)
ns_x = ns_x / ns_len
ns_y = ns_y / ns_len
ns_z = ns_z / ns_len
print("求副法线方向归一化 (ns) -- ", ns_x, ns_y, ns_z, "\n", sep="\n")


ksi0, width, t = symbols("ksi0 width t")
print("构建变量 ksi0 width", ksi0, width, "\n")


Lx = (
    x.subs(ksi, ksi0)
    + width / 2 * ns_x.subs(ksi, ksi0)
    - x.subs(ksi, ksi0 + 2 * PI + t)
    + width / 2 * ns_x.subs(ksi, ksi0 + 2 * PI + t)
)
Ly = (
    y.subs(ksi, ksi0)
    + width / 2 * ns_y.subs(ksi, ksi0)
    - y.subs(ksi, ksi0 + 2 * PI + t)
    + width / 2 * ns_y.subs(ksi, ksi0 + 2 * PI + t)
)
Lz = (
    z.subs(ksi, ksi0)
    + width / 2 * ns_z.subs(ksi, ksi0)
    - z.subs(ksi, ksi0 + 2 * PI + t)
    + width / 2 * ns_z.subs(ksi, ksi0 + 2 * PI + t)
)

print("计算 L -- ", Lx, Ly, Lz, "\n", sep="\n")

L_len2 = Lx ** 2 + Ly ** +(Lz ** 2)
print("计算 L^2 -- ", L_len2, "\n")

L_len2dt =  L_len2.diff(t)
print("计算 dL^2/dt -- ", L_len2dt, "\n")


print("============================")

print(solve(L_len2dt,t))