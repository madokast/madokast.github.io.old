from cctpy import M,MM,LocalCoordinateSystem,QS,Plot3,P3

length = 0.2 * M
aper = 30 * MM
g = 10.0
L = 0
lc = LocalCoordinateSystem(P3(),-P3.x_direct(),P3.y_direct())
qs = QS(lc, length, g, L, aper)

Plot3.plot_qs(qs)
Plot3.show()

m = qs.magnetic_field_at_cpu(P3(10 * MM, 0.1, 0.0))
print(m)