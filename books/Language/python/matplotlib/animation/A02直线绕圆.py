import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib.animation as animation

fig = plt.figure()
ax = plt.axes()


xs = []
ys = []

for i in range(360):
    xs.append(np.linspace(0, np.cos(i / 180.0 * np.pi)))
    ys.append(np.linspace(0, np.sin(i / 180.0 * np.pi)))

ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.axis("equal")

ax.plot(
    [np.cos(t) for t in np.linspace(0, 2 * np.pi, 360)],
    [np.sin(t) for t in np.linspace(0, 2 * np.pi, 360)],
    "k-",
)

(point_ani,) = plt.plot(xs[0], ys[0], "r-")


def update(index):
    point_ani.set_data(xs[index], ys[index])
    return (point_ani,)


ani = animation.FuncAnimation(fig, update, np.arange(0, 360, 2), interval=20, blit=True)
ani.save('line_round.gif', writer='imagemagick', fps=50)
plt.show()