import numpy as np

import matplotlib.pyplot as plt

from scipy.interpolate import interpn



def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :

    """

    Scatter plot colored by 2d histogram

    """

    if ax is None :

        fig , ax = plt.subplots()

    data , x_e, y_e = np.histogram2d( x, y, bins = bins)

    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )



    # Sort the points by density, so that the densest points are plotted last

    if sort :

        idx = z.argsort()

        x, y, z = x[idx], y[idx], z[idx]



    ax.scatter( x, y, c=z, **kwargs )

    return ax





if "__main__" == __name__ :

    data = np.loadtxt("out.txt")
    x = data[:,0]
    y = data[:,1]
    print(data)
    print(x)
    print(y)


    # x = np.random.normal(size=100000)

    # y = x * 3 + np.random.normal(size=100000)

    density_scatter( x, y, bins = [30,30] )

    lim = 15
    plt.xlim(-lim,lim)
    plt.ylim(-lim,lim)
    font_label = {
        "family": "Times New Roman",
        "weight": "normal",
        "size": 50,
    }
    plt.xlabel(xlabel='x/mm', fontdict=font_label)
    plt.ylabel(ylabel='y/mm', fontdict=font_label)
    plt.title(label='', fontdict=font_label)

    plt.xticks(fontproperties="Times New Roman", size=50)
    plt.yticks(fontproperties="Times New Roman", size=50)

    plt.show()