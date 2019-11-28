import matplotlib as mpl
import numpy as np
import os
pwd = os.getcwd()

def plotBackground(ax, x0, y0, colour):
    cmap=makeCmap(colour,'backgroundColour')
    yMin=y0*np.power(x0/1e20,1)
    yMax=y0*np.power(x0/1e-20,1)
    nColours=25
    change=0.01
    for i in range(nColours):
        ax.fill([1e-20,1e-20,1e20],[yMax*change**i, yMin*change**i, yMin*change**i], c=cmap((i+1)/nColours), zorder=i + 1 - 100)

def makeCmap(hexColour, name, zeroColour='#FFFFFF'):
    r0,g0,b0=mpl.colors.hex2color(zeroColour)
    r,g,b=mpl.colors.hex2color(hexColour)
    cdict = {'red':   ((0.0, r0, r0),
                      (1.0, r, r)),
             'green': ((0.0, g0, g0),
                      (1.0, g, g)),
             'blue':  ((0.0, b0, b0),
                      (1.0, b, b))
            }
    cmap = mpl.colors.LinearSegmentedColormap(name, cdict)
    return cmap


def load_data(fn):
    x, xw, _, y = np.genfromtxt(pwd + '/../data/' + fn).T
    idxs = np.nonzero(y)
    x = x[idxs]
    y = y[idxs]
    return x, y

def plot_hist(ax, data):
    label, fn, cut, c, ls, x, y = data
    x = 10**x
    if cut > 0: ax.plot(x, y, c=c, ls=ls, alpha=0.3)
    ax.plot(x[x > cut], y[x > cut], c=c, ls=ls, label=label)
