import numpy as np

def plot_background(ax, xlo, xhi, ylo, yhi, logslope, colour):
    xs = np.linspace(xlo, xhi, xhi - xlo + 1)
    ys = np.linspace(ylo, yhi, yhi - ylo + 1)
    lines = [10 ** (logslope * xs + yy + xlo) for yy in ys]
    for i in range(len(lines[:-1])):
        ax.fill_between(10 ** xs, lines[i], lines[i + 1], lw=0, color=colour, alpha= 1 - i / len(lines[:-1]))

def plot_hist(ax, data):
    label, fn, cut, c, ls, x, y = data
    x = 10**x
    if cut > 0: ax.plot(x, y, c=c, ls=ls, alpha=0.3)
    ax.plot(x[x > cut], y[x > cut], c=c, ls=ls, label=label)
