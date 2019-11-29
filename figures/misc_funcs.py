import numpy as np
import os
pwd = os.getcwd()

def load_data(fn):
    x, _, _, y = np.genfromtxt(pwd + '/../data/' + fn).T
    idxs = np.nonzero(y)
    x = x[idxs]
    y = y[idxs]
    return x, y
