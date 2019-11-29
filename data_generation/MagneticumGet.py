#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 10:47:10 2019

@author: jsb
"""
import numpy as np
import matplotlib.pyplot as plt
base = '/home/jsb/HackWeek/'
run = 'Magneticum'
magnet = np.loadtxt(base+'/Magneticum/magneticum_range.csv', delimiter=',')
logmag = np.log(magnet)

x = magnet[:,0]
trydat = ((magnet[:,1]*0.23))/(1e18)

massnumdens = np.vstack((np.log10(x), np.full_like(x, 0.1), trydat, trydat/0.1)).T
np.savetxt('BigPlotData_'+run+'_mvir.txt', massnumdens, header='[0] log10 Mass / Msun, [1] Bin width log10 Mass / Msun, [2] dN/dV (pc^-3), [3] dN/dVlogM (pc^-3 log10 Mass/Msun)')

fig, ax = plt.subplots()
ax.loglog(x, trydat)