#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 16:12:36 2019

@author: jsb
"""

import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

#%%
base = '/home/jsb/HackWeek/'
run = 'Auriga'

f = h5.File(base+run+'/snapshot_reduced_halo_6_063.hdf5')
dmpos = f['PartType1/Coordinates'][:]
r = np.sqrt(np.sum(dmpos**2, axis=1))
halodat = np.loadtxt(base+run+'/halos_0.0.ascii')
m200 = halodat[:,27]
m200 = m200[m200 < 1e9]

bins = np.arange(5, 16, 0.05)
bs = 427.66*1e3
w = 1/(4./3 * np.pi *(bs/2)**3)
h, edges = np.histogram(np.log10(m200), bins=bins, weights=np.full_like(m200, w))
edges = edges[:-1] + np.diff(edges)/2

totmass = sum(f['PartType1/Masses'][:])+sum(f['PartType0/Masses'][:])+sum(f['PartType4/Masses'][:])+sum(f['PartType5/Masses'][:])
totdens = (totmass / (4./3 * np.pi *(bs/(2*1e6))**3)) / (0.6777**2)

fig, ax = plt.subplots()
ax.set_xscale('log')
ax.set_yscale('log')
ax.step(10**edges, h/(0.05*totdens/27.75))

TNG100dat = np.loadtxt(base+'TNG100data/BigPlotData_TNG100_mvir.txt')
ax.step(10**TNG100dat[:,0], TNG100dat[:,3], label='TNG100')


massnumdens = np.vstack((edges, np.full_like(edges, 0.05), h/(totdens/27.75),  h/(0.05*totdens/27.75))).T
np.savetxt('BigPlotData_'+run+'_mvir.txt', massnumdens, header='[0] log10 Mass / Msun, [1] Bin width log10 Mass / Msun, [2] dN/dV (pc^-3), [3] dN/dVlogM (pc^-3 log10 Mass/Msun)')
