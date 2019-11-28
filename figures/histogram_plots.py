import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
pwd = os.getcwd()
import sys
sys.path.insert(0,'../data_generation')

from plot_funcs import *

fig_savename = pwd + '/dn_dv_dlogm.png'

fns = ['planets_obs.txt', 'WD_Number_Density.csv', 'NS_Number_Density.csv', \
       'galaxies_obs_dwarfGal.txt', 'galaxies_obs_SDSS.txt', 'galaxies_obs_GAMA.txt', \
       'dNbydMdV_Mstar_TNG300.txt', 'dNbydMdV_Mstar_TNG100.txt', 'dNbydMdV_Mstar_EAGLE', \
       'dNbydMdV_Mstar_EAGLE25.txt', ]
fns = [pwd + '/../data/' + x for x in fns]

cuts = [0, 0, 0, 0, 1e7, 1e7, 2e7, 2.3e6, 1e7, 1e6] # completeness cuts in x-range

labs = ['Planets', 'White dwarfs', 'Neutron stars', 'Dwarf galaxies', 'Galaxies', \
        'SDSS', 'Galaxies, GAMA', 'Galaxis, Illustris TNG300', 'Galaxis, Illustris TNG100', \
        'Galaxies, EAGLE100', 'Galaxies, EAGLE25']

cs = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', \
      '#f032e6', '#bfef45', '#fabebe', '#469990', '#e6beff', '#9A6324', '#fffac8', \
      '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9']

lss = np.repeat('-', len(fns))

xlab, ylab = [r'Mass [M$_\odot]$', r'Number density, $dN / dV$ [pc$^{-3}]$']
xlo, xhi = 1e-6, 1e16
ylo, yhi = 1e-33, 1e0


gs = GridSpec(1, 1, bottom=.12, top=.95, left=.1, right=.98, hspace=0)
smallplot = True
if smallplot: fig = plt.figure(figsize=(10, 4))
else: fig = plt.figure(figsize=(11.69, 8.27))
ax1 = fig.add_subplot(gs[0])

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(xlo, xhi)
ax1.set_ylim(ylo, yhi)
ax1.set_xlabel(xlab)
ax1.set_ylabel(ylab)


## background M^{-1} contours
plotBackground(ax1, xhi*1e12, ylo, '#011627')


## theoretical limits
from theoretical_constraints import M_lim, baryons, cdm, matter, nmin, ncollapse
ax1.plot(M_lim, baryons, c='#3498DB', ls='--', label='Baryons')
ax1.plot(M_lim, cdm,  c='#2ECC71', ls='--', label='Dark matter')
ax1.plot(M_lim, matter, c='#8E44AD', ls=':', label='Baryons + DM')
ax1.axhline(nmin, ls='--', c='#F39C12', label='Minimum')
ax1.plot(M_lim, ncollapse / M_lim, ls='--', c='k', label='Self-collapse')


## observations and simulations
plot_all(ax1, fns, cuts, cs, lss, labs)


handles, labels = ax1.get_legend_handles_labels()
print('labels',labels)
'''
reorder = [0, 1, 2, 3, 4, 5]
handles = [handles[i] for i in reorder]
labels = [labels[i] for i in reorder]
'''
plt.legend(handles, labels, loc=[.4,.7], bbox_transform=ax1.transAxes, ncol=3, fontsize=6)
plt.savefig(fig_savename)
