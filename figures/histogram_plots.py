import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys
sys.path.insert(0,'../data_generation')
from plot_funcs import *

fig_savename = pwd + '/dn_dv_dlogm.png'

data = [
        ['Planets', 'planets_obs.txt', 0, '#e6194B', '-'],
        ['Transiting GK Planets', 'transitingPlanets_GK.txt', 0, '#2EC4B6', '-'],
        #['White dwarfs', 'WD_Number_Density.csv', 0, '#3cb44b', '-'],
        #['Neutron stars', 'NS_Number_Density.csv', 0, '#ffe119', '-'],
        ['Stars, obs.', 'MADE_dNdVdlogM_M.csv', 0, '#4363d8', '-'],
        ['Stars, IMF', 'imf_dNdVdlogM_M.csv', 0, '#f58231', '-'],
        ['Dwarf galaxies', 'galaxies_obs_dwarfGal.txt', 0, '#911eb4', '-'],
        ['Galaxies, SDSS', 'galaxies_obs_SDSS.txt', 1e7, '#42d4f4', '-'],
        ['SDSS, Thanjavur+2016', 'sdss_thanjavur2016.txt', 0, '#f032e6', '-'],
        ['Galaxies, GAMA', 'galaxies_obs_GAMA.txt', 1e7, '#bfef45', '-'],
        ['Galaxies, Illustris TNG300', 'illustris/BigPlotData_TNG300_mstar.txt', 2e7, '#fabebe', '-'],
        ['Galaxies, Illustris TNG100', 'illustris/BigPlotData_TNG100_mstar.txt', 2.3e6, '#469990', '-'],
        ['Galaxies, EAGLE100', 'eagle/BigPlotData_EAGLE_mstar.txt', 1e7, '#e6beff', '-'],
        ['Galaxies, EAGLE25', 'eagle/BigPlotData_EAGLE25_mstar.txt', 1e6, '#9A6324', '-'],
        ['Auriga, M$_{vir}$', 'BigPlotData_Auriga_mvir.txt', 0, '#fffac8', '-'],
        ['TNG300, M$_{vir}$', 'illustris/BigPlotData_TNG300_mvir.txt', 1.2e10, '#800000', '-'],
        ['TNG100, M$_{vir}$', 'illustris/BigPlotData_TNG100_mvir.txt', 1.5e9, '#aaffc3', '-'],
        ['EAGLE100, M$_{vir}$', 'eagle/BigPlotData_EAGLE_mvir.txt', 3e8, '#808000', '-'],
        ['EAGLE25, M$_{vir}$', 'eagle/BigPlotData_EAGLE25_mvir.txt', 3e7, '#ffd8b1', '-'],
        ['Magneticum, M$_{vir}$', 'BigPlotData_Magneticum_mvir.txt', 0, '#000075', '-']
       ]

xlab, ylab = [r'Mass [M$_\odot]$', r'Number density, dN / dV d log$_{10}$(M) [pc$^{-3}]$']
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
from theoretical_constraints import M_lim, baryons, cdm, matter, nmin, ncollapse, m_nps, nps
ax1.plot(M_lim, baryons, c='#3498DB', ls='--', label='Baryons')
ax1.plot(M_lim, cdm,  c='#2ECC71', ls='--', label='Dark matter')
ax1.plot(M_lim, matter, c='#8E44AD', ls=':', label='Baryons + DM')
ax1.axhline(nmin, ls='--', c='#F39C12', label='Minimum')
ax1.plot(M_lim, ncollapse / M_lim, ls='--', c='k', label='Self-collapse')
ax1.plot(m_nps[:-1], nps, ls='--', c='g', label='Press–Schechter')


## observations and simulations
for i in range(len(data)):
    print('loading', data[i])
    x, y = load_data(data[i][1])
    data[i].extend([x,y])
    plot_hist(ax1, data[i])

handles, labels = ax1.get_legend_handles_labels()
print('labels',labels)
'''
reorder = [0, 1, 2, 3, 4, 5]
handles = [handles[i] for i in reorder]
labels = [labels[i] for i in reorder]
'''
plt.legend(handles, labels, loc=[.4,.7], bbox_transform=ax1.transAxes, ncol=3, fontsize=6)
plt.savefig(fig_savename)
