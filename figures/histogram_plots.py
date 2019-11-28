import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
pwd = os.getcwd()

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

gs = GridSpec(1, 1, bottom=.12, top=.95, left=.1, right=.98, hspace=0)
smallplot = True
if smallplot: fig = plt.figure(figsize=(10, 4))
else: fig = plt.figure(figsize=(11.69, 8.27))
ax1 = fig.add_subplot(gs[0])

xlo, xhi = 1e-6, 1e16
ylo, yhi = 1e-33, 1e0
M_lim = np.array([xlo, xhi])

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(xlo, xhi)
ax1.set_ylim(ylo, yhi)
ax1.set_xlabel(r'Mass [M$_\odot]$')
ax1.set_ylabel(r'Number density, $dN / dV$ [pc$^{-3}]$')


## background M^{-1} contours
plotBackground(ax1, xhi*1e12, ylo, '#011627')


## Planck cosmology limits
Omh2 = 0.143
Och2 = 0.12
Obh2 = 0.0237
H0 = 67.27
Ostar = 0.003

mass_density_gcm3 = 9.9e-30
mass_density_Mspc3 = mass_density_gcm3 * (100*3e16)**3 / (1000*2e30)

matter = (mass_density_Mspc3 * Omh2 / (H0*0.01)**2) / M_lim
baryons = (mass_density_Mspc3 * Obh2 / (H0*0.01)**2) / M_lim
cdm = (mass_density_Mspc3 * Och2 / (H0*0.01)**2) / M_lim

ax1.plot(M_lim, baryons, c='#3498DB', ls='--', label='Baryons')
ax1.plot(M_lim, cdm,  c='#2ECC71', ls='--', label='Dark matter')
ax1.plot(M_lim, matter, c='#8E44AD', ls=':', label='Baryons + DM')


## additional theoretical limits
nmin = 1e-32
ns = 2.6e38 / M_lim**3
nobs = 6.0953e-81 * M_lim**-1.5
ncollapse = 1.6e-7 / M_lim**2
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
