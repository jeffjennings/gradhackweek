import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys
sys.path.insert(0,'../data_generation')
from plot_funcs import plot_background, plot_hist
from misc_funcs import load_data
from datasets import data_to_plot

#dataset = 'dn_dv_dlogm'
dataset = 'dn_dv_dlogr'

data, fig_savename, xlab, ylab, xlo, xhi, ylo, yhi = data_to_plot(dataset)

gs = GridSpec(1, 1, bottom=.12, top=.95, left=.1, right=.98, hspace=0)
smallplot = True
if smallplot: fig = plt.figure(figsize=(10, 4))
else: fig = plt.figure(figsize=(11.69, 8.27))
ax1 = fig.add_subplot(gs[0])

## background contours
logslope = -1 # loglog slope of contours
#plot_background(ax1, xlo, xhi, ylo + 20, yhi, logslope, 'k')


## theoretical limits
if dataset == 'dn_dv_dlogm':
    from theoretical_constraints import M_lim, baryons, cdm, matter, nmin, ncollapse, m_nps, nps
    ax1.plot(M_lim, baryons, c='#3498DB', ls='--', label='Baryons')
    ax1.plot(M_lim, cdm,  c='#2ECC71', ls='--', label='Dark matter')
    ax1.plot(M_lim, matter, c='#8E44AD', ls=':', label='Baryons + DM')
    ax1.plot(m_nps[:-1], nps, ls=':', c='g', label='Press–Schechter')
    #ax1.plot(M_lim, ncollapse / M_lim, ls='--', c='k', alpha=.3, label='Self-collapse')
    ax1.axhline(nmin, ls='--', c='#F39C12', label='1 per universe')


## observations and simulations
for i in range(len(data)):
    print('loading, plotting: ', dataset, '', data[i][0])
    x, y = load_data(data[i][1])
    data[i].extend([x,y])
    print('y',y)
    plot_hist(ax1, data[i])

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(10 ** xlo, 10 ** xhi)
ax1.set_ylim(10 ** ylo, 10 ** yhi)
ax1.set_xlabel(xlab)
ax1.set_ylabel(ylab)

handles, labels = ax1.get_legend_handles_labels()
'''
print('labels',labels)
reorder = [0, 1, 2, 3, 4, 5]
handles = [handles[i] for i in reorder]
labels = [labels[i] for i in reorder]
'''
plt.legend(handles, labels, loc=[.45,.7], bbox_transform=ax1.transAxes, ncol=3, fontsize=6)
plt.savefig(fig_savename, dpi=300)
