import numpy as np
import matplotlib.pyplot as plt
import os
pwd = os.getcwd()
plt.ion()

from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive
pl = NasaExoplanetArchive.get_confirmed_planets_table()
pl.keys()

jup_to_sun = 1.89813e27 / 1.989e30

m = pl['pl_bmassj'].value # [M_J]
mu = pl['pl_bmassjerr1'].value # upper bound
ml = pl['pl_bmassjerr2'].value # lower
pl['pl_bmassjlim'].value # whether the mass is M sin i

d = pl['st_dist'].value # [pc]
du = pl['st_disterr1'].value # upper
dl = pl['st_disterr2'].value # lower

#m *= jup_to_sun # [M_\odot]
#mu *= jup_to_sun
#ml *= jup_to_sun

bins = np.logspace(-3, np.log10(50), 100)
hist, bin_edges = np.histogram(m, bins)
hist, bin_edges
print(np.sum(hist), 'planets in bins.\n', len(m), 'planets in dataset.')

plt.figure()
plt.plot(m, '.')
plt.yscale('log')
plt.ylabel('Mass [M$_J$]')

plt.figure()
plt.plot(d, m, '.')
plt.yscale('log')
plt.xlabel('Distance [pc]')
plt.ylabel('Mass [M$_J$]')

plt.figure()
plt.hist(m * jup_to_sun, bins * jup_to_sun)
plt.xscale('log')
plt.xlabel(r'Mass [M$_\odot$]')
plt.ylabel('N')

dmdn = hist / (np.diff(bin_edges) * jup_to_sun)
cut_dist = 1000 # [pc]
dmdndv = dmdn / cut_dist**3


idxs = np.nonzero(dmdndv)
bins_new = bins[idxs]
dmdndv_new = dmdndv[idxs]

plt.figure()
plt.plot(bins_new * jup_to_sun, dmdndv_new, '-')
plt.plot(bins_new * jup_to_sun, dmdndv_new, 'r.')

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Mass [M$_\odot$]')
plt.ylabel(r'Number density $dN / (dM dV)$ [M$_\odot^{-1}$ pc$^{-3}]$')
plt.tight_layout()

np.savetxt(pwd + '../data/planets_obs.txt', np.array([bins_new * jup_to_sun, dmdndv_new]).T, header=r"Planets.   M [M_\odot] 	 dN / (dM dV) [M_\odot^-1 pc^-3]")

# 1) normalize bin counts by dividing by bin width
# 2) Need distance out to which most these planets come from (distance out to which sample is ~complete)
# 3) divide normalized bin counts by this distance^3 --> dN / (dM dV)
