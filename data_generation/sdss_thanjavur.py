""""
sdss_thanjavur.py

Quick script to generate SDSS data following Thanjavur+2016
"""
PHISTAR = 1.911e-3
MSTAR = 11.116
ALPHA = -1.145

import numpy as np

def schechter_fit(m, phistar=PHISTAR,
        mstar=MSTAR, alpha=ALPHA):
    return 1e-18*(np.log(10)
            *(phistar*(np.power(10, (m-mstar)*(alpha+1))))
            *np.exp(-np.power(10, (m-mstar))))

bins = np.vstack((np.arange(9, 12.04, 0.05),
    np.arange(9.05, 12.06, 0.05)))
bin_centres = np.mean(bins, axis=0)
bin_widths = np.diff(bins, axis=0)[0,:]

dndv = schechter_fit(m=bin_centres)
dndvdm = np.divide(dndv, bin_centres)

out = np.stack((bin_centres, bin_widths, dndv, dndvdm),
        axis=1)

np.savetxt('../data/sdss_thanjavur2016.txt', out,
        header='#bin_centre bin_width dndV dn(dVdlogM)')
