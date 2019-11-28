import numpy as np

## Planck cosmology limits
Omh2 = 0.143
Och2 = 0.12
Obh2 = 0.0237
H0 = 67.27
Ostar = 0.003

M_lim = np.array([1e-40, 1e25])

mass_density_gcm3 = 9.9e-30
mass_density_Mspc3 = mass_density_gcm3 * (100*3e16)**3 / (1000*2e30)

matter = (mass_density_Mspc3 * Omh2 / (H0*0.01)**2) / M_lim
baryons = (mass_density_Mspc3 * Obh2 / (H0*0.01)**2) / M_lim
cdm = (mass_density_Mspc3 * Och2 / (H0*0.01)**2) / M_lim


## additional theoretical limits
nmin = 1e-32
ns = 2.6e38 / M_lim**3
nobs = 6.0953e-81 * M_lim**-1.5
ncollapse = 1.6e-7 / M_lim**2
