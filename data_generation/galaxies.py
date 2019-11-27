import sys

import numpy as np

from astropy.cosmology import FlatLambdaCDM

from astropy.io import fits

if len(sys.argv) == 2:
    dataset = sys.argv[1]
else:
    dataset = "dwarfGal"

def get_data(dataset="SDSS"):
    # Cosmological values from Planck, as used by IllustrisTNG
    Omega_m = 0.3089
    h = 0.6774
    cosmo = FlatLambdaCDM(H0=h*100.0, Om0=Omega_m)

    if dataset == "SDSS":
        # Columns in table are tot_mass, tot_mass_width, z, z_err, objID, specObjID, petroR50_r, R_hlr (for 202842 objects)
        sdssdr7 = np.genfromtxt("../data/sdssdr7_mass_radius.csv", delimiter=',', skip_header=1)

        #SDSS data
        tot_mass = sdssdr7[:, 0] # log(M) with M in M_sun
        # tot_mass_width = sdssdr7[:, 1]
        z = sdssdr7[:, 2] # redshift
        # z_err = sdssdr7[:, 3]
        # petroR50_r = sdssdr7[:, 6]
        # R_hlr = sdssdr7[:, 7]
    elif dataset == "GAMA":
        gamaData = fits.open('../data/GAMA.fits')[1].data
        tot_mass = gamaData['logmstar'] # log(M) with M in M_sun
        z = gamaData['Z'] # redshift

        mass_limit = np.where((tot_mass > 5.9) * (z < 0.08))
        tot_mass = tot_mass[mass_limit]
        z = z[mass_limit]
    elif dataset == "dwarfGal":
        dg = fits.open('../data/dwarfGal.fits')[1].data
        tot_mass = np.log10(1e6 * dg["Mass"]) # log(M) with M in M_sun
        distance = dg["D_MW_"]/1e3 # distance in Mpc
        # radius = dg["R1"] # in arcmin

        val_mass = np.where(np.isfinite(tot_mass))

        tot_mass = tot_mass[val_mass]
        distance = distance[val_mass]
        # radius = radius[val_mass]

    if dataset in ["SDSS", "GAMA"]:
        distance = cosmo.comoving_distance(z).to("Mpc").value

    return (tot_mass, distance)

# Bin data

def get_dist_bins(distance, dataset="SDSS"):
    _, dist_bin_edges = np.histogram(distance, bins=50)
    # _, dist_bin_edges, _ = ax_dist.hist(distance, bins=50, color='k', histtype="step")

    # dist_bins = [(-np.inf, np.inf), (-np.inf, 100), (100, 200), (200, 300), (300, 400)]

    dist_bins = [(-np.inf, np.inf)]
    labels = ["All"]

    if dataset == "SDSS":
        dbedges = np.array([0, 50, 100, 150, 200, 250, 275, 300, 320, 350])
    elif dataset == "GAMA":
        dbedges = np.arange(0, 400, 50)
    elif dataset == "dwarfGal":
        dbedges = np.zeros(1)

    for di in range(dbedges.size - 1):
        dist_bins.append((dbedges[di], dbedges[di+1]))

        if dbedges[di] == 0:
            labels.append(r"$D<{:.0f} \, \mathrm{{Mpc}}$".format(dbedges[di+1]))
        elif dbedges[di+1] > np.max(distance):
            labels.append(r"$D>{:.0f} \, \mathrm{{Mpc}}$".format(dbedges[di]))
        else:
            labels.append(r"${:.0f} < D \leq {:.0f} \, \mathrm{{Mpc}}$".format(dbedges[di], dbedges[di+1]))
    
    return (dist_bin_edges, dist_bins, labels)

def get_mass_bins(dataset="SDSS"):
    if dataset in ["SDSS", "GAMA"]:
        # N_bins = int(1e4)
        # M_bin_edges = np.arange(np.round(min_mass, 1)-0.1, np.round(max_mass, 1)+0.1, 0.1)
        # M_bin_edges = np.log10(np.logspace(np.round(min_mass, 1)-0.1, np.round(max_mass, 1)+0.1, 1e4))
        M_bin_edges = np.linspace(5, np.log10(5e15), 200)
        # M_bin_edges_Ill = np.loadtxt("../data/illustris/mbins.txt")

        M_bin_centers = 0.5*(M_bin_edges[:-1] + M_bin_edges[1:])
    elif dataset == "dwarfGal":
        M_bin_edges = np.linspace(2.5, 9.5, 5)

        M_bin_centers = 0.5*(M_bin_edges[:-1] + M_bin_edges[1:])


    return (M_bin_edges, M_bin_centers)

def get_dndmdv(tot_mass, M_bin_edges, M_bin_centers, distance, dist_bin, dataset="SDSS"):
    dist_mask = (distance > dist_bin[0]) * (distance <= dist_bin[1])
    max_distance = np.max(distance[dist_mask])

    N, _ = np.histogram(tot_mass[dist_mask], bins=M_bin_edges)

    N_err = np.sqrt(N)

    dNdM = N/10**(M_bin_centers)
    dNdM_err = N_err/10**(M_bin_centers)

    if dataset == "SDSS":
        sky_fraction = 8032.0/(4.0*np.pi*(180.0/np.pi)**2) #8032.0 square degrees for SDSS
    elif dataset == "GAMA":
        sky_fraction = 296.158/(4.0*np.pi*(180.0/np.pi)**2)
    elif dataset == "dwarfGal":
        sky_fraction = 1.0
        # Use a volume of 3 Mpc
        max_distance = 3.0

    if np.isneginf(dist_bin[0]) and np.isposinf(dist_bin[1]):
        dV = 4.0/3.0 * np.pi * (max_distance * 1e6)**3
    else:
        assert not (np.isneginf(dist_bin[0]) or np.isposinf(dist_bin[1]))
        dV = 4.0/3.0 * np.pi * ((np.min([dist_bin[1], max_distance]) * 1e6)**3 - (dist_bin[0] * 1e6)**3)

    dNdV = N/(dV*sky_fraction)
    dNdV_err = N_err/(dV*sky_fraction)

    dNdMdV = dNdM/(dV*sky_fraction)
    dNdMdV_err = dNdM_err/(dV*sky_fraction)

    return (dNdM, dNdM_err, dNdV, dNdV_err, dNdMdV, dNdMdV_err)

tot_mass, distance = get_data(dataset=dataset)

# plt.hist(tot_mass, bins=50)

min_mass = np.min(tot_mass)
max_mass = np.max(tot_mass)
# print("Mass range:", min_mass, max_mass)

dist_bin_edges, dist_bins, labels = get_dist_bins(distance, dataset=dataset)

M_bin_edges, M_bin_centers = get_mass_bins(dataset=dataset)

if dataset == "dwarfGal":
    dNdM, dNdM_err, dNdV, dNdV_err, dNdMdV, dNdMdV_err = get_dndmdv(tot_mass, M_bin_edges, M_bin_centers, distance, dist_bin=(-np.inf, np.inf), dataset=dataset)
elif dataset in ["SDSS", "GAMA"]:
    dNdV = np.tile(-np.inf, M_bin_centers.size)
    dNdMdV = np.tile(-np.inf, M_bin_centers.size)

    for dbi, dist_bin in enumerate(dist_bins[1:]):
        dNdM, dNdM_err, dNdV, dNdV_err, dNdMdV_bin, dNdMdV_err = get_dndmdv(tot_mass, M_bin_edges, M_bin_centers, distance, dist_bin, dataset=dataset)

        dNdV = np.where(dNdV_err > dNdV, dNdV_err, dNdV)
        dNdMdV = np.where(dNdMdV_bin > dNdMdV, dNdMdV_bin, dNdMdV)

data = np.array([M_bin_centers, dNdV, dNdMdV])

footerText = "/Galaxies (SDSS)/'#2ecc71'/'-'/"
np.savetxt("../data/galaxies_obs_" + dataset + ".txt", data.T,
            fmt='%1.3e \t', header="M (M_s) \t dN/dV (pc^-3) \t dN/dMdV (M_s^-1 pc^-3)", footer=footerText)