import numpy as np

from astropy.cosmology import FlatLambdaCDM

from astropy.io import fits

def get_data(dataset="SDSS"):
    # Cosmological values from Planck, as used by IllustrisTNG
    Omega_m = 0.3089
    h = 0.6774
    cosmo = FlatLambdaCDM(H0=h*100.0, Om0=Omega_m)

    if dataset == "SDSS":
        # Columns in table are tot_mass, tot_mass_width, z, z_err, objID, specObjID, petroR50_r, R_hlr (for 202842 objects)
        sdssdr7 = np.genfromtxt("../data/sdssdr7_mass_radius.csv", delimiter=',', skip_header=1)

        #SDSS data
        tot_mass = sdssdr7[:, 0]
        # tot_mass_width = sdssdr7[:, 1]
        z = sdssdr7[:, 2]
        # z_err = sdssdr7[:, 3]
        # petroR50_r = sdssdr7[:, 6]
        # R_hlr = sdssdr7[:, 7]
    elif dataset == "GAMA":
        gamaData = fits.open('../data/GAMA.fits')[1].data
        tot_mass = gamaData['logmstar']
        z = gamaData['Z']

        mass_limit = np.where((tot_mass > 5.9) * (z < 0.08))
        tot_mass = tot_mass[mass_limit]
        z = z[mass_limit]

    distance = cosmo.comoving_distance(z).to("Mpc").value

    return (tot_mass, distance)

tot_mass, distance = get_data(dataset="SDSS")

# plt.hist(tot_mass, bins=50)

min_mass = np.min(tot_mass)
max_mass = np.max(tot_mass)
# print("Mass range:", min_mass, max_mass)

# Bin data

def get_dist_bins(distance):
    _, dist_bin_edges = np.histogram(distance, bins=50)
    # _, dist_bin_edges, _ = ax_dist.hist(distance, bins=50, color='k', histtype="step")

    # dist_bins = [(-np.inf, np.inf), (-np.inf, 100), (100, 200), (200, 300), (300, 400)]

    dist_bins = [(-np.inf, np.inf)]
    labels = ["All"]

    # dbedges = np.arange(0, 400, 50)
    dbedges = np.array([0, 50, 100, 150, 200, 250, 275, 300, 320, 350])
    n_dist_bins = dbedges.size - 1

    for di in range(n_dist_bins):
        dist_bins.append((dbedges[di], dbedges[di+1]))

        if dbedges[di] == 0:
            labels.append(r"$D<{:.0f} \, \mathrm{{Mpc}}$".format(dbedges[di+1]))
        elif dbedges[di+1] > np.max(distance):
            labels.append(r"$D>{:.0f} \, \mathrm{{Mpc}}$".format(dbedges[di]))
        else:
            labels.append(r"${:.0f} < D \leq {:.0f} \, \mathrm{{Mpc}}$".format(dbedges[di], dbedges[di+1]))
    
    return (dist_bin_edges, dist_bins, labels)

dist_bin_edges, dist_bins, labels = get_dist_bins(distance)

def get_mass_bins():
    # N_bins = int(1e4)
    # M_bin_edges = np.arange(np.round(min_mass, 1)-0.1, np.round(max_mass, 1)+0.1, 0.1)
    # M_bin_edges = np.log10(np.logspace(np.round(min_mass, 1)-0.1, np.round(max_mass, 1)+0.1, 1e4))
    M_bin_edges = np.linspace(5, np.log10(5e15), 200)
    # M_bin_edges_Ill = np.loadtxt("../data/illustris/mbins.txt")

    M_bin_centers = 0.5*(M_bin_edges[:-1] + M_bin_edges[1:])

    return (M_bin_edges, M_bin_centers)

M_bin_edges, M_bin_centers = get_mass_bins()

def get_dndmdv(tot_mass, M_bin_edges, dist_bin):
    dist_mask = (distance > dist_bin[0]) * (distance <= dist_bin[1])
    max_distance = np.max(distance[dist_mask])

    N, _ = np.histogram(tot_mass[dist_mask], bins=M_bin_edges)

    N_err = np.sqrt(N)

    dNdM = N/10**(M_bin_centers)
    dNdM_err = N_err/10**(M_bin_centers)

    sky_fraction = 8032.0/(4.0*np.pi*(180.0/np.pi)**2) #8032.0 for SDSS

    if np.isneginf(dist_bin[0]) and np.isposinf(dist_bin[1]):
        dV = 4.0/3.0 * np.pi * (max_distance * 1e6)**3
    else:
        assert not (np.isneginf(dist_bin[0]) or np.isposinf(dist_bin[1]))
        dV = 4.0/3.0 * np.pi * ((np.min([dist_bin[1], max_distance]) * 1e6)**3 - (dist_bin[0] * 1e6)**3)

    dNdMdV = dNdM/(dV*sky_fraction)
    dNdMdV_err = dNdM_err/(dV*sky_fraction)

    return (dNdM, dNdM_err, dNdMdV, dNdMdV_err)

# dNdM, dNdM_err, dNdMdV, dNdMdV_err = get_dndmdv(tot_mass, M_bin_edges, dist_bin=(-np.inf, np.inf))

dNdMdV = np.tile(-np.inf, M_bin_centers.size)

for dbi, dist_bin in enumerate(dist_bins[1:]):
    dNdM, dNdM_err, dNdMdV_bin, dNdMdV_err = get_dndmdv(tot_mass, M_bin_edges, dist_bin)

    dNdMdV = np.where(dNdMdV_bin > dNdMdV, dNdMdV_bin, dNdMdV)

data = np.array([M_bin_centers, dNdMdV])

footerText = "/Galaxies (SDSS)/'#2ecc71'/'-'/"
np.savetxt("../data/galaxies_obs.txt", data.T, fmt='%1.3e \t', header="M (M_s) \t dN/dMdV (M_s^-1 pc^-3)", footer=footerText)