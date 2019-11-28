import sys

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("ticks")

from galaxies import get_data, get_dist_bins, get_mass_bins, get_dndmdv

if __name__ == '__main__':
    if len(sys.argv) > 1:
        datasets = [*sys.argv[1:]]
    else:
        datasets = ["SDSS", "GAMA", "dwarfGal"]
else:
    datasets = []

for dataset in datasets:
    tot_mass, distance, z = get_data(dataset=dataset)
    dist_bin_edges, dist_bins, z_bins, labels = get_dist_bins(distance, z=z, dataset=dataset)
    M_bin_edges, M_bin_centers = get_mass_bins(dataset=dataset)

    fig_dist = plt.figure()
    ax_dist = fig_dist.add_subplot(111)

    fig_mass = plt.figure()
    ax_mass = fig_mass.add_subplot(111)

    colors = sns.color_palette("Set2", len(dist_bins))

    if dataset in ["SDSS", "GAMA"]:
        ax_dist.hist([distance[(distance > dist_bin[0]) * (distance <= dist_bin[1])] for dist_bin in dist_bins[1:]],
                        bins=dist_bin_edges, histtype="step", stacked=True, color=colors[1:], label=labels[1:])
    elif dataset == "dwarfGal":
        ax_dist.hist(distance, bins=dist_bin_edges, histtype="step", stacked=False, color=colors[0], label=labels[0])

    for dbi, dist_bin in enumerate(dist_bins):
        dist_mask = (distance > dist_bin[0]) * (distance <= dist_bin[1])

        if dbi > 0:
            ax_dist.annotate(s=r"$N={:d}$".format(np.sum(dist_mask)), xy=(np.mean(dist_bin), 0), xytext=(0, 4+8*(1-(-1)**(dbi-1))),
                            xycoords=matplotlib.transforms.blended_transform_factory(ax_dist.transData, ax_dist.transAxes),
                            textcoords="offset points", color=colors[dbi], verticalalignment="bottom", horizontalalignment="center")

        N, bin_edges = np.histogram(tot_mass[dist_mask], bins=50)
        bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])

        ax_mass.errorbar(bin_centers, N, yerr=np.sqrt(N), color=colors[dbi], label=labels[dbi])

    ax_dist.legend(loc="upper left")

    ax_dist.set_xlabel(r"Distance ($\mathrm{Mpc}$)")
    ax_dist.set_ylabel(r"Number")

    ax_mass.legend(loc="upper left")

    ax_mass.set_xlabel(r"Mass ($\mathrm{M_\odot}$)")
    ax_mass.set_ylabel(r"Number")

    fig_dist.savefig("../figures/Distance_distribution_" + dataset + ".pdf", rasterized=True, bbox_inches="tight", pad_inches=0.5)
    fig_mass.savefig("../figures/Mass_distribution_" + dataset + ".pdf", rasterized=True, bbox_inches="tight", pad_inches=0.5)
    # plt.show()

    # Plot results

    fig_dndm = plt.figure()
    ax_dndm = fig_dndm.add_subplot(111)

    fig_dndv = plt.figure()
    ax_dndv = fig_dndv.add_subplot(111)

    fig_dndmdv = plt.figure()
    ax_dndmdv = fig_dndmdv.add_subplot(111)

    if dataset in ["SDSS", "GAMA"]:
        dNdV_final = np.tile(-np.inf, M_bin_centers.size)
        dNdMdV_final = np.tile(-np.inf, M_bin_centers.size)

    for dbi, dist_bin in enumerate(dist_bins):
        dNdM, dNdM_err, dNdV, dNdV_err, dNdMdV, dNdMdV_err = get_dndmdv(tot_mass, M_bin_edges, distance,
                                                                    dist_bin, z_bins[dbi], z=z, dataset=dataset)

        ax_dndm.errorbar(M_bin_centers, dNdM, linewidth=1.0, yerr=dNdM_err, color=colors[dbi], label=labels[dbi])
        ax_dndv.errorbar(M_bin_centers, dNdV, linewidth=1.0, yerr=dNdV_err, color=colors[dbi], label=labels[dbi])
        ax_dndmdv.errorbar(M_bin_centers, dNdMdV, linewidth=1.0, yerr=dNdMdV_err, color=colors[dbi], label=labels[dbi])

        if dataset in ["SDSS", "GAMA"]:
            dNdV_final = np.where(dNdV > dNdV_final, dNdV, dNdV_final)
            dNdMdV_final = np.where(dNdMdV > dNdMdV_final, dNdMdV, dNdMdV_final)

    if dataset in ["SDSS", "GAMA"]:
        ax_dndv.errorbar(M_bin_centers, dNdV_final, color='k', label="Final")
        ax_dndmdv.errorbar(M_bin_centers, dNdMdV_final, color='k', label="Final")

    ax_dndm.legend(loc="upper right")
    ax_dndm.set_yscale("log")

    ax_dndm.set_xlabel(r"Mass ($\mathrm{M_\odot}$)")
    ax_dndm.set_ylabel(r"Mass density ($\log \mathrm{M_\odot}^{-1}$)")

    ax_dndv.legend(loc="upper right")
    ax_dndv.set_yscale("log")

    ax_dndv.set_xlabel(r"Mass ($\mathrm{M_\odot}$)")
    ax_dndv.set_ylabel(r"Number density ($\mathrm{pc^{-3}}$)")

    ax_dndmdv.legend(loc="upper right")
    ax_dndmdv.set_yscale("log")

    ax_dndmdv.set_xlabel(r"Mass ($\mathrm{M_\odot}$)")
    ax_dndmdv.set_ylabel(r"Number density ($\log \mathrm{M_\odot}^{-1} \mathrm{pc^{-3}}$)")

    fig_dndm.savefig("../figures/Mass_density_" + dataset + ".pdf", rasterized=True, bbox_inches="tight", pad_inches=0.5)
    fig_dndv.savefig("../figures/Number_density_" + dataset + ".pdf", rasterized=True, bbox_inches="tight", pad_inches=0.5)
    fig_dndmdv.savefig("../figures/Number_mass_density_" + dataset + ".pdf", rasterized=True, bbox_inches="tight", pad_inches=0.5)
    # plt.show()