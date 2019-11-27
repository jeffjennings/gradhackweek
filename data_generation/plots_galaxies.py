import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("ticks")

from galaxies import get_data, get_dist_bins, get_mass_bins, get_dndmdv

dataset = "GAMA"

tot_mass, distance = get_data(dataset=dataset)
dist_bin_edges, dist_bins, labels = get_dist_bins(distance, dataset=dataset)
M_bin_edges, M_bin_centers = get_mass_bins()

fig_dist = plt.figure()
ax_dist = fig_dist.add_subplot(111)

fig_mass = plt.figure()
ax_mass = fig_mass.add_subplot(111)

colors = sns.color_palette("Set2", len(dist_bins))

ax_dist.hist([distance[(distance > dist_bin[0]) * (distance <= dist_bin[1])] for dist_bin in dist_bins[1:]],
                bins=dist_bin_edges, histtype="step", stacked=True, color=colors[1:], label=labels[1:])

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

fig_dndmdv = plt.figure()
ax_dndmdv = fig_dndmdv.add_subplot(111)

for dbi, dist_bin in enumerate(dist_bins):
    dNdM, dNdM_err, dNdMdV, dNdMdV_err = get_dndmdv(tot_mass, M_bin_edges, M_bin_centers, distance, dist_bin, dataset=dataset)

    ax_dndm.errorbar(M_bin_centers, dNdM, yerr=dNdM_err, color=colors[dbi], label=labels[dbi])
    ax_dndmdv.errorbar(M_bin_centers, dNdMdV, yerr=dNdMdV_err, color=colors[dbi], label=labels[dbi])

ax_dndm.legend(loc="upper right")
ax_dndm.set_yscale("log")

ax_dndm.set_xlabel(r"Mass ($\mathrm{M_\odot}$)")
ax_dndm.set_ylabel(r"Mass density ($\mathrm{M_\odot^{-1}}$)")

ax_dndmdv.legend(loc="upper right")
ax_dndmdv.set_yscale("log")

ax_dndmdv.set_xlabel(r"Mass ($\mathrm{M_\odot}$)")
ax_dndmdv.set_ylabel(r"Number density ($\mathrm{M_\odot^{-1} pc^{-3}}$)")

fig_dndm.savefig("../figures/Mass_density_" + dataset + ".pdf", rasterized=True, bbox_inches="tight", pad_inches=0.5)
fig_dndmdv.savefig("../figures/Number_density_" + dataset + ".pdf", rasterized=True, bbox_inches="tight", pad_inches=0.5)
# plt.show()