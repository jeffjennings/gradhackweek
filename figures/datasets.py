import os
pwd = os.getcwd()

def data_to_plot(set):
    """
    The form of each entry in data is:
    [Legend label, datafile, x-value below which data is plotted with low alpha, line color, linestyle]
    """
    if set == 'dn_dv_dlogm':
        data = [
                ['Planets', 'planets_obs.txt', 0, '#e6194B', '-'],
                ['Transiting GK Planets', 'transitingPlanets_GK.txt', 0, '#2EC4B6', '-'],
                #['White dwarfs', 'WD_Number_Density.csv', 0, '#3cb44b', '-'],
                #['Neutron stars', 'NS_Number_Density.csv', 0, '#ffe119', '-'],
                ['Stars, obs.', 'MADE_dNdVdlogM_M.csv', 0, '#4363d8', '-'],
                ['Stars, IMF', 'imf_dNdVdlogM_M.csv', 0, '#f58231', '-'],
                ['Dwarf galaxies', 'galaxies_obs_mass_dwarfGal.txt', 0, '#911eb4', '-'],
                ['Galaxies, SDSS', 'galaxies_obs_mass_SDSS.txt', 1e7, '#42d4f4', '-'],
                ['SDSS, Thanjavur+2016', 'sdss_thanjavur2016.txt', 0, '#f032e6', '-'],
                ['Galaxies, GAMA', 'galaxies_obs_mass_GAMA.txt', 1e7, '#bfef45', '-'],
                ['Galaxies, Illustris TNG300', 'illustris/BigPlotData_TNG300_mstar.txt', 2e7, '#fabebe', '-'],
                ['Galaxies, Illustris TNG100', 'illustris/BigPlotData_TNG100_mstar.txt', 2.3e6, '#469990', '-'],
                ['Galaxies, EAGLE100', 'eagle/BigPlotData_EAGLE_mstar.txt', 1e7, '#e6beff', '-'],
                ['Galaxies, EAGLE25', 'eagle/BigPlotData_EAGLE25_mstar.txt', 1e6, '#9A6324', '-'],
                ['Auriga, M$_{vir}$', 'BigPlotData_Auriga_mvir.txt', 0, '#fffac8', '-'],
                ['TNG300, M$_{vir}$', 'illustris/BigPlotData_TNG300_mvir.txt', 1.2e10, '#800000', '-'],
                ['TNG100, M$_{vir}$', 'illustris/BigPlotData_TNG100_mvir.txt', 1.5e9, '#aaffc3', '-'],

                ['TNG300, M$_{BH}', 'illustris/BigPlotData_TNG300_mbh.txt', 1e7, 'k', '-'],
                ['TNG100, M$_{BH}', 'illustris/BigPlotData_TNG100_mbh.txt', 1e7, 'k', '-'],

                ['EAGLE100, M$_{vir}$', 'eagle/BigPlotData_EAGLE_mvir.txt', 3e8, '#808000', '-'],
                ['EAGLE25, M$_{vir}$', 'eagle/BigPlotData_EAGLE25_mvir.txt', 3e7, '#ffd8b1', '-'],
                ['Magneticum, M$_{vir}$', 'BigPlotData_Magneticum_mvir.txt', 0, '#000075', '-']
               ]

        fig_savename = pwd + '/dn_dv_dlogm.png'
        xlab, ylab = [r'Mass [M$_\odot]$', r'Number density, dN / dV d log$_{10}$(M) [pc$^{-3}]$']
        xlo, xhi = -6, 16 # figure x-axis bounds in log10(x)
        ylo, yhi = -33, 0

    if set == 'dn_dv_dlogr':
        data = [
                #['Planets', 'planets_obs.txt', 0, '#e6194B', '-'],
                #['Transiting GK Planets', 'transitingPlanets_GK.txt', 0, '#2EC4B6', '-'],
                #['White dwarfs', 'WD_Number_Density.csv', 0, '#3cb44b', '-'],
                #['Neutron stars', 'NS_Number_Density.csv', 0, '#ffe119', '-'],
                #['Stars, obs.', 'MADE_dNdVdlogM_M.csv', 0, '#4363d8', '-'],
                #['Stars, IMF', 'imf_dNdVdlogM_M.csv', 0, '#f58231', '-'],
                ['Dwarf galaxies', 'galaxies_obs_radius_dwarfGal.txt', 0, '#911eb4', '-'],
                ['Galaxies, SDSS', 'galaxies_obs_radius_SDSS.txt', 0, '#42d4f4', '-'],
                #['SDSS, Thanjavur+2016', 'sdss_thanjavur2016.txt', 0, '#f032e6', '-'],
                ['Galaxies, GAMA', 'galaxies_obs_radius_GAMA.txt', 0, '#bfef45', '-'],
                ['Galaxies, Illustris TNG300', 'illustris/BigPlotRadiusData_TNG300_rhalfmassstar.txt', 0, '#fabebe', '-'],
                ['Galaxies, Illustris TNG100', 'illustris/BigPlotRadiusData_TNG100_rhalfmassstar.txt', 0, '#469990', '-'],
                ['Galaxies, EAGLE100', 'eagle/BigPlotRadiusData_EAGLE_halfstarrad.txt', 0, '#e6beff', '-'],
                #['Galaxies, EAGLE25', 'eagle/BigPlotData_EAGLE25_mstar.txt', 1e6, '#9A6324', '-'],
                #['Auriga, M$_{vir}$', 'BigPlotData_Auriga_mvir.txt', 0, '#fffac8', '-'],
                ['TNG300, R$_{vir}$', 'illustris/BigPlotRadiusData_TNG300_rvir.txt', 0, '#800000', '-'],
                ['TNG100, R$_{vir}$', 'illustris/BigPlotRadiusData_TNG100_rvir.txt', 0, '#aaffc3', '-'],
                ['EAGLE100, R$_{vir}$', 'eagle/BigPlotRadiusData_EAGLE_rvir.txt', 3e8, '#808000', '-'],
                #['EAGLE25, M$_{vir}$', 'eagle/BigPlotData_EAGLE25_mvir.txt', 3e7, '#ffd8b1', '-'],
                #['Magneticum, M$_{vir}$', 'BigPlotData_Magneticum_mvir.txt', 0, '#000075', '-']
               ]

        fig_savename = pwd + '/dn_dv_dlogr.png'
        xlab, ylab = [r'Radius [pc]', r'Number density, dN / dV d log$_{10}$(R) [pc$^{-3}]$']
        xlo, xhi = -10, 10 # figure x-axis bounds in log10(x)
        ylo, yhi = -25, -10

    if set == 'dn_dv_dlogrho':
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

        fig_savename = pwd + '/dn_dv_dlogrho.png'
        xlab, ylab = [r'Object mass density [$\rho_\odot$]$', r'Number density, dN / dV d log$_{10}(\rho) [pc$^{-3}]$']
        xlo, xhi = -20, 20 # figure x-axis bounds in log10(x)
        ylo, yhi = -20, 20

    if set == 'm_vs_r':
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

        fig_savename = pwd + '/m_vs_r.png'
        xlab, ylab = [r'Radius [R$_\odot$]', r'Mass [M$_\odot$]']
        xlo, xhi = -20, 20 # figure x-axis bounds in log10(x)
        ylo, yhi = -20, 20

    return data, fig_savename, xlab, ylab, xlo, xhi, ylo, yhi
