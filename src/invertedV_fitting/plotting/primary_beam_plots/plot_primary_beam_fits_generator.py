# --- primaryBeamFits_Plotting.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: get the data from the primary beam fits and plot the data WITHOUT regenerating all the fits again
import numpy as np

from src.invertedV_fitting.primary_beam_fit.primary_beam_fit_toggles import PrimaryBeamToggles
from src.invertedV_fitting.plotting.primary_beam_plots.plot_primary_beam_fits_toggles import PlotPrimaryBeamFitsToggles
from src.invertedV_fitting.primary_beam_fit.primary_beam_fit_classes import PrimaryBeamClasses


def plot_primary_beam_fits_generator():

    from src.invertedV_fitting.user_toggles.user_toggles import UserToggles
    from src.invertedV_fitting.primary_beam_fit.primary_beam_fit_classes import PrimaryBeamClasses
    import spaceToolsLib as stl
    from copy import deepcopy
    import matplotlib.pyplot as plt
    from tqdm import tqdm
    import os, shutil
    import datetime as dt
    from glob import glob


    ###########################
    # --- LOAD FIT THE DATA ---
    ###########################
    data_dict_beamFits = stl.loadDictFromFile(f'{UserToggles.run_folder_path}/primary_beam_fit.cdf')
    diffNFlux = data_dict_beamFits['diffNFlux'][0]
    counts = data_dict_beamFits['counts'][0]
    counts_std = data_dict_beamFits['counts_std'][0]
    epoch = deepcopy(data_dict_beamFits[f'{UserToggles.Epoch_key}'][0])
    energy = deepcopy(data_dict_beamFits['Energy'][0])
    pitch_angle = deepcopy(data_dict_beamFits['pitch_angle'][0])

    # check if output folder is there. If not, create it
    output_folder = f'{UserToggles.run_folder_path}/primary_beam_fits/'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # --- PLOTTING STYLES ---
    plt.style.use(r'/home/connor/PycharmProjects/invertedV_fitting/src/invertedV_fitting/plotting/primary_beam_plots/primary_beam_plots_style_sheet.mplstyle')
    my_cmap = stl.apl_rainbow_black0_cmap()
    my_cmap.set_bad(color=(0, 0, 0))
    cbarMin, cbarMax = 1E5, 1E7
    cbarTickLabelSize = 14
    cbar_Fontsize = 20

    ########################
    # --- MAKE THE PLOTS ---
    ########################
    # loop over time and fit the data at each point

    for tmeIdx in tqdm(range(len(epoch))):

        fig, ax = plt.subplots()
        fig.set_figheight(8)
        fig.set_figwidth(12)
        fig.suptitle(
            f'{epoch[tmeIdx]} UTC' +'\n'
        )

        # Raw Data - with error bars
        try:
            error = PrimaryBeamClasses().calc_jN_error(counts_val=counts[tmeIdx], energy_value=energy, pitch_angles=pitch_angle)
            yerr = [ np.clip(diffNFlux[tmeIdx]-error,min=0), np.clip(diffNFlux[tmeIdx]+error,min=0)]
            ax.plot(energy,diffNFlux[tmeIdx],zorder=0)
            # ax.errorbar(x=energy, y=diffNFlux[tmeIdx], yerr=yerr, capsize=4, color='tab:blue',zorder=0)
        except Exception as e:
            print(diffNFlux[tmeIdx]-error)

        ax.scatter(energy, diffNFlux[tmeIdx], zorder=0)

        # plot the phi0 guess
        phi0_guess = data_dict_beamFits['phi'][0][tmeIdx]
        ax.axvline(x=phi0_guess,color='red',linestyle='-')

        # plot the fit itself
        data_idx = int(data_dict_beamFits['find_fit_data_engy_idx'][0][tmeIdx])

        if data_idx != 0:
            # (a) get the fit data again
            xData_fit = energy[:data_idx+1]

            # (b) create a psuedo dataset
            fit_func, kwargs_dict = PrimaryBeamClasses().form_fit_params(phi0_guess)
            xData_plot = np.linspace(xData_fit.min(), xData_fit.max(),1000)
            params = [data_dict_beamFits['n'][0][tmeIdx],data_dict_beamFits['Te'][0][tmeIdx],data_dict_beamFits['phi'][0][tmeIdx]]

            if PrimaryBeamToggles.fit_dist == 'kappa':
                params +=  [data_dict_beamFits['kappa'][0][tmeIdx]]
            yData_plot = fit_func(xData_plot, *params)

            # (c) Plot the fitted data
            # legend_label = r'n$_{0}$ = ' + f'{round(params[0],2)}' +'cm$^{-3}$\n' + r'T$_{e}$ = ' + f'{round(params[1],2)} eV\n' + r'$\phi_{0}$ = ' + f'{round(params[2],2)} eV\n' + r'$\chi^{2}$ = ' + f'{round(data_dict_beamFits['chi2'][0][tmeIdx],3)}'
            legend_label = r'n$_{0}$ = ' + f'{round(params[0], 2)}' + 'cm$^{-3}$\n' + r'T$_{e}$ = ' + f'{round(params[1], 2)} eV\n' + r'$\phi_{0}$ = ' + f'{round(params[2], 2)} eV\n' + r'$\chi^{2}$ = '
            if PrimaryBeamToggles.fit_dist == 'kappa':
                legend_label += r'$\kappa$ = ' + f'{round(params[3],2)}\n'

            ax.plot(xData_plot, yData_plot, color='red', label=legend_label,zorder=1)


        # add the baseline counts level
        noise_level = PrimaryBeamClasses().calc_noise_level(counts_level= PlotPrimaryBeamFitsToggles.counts_level,
                                                                energy_value=energy,
                                                                pitch_angles=pitch_angle)
        ax.plot(energy, noise_level ,color='black', label=f'{PlotPrimaryBeamFitsToggles.counts_level}-count noise')

        # plot formatting
        ax.set_ylabel('[cm$^{-2}$-s$^{-1}$-str$^{-1}$-eV$^{-1}$]')
        ax.set_xlabel('Energy [eV]')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_ylim(1E3, 1E7)
        ax.legend()

        # Save the plot
        plt.tight_layout()
        fig.savefig(
            f'{output_folder}/fit_data_{tmeIdx}.png'
        )
        plt.close()

