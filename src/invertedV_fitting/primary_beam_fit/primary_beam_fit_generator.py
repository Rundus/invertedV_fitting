# --- primary_beam_fit_generator.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: using the method outline in Kaeppler's thesis, we can fit inverted-V distributions
# to get estimate the magnetospheric temperature, density and electrostatic potential that accelerated
# our particles

from timebudget import timebudget

@timebudget
def primary_beam_fit_generator():

    # --- imports ---
    import spaceToolsLib as stl
    import numpy as np
    from tqdm import tqdm
    from glob import glob
    from src.invertedV_fitting.user_toggles import FitDataToggles,PrimaryBeamFitToggles, FileToggles
    from src.invertedV_fitting.primary_beam_fit.primary_beam_fit_classes import PrimaryBeamClasses
    from scipy.optimize import curve_fit
    import datetime as dt

    ###############################
    # --- LOADING THE FLUX DATA ---
    ###############################
    data_dict = stl.loadDictFromFile(glob(f'{FileToggles.RUN_PATH}/fit_data/*.cdf*')[0])
    epoch = data_dict[f'{FitDataToggles.epoch_key}'][0].copy()
    energy = data_dict[f'{FitDataToggles.energy_key}'][0].copy()
    pitch_angle = data_dict[f'{FitDataToggles.pitch_angle_key}'][0].copy()

    # Collect the relevant data
    diffEFlux = data_dict[f'{FitDataToggles.differential_energy_flux_key}'][0].copy()
    diffNFlux = np.array([np.divide(diffEFlux[tmeIdx].T, energy).T for tmeIdx in range(len(epoch))])
    counts = data_dict[f'{FitDataToggles.counts_key}'][0].copy()
    fit_data = [counts, diffNFlux, diffEFlux, ]

    # Standardize the format of the data [Epoch, pitch angle, Energy] where Pitch angle goes low to high, and energy high to low
    fit_obj = PrimaryBeamClasses(epoch.copy(),energy.copy(),pitch_angle.copy())
    for i in range(3):
        spectra_out, energy_out, pitch_out, info = fit_obj.standardize_flux_array(fit_data[i])
        fit_data[i] = spectra_out
        energy = energy_out
        pitch_angle = pitch_out


    ######################################
    # --- PREPARE THE DATA FOR FITTING ---
    ######################################

    # --- [0] REDUCE THE DATA TO THE FIT REGION ---
    fit_window = [
        dt.datetime(int(FileToggles.data_year), int(FileToggles.data_month), int(FileToggles.data_day), int(PrimaryBeamFitToggles.start_hour), int(PrimaryBeamFitToggles.start_minute)),
        dt.datetime(int(FileToggles.data_year), int(FileToggles.data_month), int(FileToggles.data_day), int(PrimaryBeamFitToggles.end_hour), int(PrimaryBeamFitToggles.end_minute)),
    ]
    low_idx, high_idx = np.abs(epoch - fit_window[0]).argmin(), np.abs(epoch - fit_window[1]).argmin()
    fit_data = [thing[low_idx:high_idx + 1] for thing in fit_data]
    fit_obj.epoch = epoch[low_idx:high_idx + 1]

    # --- Clean up and average the data based on user inputs ---
    for i in range(3):

        # --- [1] Remove bad values in the data ---
        fit_data[i][fit_data[i] < 0] = 0 # mask counts, diffNFlux, diffEFlux data <0

        # --- [2] Average the data over the desired pitch angle range ---
        fit_data[i] = fit_obj.pitch_average(fit_data[i],groups=PrimaryBeamFitToggles.pitch_angles_to_fit)

        # --- [3] Average the data over the desired time range with the cadence desired ---
        spectra_tme_avg,epoch_avg = fit_obj.time_average(fit_data[i],n_avg=PrimaryBeamFitToggles.N_time_avg)
        fit_data[i] = spectra_tme_avg
        epoch = epoch_avg

        # --- [4] Determine the error in each counts measurement ---
        fit_data[i][fit_data[i]<1] = 0 # clamp all values less than 1 to zero.

    # Note: the error in the averaged counts is: deltaN = (1/Num_of_ptchs_avged) * sqrt(N_ptch0 + N_ptch1 + ...) for a given time/energy
    counts_error = np.sqrt(fit_data[0])

    ######################################
    # --- PREPARE THE OUTPUT DATA DICT ---
    ######################################

    data_dict_output = {
        # 'Te': [np.zeros(len(fit_data[0])), {'DEPEND_0': f'{FitDataToggles.epoch_key}', 'DEPEND_1': f'{FitDataToggles.pitch_angle_key}', 'UNITS': 'eV', 'LABLAXIS': 'Te'}],
        # 'n': [np.zeros(len(fit_data[0])), {'DEPEND_0': f'{FitDataToggles.epoch_key}', 'DEPEND_1': f'{FitDataToggles.pitch_angle_key}', 'UNITS': 'cm!A-3!N', 'LABLAXIS': 'ne'}],
        # 'phi': [np.zeros(len(fit_data[0])), {'DEPEND_0': f'{FitDataToggles.epoch_key}', 'DEPEND_1': f'{FitDataToggles.pitch_angle_key}', 'UNITS': 'eV', 'LABLAXIS': '&phi;!B0!N'}],
        # 'kappa': [np.zeros(len(fit_data[0])), {'DEPEND_0': f'{FitDataToggles.epoch_key}', 'DEPEND_1': f'{FitDataToggles.pitch_angle_key}', 'UNITS': None, 'LABLAXIS': '&kappa;'}],
        # 'chi2_goodness_of_fit': [np.zeros(len(fit_data[0])), {'DEPEND_0': f'{FitDataToggles.epoch_key}', 'UNITS': None, 'LABLAXIS': '&chi;!A^2!N'}],
        # 'find_fit_data_engy_idx': [np.zeros(len(diffNFlux_tmeAvg)), {}],
        # 'fitted_N_points': [np.zeros(len(fit_data[0])), {'DEPEND_0': None, 'UNITS': None, 'LABLAXIS': 'Number of Fitted Points'}],
        'fitted_diffNFlux':[fit_data[1],{'DEPEND_0': 'epoch','DEPEND_1': 'pitch_angle','DEPEND_2': 'energy','UNITS':'1/cm^2-sr-s-eV'}],
        'fitted_counts_error': [counts_error, {'DEPEND_0': 'epoch','DEPEND_1': 'pitch_angle','DEPEND_2': 'energy'}],
        'fitted_counts': [fit_data[0], {'DEPEND_0': 'epoch','DEPEND_1': 'pitch_angle','DEPEND_2': 'energy', 'UNITS': 'counts'}],
        f'pitch_angle':[np.array([np.nanmean(arr) for arr in PrimaryBeamFitToggles.pitch_angles_to_fit]),{'UNITS': 'degrees'}],
        f'energy':[fit_obj.energy,{'UNITS':'eV'}],
        f'epoch': [fit_obj.epoch, data_dict[f'{FitDataToggles.epoch_key}'][1].copy()],
    }

    ########################################
    # --- PERFORM THE INVERTED-V FITTING ---
    ########################################
    #
    # for tmeIdx in tqdm(range(len(diffNFlux_tmeAvg))):
    #
    #     try:
    #
    #         diffNFlux_slice = diffNFlux_tmeAvg[tmeIdx]
    #
    #         # [1] find the acceleration potential from the peak in diffNFlux above an energy threshold
    #         engy_thresh_idx = np.abs(energy-PrimaryBeamFitToggles.energy_thesh).argmin()
    #         xData = energy[:engy_thresh_idx+1] if PrimaryBeamFitToggles.energy_bin_direction == 0 else energy[engy_thresh_idx:]
    #         yData = diffNFlux_slice[:engy_thresh_idx+1] if PrimaryBeamFitToggles.energy_bin_direction == 0 else diffNFlux_slice[engy_thresh_idx:]
    #
    #         phi0_guess_idx = yData.argmax()
    #         phi0_guess = xData[yData.argmax()]
    #
    #         # [2] collect the xData/yData for the fit above the energy threshold
    #         # (a) Collect the data
    #         xData_fit = xData[:phi0_guess_idx+1]
    #         yData_fit = yData[:phi0_guess_idx+1]
    #
    #         # (b) Form the guesses/boundaries
    #         fit_func, kwargs_dict = fit_obj.form_fit_params(phi0_guess)
    #
    #         # [3] Marquart-Levenburg Fitting
    #         # (a) Fit the data
    #         params, cov = curve_fit(fit_func, xData_fit, yData_fit, **kwargs_dict)
    #
    #         # (b) Calculate Chi2
    #         std_devs = fit_obj.calc_jN_error(counts_val=counts_tmeAvg[tmeIdx][:phi0_guess_idx+1], energy_value=xData_fit, pitch_angles=pitch_angle)
    #         chi2 = (1 / (len(params) - 1)) * sum([(fit_func(xData_fit[i], *params) - yData_fit[i]) ** 2 / (std_devs[i] ** 2) for i in range(len(xData_fit))])
    #
    #         # [4] Refine the fit using the kaeppler method
    #
    #         # [5] Store the results
    #         # --- update the data_dict_output ---
    #         data_dict_output['n'][0][tmeIdx] = params[0]
    #         data_dict_output['Te'][0][tmeIdx] = params[1]
    #         data_dict_output['phi'][0][tmeIdx] = params[2]
    #         data_dict_output['kappa'][0][tmeIdx] = np.nan if PrimaryBeamToggles.fit_dist == 'maxwellian' else params[3]
    #         data_dict_output['chi2'][0][tmeIdx] = chi2
    #         data_dict_output['N_fitted_points'][0][tmeIdx] = len(yData_fit)
    #         data_dict_output['find_fit_data_engy_idx'][0][tmeIdx] = phi0_guess_idx
    #
    #     except Exception as e:
    #         print(e)
    #         # --- update the data_dict_output ---
    #         data_dict_output['Te'][0][tmeIdx] = np.nan
    #         data_dict_output['n'][0][tmeIdx]= np.nan
    #         data_dict_output['phi'][0][tmeIdx] = np.nan
    #         data_dict_output['kappa'][0][tmeIdx] = np.nan
    #         data_dict_output['chi2'][0][tmeIdx] = np.nan
    #         data_dict_output['N_fitted_points'][0][tmeIdx] = np.nan

    # --- --- --- --- --- ---
    # --- OUTPUT THE DATA ---
    # --- --- --- --- --- ---
    outputPath = rf'{FileToggles.RUN_PATH}/primary_beam_fits/primary_beam_fit_T{PrimaryBeamFitToggles.start_hour}{PrimaryBeamFitToggles.start_minute}_to_T{PrimaryBeamFitToggles.end_hour}{PrimaryBeamFitToggles.end_minute}.cdf'
    stl.outputDataDict(outputPath=outputPath,data_dict=data_dict_output)