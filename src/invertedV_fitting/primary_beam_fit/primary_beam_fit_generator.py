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
    from scipy.optimize import curve_fit
    import datetime as dt

    ###############################
    # --- LOADING THE FLUX DATA ---
    ###############################
    data_dict = stl.loadDictFromFile(glob(FitDataToggles.path_to_eESA_flux_data)[0])
    epoch = data_dict[f'{FitDataToggles.epoch_key}'][0]
    energy = data_dict[f'{FitDataToggles.energy_key}'][0]
    pitch_angle = data_dict[f'{FitDataToggles.pitch_angle_key}'][0]

    # Clean up the diffEFlux data
    diffEFlux = data_dict[f'{FitDataToggles.differential_energy_flux_key}'][0]

    # Clean up the Counts Data
    counts = data_dict[f'{FitDataToggles.counts_key}'][0]

    # Calculate diffNFlux
    diffNFlux = np.array([np.divide(diffEFlux[tmeIdx].T, energy).T for tmeIdx in range(len(epoch))])

    ######################################
    # --- PREPARE THE DATA FOR FITTING ---
    ######################################

    # --- [0] REDUCE THE DATA TO THE FIT REGION ---
    low_idx, high_idx = np.abs(epoch - PrimaryBeamFitToggles.datetime_low).argmin(), np.abs(epoch - PrimaryBeamFitToggles.datetime_high).argmin()
    epoch = epoch[low_idx:high_idx + 1]
    fit_data = [counts, diffNFlux, diffEFlux]
    fit_data = [thing[low_idx:high_idx + 1] for thing in fit_data]

    # --- [1] Remove bad values in the data ---
    fit_data[1][fit_data[1] < 0] = 0 # mask diffNFlux <0
    fit_data[0][fit_data[0] < 0] = 0 # mask counts <0


    # --- [2] Average the data over the desired pitch angle range ---
    dependency_indices = [[],[],[]]
    for i in range(len(FitDataToggles.dependency_structure)):
        wIdx = FitDataToggles.dependency_structure.index(FitDataToggles.dependency_structure[i])
        if 'epoch' in FitDataToggles.dependency_structure[i].lower():
            dependency_indices[wIdx] =[i for i in range(len(epoch))]
        elif 'energy' in FitDataToggles.dependency_structure[i].lower():
            dependency_indices[wIdx] =[i for i in range(len(energy))]
        elif 'pitch' in FitDataToggles.dependency_structure[i].lower():
            dependency_indices[wIdx] =[i for i in range(len(pitch_angle)) if pitch_angle[i] in PrimaryBeamFitToggles.pitch_angles_to_fit]

    diffNFlux_ptch_avg = np.nanmean(fit_data[1][*np.ix_(*dependency_indices)], axis=FitDataToggles.dependency_structure.index(FitDataToggles.pitch_angle_key))
    counts_ptch_avg = np.round(np.nanmean(fit_data[0][*np.ix_(*dependency_indices)], axis=FitDataToggles.dependency_structure.index(FitDataToggles.pitch_angle_key)))

    # --- [3] Average the data over the desired time range with the cadence desired ---
    fitting_window = [
                        dt.datetime(int(FileToggles.data_year), int(FileToggles.data_month), int(FileToggles.data_day), PrimaryBeamFitToggles.start_hour, PrimaryBeamFitToggles.start_minute),
        dt.datetime(int(FileToggles.data_year), int(FileToggles.data_month), int(FileToggles.data_day),
                    PrimaryBeamFitToggles.start_hour, PrimaryBeamFitToggles.start_minute),
                      ]
    low_idx, high_idx = np.abs(epoch - fitting_window[0]).argmin(),np.abs(epoch-fitting_window[0]).argmin()
    diffNFlux_tmeAvg = diffNFlux_ptch_avg[low_idx:high_idx+1]
    counts_tmeAvg = counts_ptch_avg[low_idx:high_idx+1]
    counts_tmeAvg[counts_tmeAvg<0] = 0

    # --- [4] Determine the error in each counts measurement ---
    # Note: the error in the averaged counts is: deltaN = (1/Num_of_ptchs_avged) * sqrt(N_ptch0 + N_ptch1 + ...) for a given time/energy
    # counts_stdDev = (1/len(PrimaryBeamToggles.pitch_angles_to_fit))*np.sqrt(np.nansum(counts[*np.ix_(*dependency_indices)], axis=UserToggles.dependency_structure.index(UserToggles.pitch_angle_key))).round()[low_idx:high_idx+1]
    counts_stdDev = np.sqrt(counts_tmeAvg)

    ######################################
    # --- PREPARE THE OUTPUT DATA DICT ---
    ######################################

    data_dict_output = {
        f'{FitDataToggles.epoch_key}': [epoch[low_idx:high_idx+1],data_dict[f'{FitDataToggles.epoch_key}'][1].copy()],
        'Te': [np.zeros(len(diffNFlux_tmeAvg)), {'DEPEND_0': f'{FitDataToggles.epoch_key}', 'DEPEND_1': f'{FitDataToggles.pitch_angle_key}', 'UNITS': 'eV', 'LABLAXIS': 'Te'}],
        'n': [np.zeros(len(diffNFlux_tmeAvg)), {'DEPEND_0': f'{FitDataToggles.epoch_key}', 'DEPEND_1': f'{FitDataToggles.pitch_angle_key}', 'UNITS': 'cm!A-3!N', 'LABLAXIS': 'ne'}],
        'phi': [np.zeros(len(diffNFlux_tmeAvg)), {'DEPEND_0': f'{FitDataToggles.epoch_key}', 'DEPEND_1': f'{FitDataToggles.pitch_angle_key}', 'UNITS': 'eV', 'LABLAXIS': '&phi;!B0!N'}],
        'kappa': [np.zeros(len(diffNFlux_tmeAvg)), {'DEPEND_0': f'{FitDataToggles.epoch_key}', 'DEPEND_1': f'{FitDataToggles.pitch_angle_key}', 'UNITS': None, 'LABLAXIS': '&kappa;'}],
        'chi2_goodness_of_fit': [np.zeros(len(diffNFlux_tmeAvg)), {'DEPEND_0': f'{FitDataToggles.epoch_key}', 'UNITS': None, 'LABLAXIS': '&chi;!A^2!N'}],
        'fitted_N_points': [np.zeros(len(diffNFlux_tmeAvg)), {'DEPEND_0': None, 'UNITS': None, 'LABLAXIS': 'Number of Fitted Points'}],
        'fitted_diffNFlux':[diffNFlux_tmeAvg,{'DEPEND_0':'Epoch','DEPEND_1':'Energy','UNITS':'1/cm^2-str-s-eV'}],
        'fitted_counts_std': [counts_stdDev, {'DEPEND_0': 'Epoch', 'DEPEND_1': 'Energy','UNITS':'counts'}],
        'fitted_counts': [counts_tmeAvg, {'DEPEND_0': 'Epoch', 'DEPEND_1': 'Energy', 'UNITS': 'counts'}],
        f'{FitDataToggles.pitch_angle_key}':data_dict[f'{FitDataToggles.pitch_angle_key}'],
        f'{FitDataToggles.energy_key}':data_dict[f'{FitDataToggles.energy_key}'],
        'find_fit_data_engy_idx' : [np.zeros(len(diffNFlux_tmeAvg)),{}],
    }

    ########################################
    # --- PERFORM THE INVERTED-V FITTING ---
    ########################################

    for tmeIdx in tqdm(range(len(diffNFlux_tmeAvg))):

        try:

            diffNFlux_slice = diffNFlux_tmeAvg[tmeIdx]

            # [1] find the acceleration potential from the peak in diffNFlux above an energy threshold
            engy_thresh_idx = np.abs(energy-PrimaryBeamFitToggles.energy_thesh).argmin()
            xData = energy[:engy_thresh_idx+1] if PrimaryBeamFitToggles.energy_bin_direction == 0 else energy[engy_thresh_idx:]
            yData = diffNFlux_slice[:engy_thresh_idx+1] if PrimaryBeamFitToggles.energy_bin_direction == 0 else diffNFlux_slice[engy_thresh_idx:]

            phi0_guess_idx = yData.argmax()
            phi0_guess = xData[yData.argmax()]

            # [2] collect the xData/yData for the fit above the energy threshold
            # (a) Collect the data
            xData_fit = xData[:phi0_guess_idx+1]
            yData_fit = yData[:phi0_guess_idx+1]

            # (b) Form the guesses/boundaries
            fit_func, kwargs_dict = PrimaryBeamClasses().form_fit_params(phi0_guess)

            # [3] Marquart-Levenburg Fitting
            # (a) Fit the data
            params, cov = curve_fit(fit_func, xData_fit, yData_fit, **kwargs_dict)

            # (b) Calculate Chi2
            std_devs = PrimaryBeamClasses().calc_jN_error(counts_val=counts_tmeAvg[tmeIdx][:phi0_guess_idx+1], energy_value=xData_fit, pitch_angles=pitch_angle)
            chi2 = (1 / (len(params) - 1)) * sum([(fit_func(xData_fit[i], *params) - yData_fit[i]) ** 2 / (std_devs[i] ** 2) for i in range(len(xData_fit))])

            # [4] Refine the fit using the kaeppler method

            # [5] Store the results
            # --- update the data_dict_output ---
            data_dict_output['n'][0][tmeIdx] = params[0]
            data_dict_output['Te'][0][tmeIdx] = params[1]
            data_dict_output['phi'][0][tmeIdx] = params[2]
            data_dict_output['kappa'][0][tmeIdx] = np.nan if PrimaryBeamToggles.fit_dist == 'maxwellian' else params[3]
            data_dict_output['chi2'][0][tmeIdx] = chi2
            data_dict_output['N_fitted_points'][0][tmeIdx] = len(yData_fit)
            data_dict_output['find_fit_data_engy_idx'][0][tmeIdx] = phi0_guess_idx

        except Exception as e:
            print(e)
            # --- update the data_dict_output ---
            data_dict_output['Te'][0][tmeIdx] = np.nan
            data_dict_output['n'][0][tmeIdx]= np.nan
            data_dict_output['phi'][0][tmeIdx] = np.nan
            data_dict_output['kappa'][0][tmeIdx] = np.nan
            data_dict_output['chi2'][0][tmeIdx] = np.nan
            data_dict_output['N_fitted_points'][0][tmeIdx] = np.nan

    # --- --- --- --- --- ---
    # --- OUTPUT THE DATA ---
    # --- --- --- --- --- ---
    outputPath = rf'{FileToggles.RUN_PATH}/primary_beam_fits/primary_beam_fit_{}_to_{}.cdf'
    stl.outputDataDict(outputPath=outputPath,data_dict=data_dict_output)