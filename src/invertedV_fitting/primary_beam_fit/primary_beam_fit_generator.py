# --- primary_beam_fit_generator.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: using the method outline in Kaeppler's thesis, we can fit inverted-V distributions
# to get estimate the magnetospheric temperature, density and electrostatic potential that accelerated
# our particles

# TODO: Implement Time-averaging into fit routine

from timebudget import timebudget

@timebudget
def primary_beam_fit_generator():

    # --- imports ---
    import spaceToolsLib as stl
    import numpy as np
    from tqdm import tqdm
    from copy import deepcopy
    from glob import glob
    from src.invertedV_fitting.user_toggles.user_toggles import UserToggles
    from src.invertedV_fitting.primary_beam_fit.primary_beam_fit_toggles import PrimaryBeamToggles
    from src.invertedV_fitting.primary_beam_fit.primary_beam_fit_classes import PrimaryBeamClasses
    from scipy.optimize import curve_fit


    ###############################
    # --- LOADING THE FLUX DATA ---
    ###############################
    data_dict_flux = stl.loadDictFromFile(glob(UserToggles.path_to_eESA_flux_data)[0])
    data_dict_counts = stl.loadDictFromFile(glob(UserToggles.path_to_eESA_counts_data)[0])
    epoch = deepcopy(data_dict_flux[f'{UserToggles.Epoch_key}'][0])
    energy = deepcopy(data_dict_flux[f'{UserToggles.energy_key}'][0])
    pitch_angle = deepcopy(data_dict_flux[f'{UserToggles.pitch_angle_key}'][0])

    diffEFlux = deepcopy(data_dict_flux[f'{UserToggles.differential_energy_flux_key}'][0])
    diffEFlux[diffEFlux<0] = 0
    diffNFlux = np.array([np.divide(diffEFlux[tmeIdx].T,energy).T for tmeIdx in range(len(epoch)) ])

    counts = deepcopy(data_dict_counts[f'{UserToggles.counts_key}'][0])
    counts[counts<0] = 0

    ################################
    # --- PREPARE THE INPUT DATA ---
    ################################

    # [1] Average the data over the desired pitch angle range
    dependency_indices = [[],[],[]]
    for i in range(len(UserToggles.dependency_structure)):
        wIdx = UserToggles.dependency_structure.index(UserToggles.dependency_structure[i])
        if 'epoch' in UserToggles.dependency_structure[i].lower():
            dependency_indices[wIdx] =[i for i in range(len(epoch))]
        elif 'energy' in UserToggles.dependency_structure[i].lower():
            dependency_indices[wIdx] =[i for i in range(len(energy))]
        elif 'pitch' in UserToggles.dependency_structure[i].lower():
            dependency_indices[wIdx] =[i for i in range(len(pitch_angle)) if pitch_angle[i] in PrimaryBeamToggles.pitch_angles_to_fit]

    diffNFlux_ptchAvg = np.nanmean(diffNFlux[*np.ix_(*dependency_indices)], axis=UserToggles.dependency_structure.index(UserToggles.pitch_angle_key))
    counts_ptchAvg = np.round(np.nanmean(counts[*np.ix_(*dependency_indices)], axis=UserToggles.dependency_structure.index(UserToggles.pitch_angle_key)))

    # [2] Average the data over the desired time range with the cadence desired
    low_idx, high_idx = np.abs(epoch - UserToggles.datetime_low).argmin(),np.abs(epoch-UserToggles.datetime_high).argmin()
    diffNFlux_tmeAvg = deepcopy(diffNFlux_ptchAvg)[low_idx:high_idx+1]
    counts_tmeAvg = deepcopy(counts_ptchAvg)[low_idx:high_idx+1]
    counts_tmeAvg[counts_tmeAvg<0] = 0

    # [3] Determine the error in each counts measurement
    # Note: the error in the averaged counts is: deltaN = (1/Num_of_ptchs_avged) * sqrt(N_ptch0 + N_ptch1 + ...) for a given time/energy
    # counts_stdDev = (1/len(PrimaryBeamToggles.pitch_angles_to_fit))*np.sqrt(np.nansum(counts[*np.ix_(*dependency_indices)], axis=UserToggles.dependency_structure.index(UserToggles.pitch_angle_key))).round()[low_idx:high_idx+1]
    counts_stdDev = np.sqrt(counts_tmeAvg)


    #################################
    # --- PREPARE THE OUTPUT DATA ---
    #################################

    data_dict_output = {
        'Epoch': [epoch[low_idx:high_idx+1],deepcopy(data_dict_flux[f'{UserToggles.Epoch_key}'][1])],
        'Te': [np.zeros(len(diffNFlux_tmeAvg)), {'DEPEND_0': f'{UserToggles.Epoch_key}', 'DEPEND_1': f'{UserToggles.pitch_angle_key}', 'UNITS': 'eV', 'LABLAXIS': 'Te'}],
        'n': [np.zeros(len(diffNFlux_tmeAvg)), {'DEPEND_0': f'{UserToggles.Epoch_key}', 'DEPEND_1': f'{UserToggles.pitch_angle_key}', 'UNITS': 'cm!A-3!N', 'LABLAXIS': 'ne'}],
        'phi': [np.zeros(len(diffNFlux_tmeAvg)), {'DEPEND_0': f'{UserToggles.Epoch_key}', 'DEPEND_1': f'{UserToggles.pitch_angle_key}', 'UNITS': 'eV', 'LABLAXIS': '&phi;!B0!N'}],
        'kappa': [np.zeros(len(diffNFlux_tmeAvg)), {'DEPEND_0': f'{UserToggles.Epoch_key}', 'DEPEND_1': f'{UserToggles.pitch_angle_key}', 'UNITS': None, 'LABLAXIS': '&kappa;'}],
        'chi2': [np.zeros(len(diffNFlux_tmeAvg)), {'DEPEND_0': f'{UserToggles.Epoch_key}', 'UNITS': None, 'LABLAXIS': '&chi;!A^2!N'}],
        'N_fitted_points': [np.zeros(len(diffNFlux_tmeAvg)), {'DEPEND_0': None, 'UNITS': None, 'LABLAXIS': 'Number of Fitted Points'}],
        'diffNFlux':[diffNFlux_tmeAvg,{'DEPEND_0':'Epoch','DEPEND_1':'Energy','UNITS':'1/cm^2-str-s-eV'}],
        'counts_std': [counts_stdDev, {'DEPEND_0': 'Epoch', 'DEPEND_1': 'Energy','UNITS':'counts'}],
        'counts': [counts_tmeAvg, {'DEPEND_0': 'Epoch', 'DEPEND_1': 'Energy', 'UNITS': 'counts'}],
        'pitch_angle':deepcopy(data_dict_flux[f'{UserToggles.pitch_angle_key}']),
        'Energy':deepcopy(data_dict_flux[f'{UserToggles.energy_key}']),
        'find_fit_data_engy_idx' : [np.zeros(len(diffNFlux_tmeAvg)),{}]
    }

    ##################################
    # --- LOOP THROUGH DATA TO FIT ---
    ##################################

    for tmeIdx in tqdm(range(len(diffNFlux_tmeAvg))):

        try:

            diffNFlux_slice = diffNFlux_tmeAvg[tmeIdx]

            # [1] find the acceleration potential from the peak in diffNFlux above an energy threshold
            engy_thresh_idx = np.abs(energy-PrimaryBeamToggles.energy_thesh).argmin()
            xData = energy[:engy_thresh_idx+1] if UserToggles.energy_bin_direction == 0 else energy[engy_thresh_idx:]
            yData = diffNFlux_slice[:engy_thresh_idx+1] if UserToggles.energy_bin_direction == 0 else diffNFlux_slice[engy_thresh_idx:]

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
    outputPath = rf'{UserToggles.run_folder_path}/primary_beam_fit.cdf'
    stl.outputDataDict(outputPath=outputPath,data_dict=data_dict_output)