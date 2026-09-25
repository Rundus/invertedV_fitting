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
        dt.datetime(int(FileToggles.data_year), int(FileToggles.data_month), int(FileToggles.data_day), int(PrimaryBeamFitToggles.start_hour), int(PrimaryBeamFitToggles.start_minute),int(PrimaryBeamFitToggles.start_second)),
        dt.datetime(int(FileToggles.data_year), int(FileToggles.data_month), int(FileToggles.data_day), int(PrimaryBeamFitToggles.end_hour), int(PrimaryBeamFitToggles.end_minute),int(PrimaryBeamFitToggles.end_second)),
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
        spectra_tme_avg, epoch_avg = fit_obj.time_average(fit_data[i],n_avg=PrimaryBeamFitToggles.N_time_avg)
        fit_data[i] = spectra_tme_avg
        fit_obj.epoch = epoch_avg

        # --- [4] Determine the error in each counts measurement ---
        fit_data[i][fit_data[i]<1] = 0 # clamp all values less than 1 to zero.

    # Note: the error in the averaged counts is: deltaN = (1/Num_of_ptchs_avged) * sqrt(N_ptch0 + N_ptch1 + ...) for a given time/energy
    counts_error = np.sqrt(fit_data[0])

    # Average the desired geometric factors together
    ptch_idxs = np.array([[int(np.abs(pitch_angle - val).argmin()) for val in arr] for arr in PrimaryBeamFitToggles.pitch_angles_to_fit])
    geoFactor_avg = np.nanmean([FitDataToggles.geoFactor[idxs] for idxs in ptch_idxs],axis=1)
    fit_obj.geometric_factors = geoFactor_avg


    ######################################
    # --- PREPARE THE OUTPUT DATA DICT ---
    ######################################
    output_shape = np.zeros_like(fit_data[0][:,0,:]) # Shape (N_epoch, M_fitted_pitch)

    data_dict_output = {
        'Te': [output_shape.copy(), {'DEPEND_0': 'epoch', 'DEPEND_1': 'pitch_angle', 'UNITS': 'eV', 'LABLAXIS': 'Te'}],
        'n': [output_shape.copy(), {'DEPEND_0': 'epoch', 'DEPEND_1': 'pitch_angle', 'UNITS': 'cm!A-3!N', 'LABLAXIS': 'ne'}],
        'phi': [output_shape.copy(), {'DEPEND_0': 'epoch', 'DEPEND_1': 'pitch_angle', 'UNITS': 'eV', 'LABLAXIS': '&Phi;!B0!N'}],
        'chi2_goodness_of_fit': [output_shape.copy(), {'DEPEND_0': 'epoch', 'DEPEND_1': 'pitch_angle', 'UNITS': None, 'LABLAXIS': '&chi;!A2!N'}],
        'fitted_N_points': [output_shape.copy(), {'DEPEND_0': 'epoch', 'DEPEND_1': 'pitch_angle', 'LABLAXIS': 'Number of Fitted Points'}],
        'kappa': [output_shape.copy(), {'DEPEND_0': 'epoch', 'DEPEND_1': 'pitch_angle', 'UNITS': None, 'LABLAXIS': '&kappa;'}],
        'fitted_diffNFlux':[fit_data[1],{'DEPEND_0': 'epoch','DEPEND_2': 'pitch_angle','DEPEND_1': 'energy','UNITS':'1/cm^2-sr-s-eV'}],
        'fitted_diffNFlux_error': [output_shape, {'DEPEND_0': 'epoch', 'DEPEND_2': 'pitch_angle', 'DEPEND_1': 'energy', 'UNITS': '1/cm^2-sr-s-eV'}],
        'fitted_counts_error': [np.zeros_like(fit_data[0]), {'DEPEND_0': 'epoch','DEPEND_2': 'pitch_angle','DEPEND_1': 'energy'}],
        'fitted_counts': [fit_data[0], {'DEPEND_0': 'epoch','DEPEND_2': 'pitch_angle','DEPEND_1': 'energy', 'UNITS': 'counts'}],
        f'pitch_angle':[np.array([np.nanmean(arr) for arr in PrimaryBeamFitToggles.pitch_angles_to_fit]),{'UNITS': 'degrees'}],
        f'energy':[fit_obj.energy,{'UNITS':'eV'}],
        f'epoch': [fit_obj.epoch, data_dict[f'{FitDataToggles.epoch_key}'][1].copy()],
    }

    ########################################
    # --- PERFORM THE INVERTED-V FITTING ---
    ########################################

    for tmeIdx in tqdm(range(len(data_dict_output['epoch'][0]))):
        for ptchIdx in range(len(data_dict_output['pitch_angle'][0])):

            try:

                spectra_fit = data_dict_output['fitted_diffNFlux'][0][tmeIdx,:,ptchIdx].copy()
                counts_fit = data_dict_output['fitted_counts'][0][tmeIdx,:,ptchIdx].copy()
                energy_fit = data_dict_output['energy'][0].copy()


                # --- [1] find the acceleration potential from the peak in diffNFlux above an energy threshold ---

                # Reduce fitted data to that above energy threshold
                engy_thresh_idx = np.abs(energy_fit-PrimaryBeamFitToggles.energy_thesh).argmin()
                xData = energy_fit[:engy_thresh_idx+1]
                yData = spectra_fit[:engy_thresh_idx+1]
                countsData = counts_fit[:engy_thresh_idx+1]

                # Remove any points where counts == 0
                bad_idxs = np.where(countsData == 0)
                xData = np.delete(xData,bad_idxs)
                yData = np.delete(yData, bad_idxs)
                countsData = np.delete(countsData,bad_idxs)

                # check if number of fitted points >= N_fit_min:
                if len(yData) < PrimaryBeamFitToggles.N_fit_min:
                    raise Exception('Too Few Fit Points')

                # Iterate: Find which choice of Phi0 gives the lowest Chi2 value:
                chi2_min = np.inf
                params_best = []
                N_fit_best = 0

                for i in range(len(yData) - PrimaryBeamFitToggles.N_fit_min):

                    iter_idx = i +PrimaryBeamFitToggles.N_fit_min
                    xData_iter = xData[:iter_idx]
                    yData_iter = yData[:iter_idx]
                    countsData_iter = countsData[:iter_idx]

                    # estimate the initial guess for characteristic inverted-V energy via the peak in diffNFlux energy location
                    phi0_guess_idx = yData_iter.argmax()
                    phi0_guess = xData[yData_iter.argmax()]

                    # Form the initial guesses/fitting function
                    fit_func, kwargs_dict = fit_obj.form_fit_params(phi0_guess)

                    # --- [2] Perform the Marquart-Levenburg Fit ---
                    # Fit the data
                    params, cov = curve_fit(fit_func, xData_iter, yData_iter, **kwargs_dict)

                    # Calculate error in each datapoint
                    std_devs = fit_obj.calc_jN_error(fit_counts=countsData_iter,
                                                     fit_energies=xData_iter,
                                                     pitch_idx = ptchIdx)

                    # Calculate Chi2
                    chi2 = (1 / (len(params) - 1)) * sum([(fit_func(xData_iter[i], *params) - yData_iter[i]) ** 2 / (std_devs[i] ** 2) for i in range(len(xData_iter))])

                    if np.all([chi2<chi2_min, chi2>1E-3]):
                        chi2_min = chi2
                        params_best = params
                        N_fit_best = len(xData_iter)

                # [4] Refine the fit using the kaeppler method
                # TODO

                # [5] Store the results
                n_fit = params_best[0]
                Te_fit = params_best[1]
                phi_fit =  params_best[2]
                kappa_fit = np.nan if PrimaryBeamFitToggles.fit_dist == 'maxwellian' else params_best[3]
                ch2_fit = chi2_min
                N_fit =  N_fit_best

            except Exception as e:
                print(e)
                n_fit = np.nan
                Te_fit = np.nan
                phi_fit = np.nan
                kappa_fit = np.nan
                ch2_fit = np.nan
                N_fit = np.nan

            # store the output
            data_dict_output['n'][0][tmeIdx, ptchIdx] = n_fit
            data_dict_output['Te'][0][tmeIdx, ptchIdx] = Te_fit
            data_dict_output['phi'][0][tmeIdx, ptchIdx] = phi_fit
            data_dict_output['kappa'][0][tmeIdx, ptchIdx] = kappa_fit
            data_dict_output['chi2_goodness_of_fit'][0][tmeIdx, ptchIdx] = ch2_fit
            data_dict_output['fitted_N_points'][0][tmeIdx, ptchIdx] = N_fit

    # --- --- --- --- --- ---
    # --- OUTPUT THE DATA ---
    # --- --- --- --- --- ---
    outputPath = rf'{FileToggles.RUN_PATH}/primary_beam_fits/primary_beam_fit_T{PrimaryBeamFitToggles.start_hour}{PrimaryBeamFitToggles.start_minute}{PrimaryBeamFitToggles.start_second}_to_{PrimaryBeamFitToggles.end_hour}{PrimaryBeamFitToggles.end_minute}{PrimaryBeamFitToggles.end_second}.cdf'
    stl.outputDataDict(outputPath=outputPath,data_dict=data_dict_output)