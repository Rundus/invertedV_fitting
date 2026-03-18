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


    # --- prepare the output ---
    data_dict_output = {
                'Te': [[], {'DEPEND_0': f'{UserToggles.Epoch_key}', 'DEPEND_1':f'{UserToggles.pitch_angle_key}',  'UNITS': 'eV', 'LABLAXIS': 'Te'}],
                 'n': [[], {'DEPEND_0': f'{UserToggles.Epoch_key}', 'DEPEND_1':f'{UserToggles.pitch_angle_key}','UNITS': 'cm!A-3!N', 'LABLAXIS': 'ne'}],
                 'phi0': [[], {'DEPEND_0': f'{UserToggles.Epoch_key}', 'DEPEND_1':f'{UserToggles.pitch_angle_key}','UNITS': 'eV', 'LABLAXIS': '%phi;!B0!N'}],
                 'kappa': [[], {'DEPEND_0': f'{UserToggles.Epoch_key}', 'DEPEND_1':f'{UserToggles.pitch_angle_key}','UNITS': None, 'LABLAXIS': '&kappa;'}],
                 'chi2': [[], {'DEPEND_0': f'{UserToggles.Epoch_key}', 'UNITS': None, 'LABLAXIS': '&chi;!A^2!N'}],
                 f'{UserToggles.Epoch_key}': [[], {'DEPEND_0': None, 'UNITS': 'ns', 'LABLAXIS': '&phi;0'}],
                 'N_fitted_points': [[], {'DEPEND_0': None, 'UNITS': None, 'LABLAXIS': 'Number of Fitted Points'}],
                 }

    ###############################
    # --- LOADING THE FLUX DATA ---
    ###############################
    data_dict_flux = stl.loadDictFromFile(glob(UserToggles.path_to_eESA_data)[0])
    epoch = data_dict_flux[f'{UserToggles.Epoch_key}'][0]
    energy = data_dict_flux[f'{UserToggles.energy_key}'][0]
    pitch_angle = data_dict_flux[f'{UserToggles.pitch_angle_key}'][0]
    diffEFlux = data_dict_flux[f'{UserToggles.differential_energy_flux_key}'][0]
    diffNFlux = np.divide(deepcopy(diffEFlux),energy,axis=UserToggles.dependency_structure.index(UserToggles.energy_key))

    ################################
    # --- PREPARE THE INPUT DATA ---
    ################################

    # [1] Average the data over the desired pitch angle range
    dependency_indices = []
    for i in range(len(UserToggles.dependency_structure)):
        if 'epoch' in UserToggles.dependency_structure[i]:
            dependency_indices.append([i for i in range(len(epoch))])
        elif 'energy' in UserToggles.dependency_structure[i]:
            dependency_indices.append([i for i in range(len(energy))])
        elif 'pitch' in UserToggles.dependency_structure[i]:
            dependency_indices.append([i for i in range(len(pitch_angle)) if pitch_angle[i] in UserToggles.pitch_angles_to_fit])

    diffNFlux_ptchAvg = np.nanmean(diffNFlux[*dependency_indices], axis=UserToggles.dependency_structure.index(UserToggles.pitch_angle_key))

    # [2] Average the data over the desired time range with the cadence desired
    diffNFlux_tmeAvg = diffNFlux_ptchAvg


    ##################################
    # --- LOOP THROUGH DATA TO FIT ---
    ##################################
    noiseData = helperFuncs().generateNoiseLevel(data_dict_diffFlux['Energy'][0], countNoiseLevel=primaryBeamToggles.countNoiseLevel)
    for loopIdx, pitchAngle in enumerate(primaryBeamToggles.wPitchsToFit):
        ####################################
        # --- FIT AVERAGED TIME SECTIONS ---
        ####################################
        # for each slice in time, loop over the data and identify the peak differentialNumberFlux (This corresponds to the
        # peak energy of the inverted-V since the location of the maximum number flux tells you what energy the low-energy BULk got accelerated to)
        # Note: The peak in the number flux is very likely the maximum value AFTER 100 eV, just find this point
        pitchIdx = np.abs(data_dict_diffFlux['Pitch_Angle'][0] - pitchAngle).argmin()

        for tmeIdx in tqdm(fit_idxs):
            try: # try fit the data
                # --- Determine the accelerated potential from the peak in diffNflux based on a threshold limit ---
                engythresh_Idx = np.abs(data_dict_diffFlux['Energy'][0] - primaryBeamToggles.engy_Thresh).argmin()  # only consider data above a certain index
                peakDiffNIdx = np.nanargmax(fitData[tmeIdx][pitchIdx][:engythresh_Idx + 1])  # only consider data above the Energy_Threshold, to avoid secondaries/backscatter

                # Determine which datapoints are good to fit
                # Description: After finding the peak in the jN spectrum, fit the data starting 1 BEFORE this peak value
                peakIdx = peakDiffNIdx  # Kaeppler added +1 here to make his ChiSquare's get better. Our are MUCH better without doing that
                dataIdxs = np.array([1 if fitData[tmeIdx][pitchIdx][j] > noiseData[j] and j <= peakIdx else 0 for j in range(len(data_dict_diffFlux['Energy'][0]))])

                ###################################
                # ---  get the Beam data subset ---
                ###################################
                fitTheseIndicies = np.where(dataIdxs == 1)[0]
                xData_fit, yData_fit, yData_fit_stdDev = np.array(data_dict_diffFlux['Energy'][0][fitTheseIndicies]), np.array(fitData[tmeIdx][pitchIdx][fitTheseIndicies]), np.array(fitData_stdDev[tmeIdx][pitchIdx][fitTheseIndicies])
                V0_guess = xData_fit[-1]

                # Only include data with non-zero points
                nonZeroIndicies = np.where(yData_fit!=0)[0]
                xData_fit, yData_fit, yData_fit_stdDev = xData_fit[nonZeroIndicies],  yData_fit[nonZeroIndicies], yData_fit_stdDev[nonZeroIndicies]

                ### Perform the fit ###
                params, ChiSquare = fitPrimaryBeam(xData_fit, yData_fit, yData_fit_stdDev, V0_guess, primaryBeamToggles, useNoGuess=primaryBeamToggles.useNoGuess)

                ### Kaeppler fit refinement ###
                if primaryBeamToggles.useFitRefinement:

                    ### Use the First fit to motivate the n0 guess ###
                    firstFitParams = deepcopy(params)  # construct the guess

                    # define the fit function with specific charge/mass
                    if primaryBeamToggles.wDistributionToFit == 'Maxwellian':
                        firstFitParams[3] = 10

                    ### Use the First fit to motivate the n0 guess ###
                    n0guess = n0GuessKaeppler2014(fitData[tmeIdx],
                                                  firstFitParams,
                                                  peakDiffNIdx,
                                                  beta=primaryBeamToggles.beta_guess,
                                                  energyRange=data_dict_diffFlux['Energy'][0])

                    ### Perform the informed fit again ###
                    newBounds = [[n0guess*(1-primaryBeamToggles.n0guess_deviation), n0guess*(1+primaryBeamToggles.n0guess_deviation)],
                                 primaryBeamToggles.Te_bounds,
                                 [(1 - primaryBeamToggles.V0_deviation) * V0_guess, (1 + primaryBeamToggles.V0_deviation) * V0_guess], #V0
                                 primaryBeamToggles.kappa_bounds # kappa
                                ]
                    params, ChiSquare = fitPrimaryBeam(xData_fit, yData_fit, yData_fit_stdDev, V0_guess, primaryBeamToggles,
                                                       specifyGuess = [n0guess, params[1], params[2], params[3]],
                                                       specifyBoundVals = newBounds,
                                                       useNoGuess = primaryBeamToggles.useNoGuess
                                                       )

                # --- update the data_dict_output ---
                data_dict_output['Te'][0][tmeIdx][pitchIdx] = params[1]
                data_dict_output['n'][0][tmeIdx][pitchIdx] = params[0]
                data_dict_output['V0'][0][tmeIdx][pitchIdx] = params[2]
                data_dict_output['kappa'][0][tmeIdx][pitchIdx] = 0 if primaryBeamToggles.wDistributionToFit == 'Maxwellian' else params[3]
                data_dict_output['ChiSquare'][0][tmeIdx][pitchIdx]= ChiSquare
                data_dict_output['dataIdxs'][0][tmeIdx][pitchIdx] = dataIdxs
                data_dict_output['numFittedPoints'][0][tmeIdx][pitchIdx] = len(yData_fit)
            except: # output all nans
                # --- update the data_dict_output ---
                data_dict_output['Te'][0][tmeIdx][pitchIdx] = np.nan
                data_dict_output['n'][0][tmeIdx][pitchIdx] = np.nan
                data_dict_output['V0'][0][tmeIdx][pitchIdx] = np.nan
                data_dict_output['kappa'][0][tmeIdx][pitchIdx] = np.nan
                data_dict_output['ChiSquare'][0][tmeIdx][pitchIdx] = np.nan
                data_dict_output['dataIdxs'][0][tmeIdx][pitchIdx] = np.nan
                data_dict_output['numFittedPoints'][0][tmeIdx][pitchIdx] = np.nan

    # --- --- --- --- --- ---
    # --- OUTPUT THE DATA ---
    # --- --- --- --- --- ---
    # Construct the Data Dict
    exampleVar = {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808,
                  'FORMAT': 'I5', 'UNITS': 'm', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                  'SCALETYP': 'linear', 'LABLAXIS': None}

    # update the data dict attrs
    for key, val in data_dict_output.items():

        # convert the data to numpy arrays
        data_dict_output[key][0] = np.array(data_dict_output[key][0])

        # update the attributes
        newAttrs = deepcopy(exampleVar)

        for subKey, subVal in data_dict_output[key][1].items():
            newAttrs[subKey] = subVal

        data_dict_output[key][1] = newAttrs

    outputPath = rf'{outputFolder}\primaryBeam_fitting_parameters.cdf'
    stl.outputCDFdata(outputPath=outputPath,data_dict=data_dict_output)