
##########################
# --- PRIMARY BEAM FIT ---
##########################
from datetime import datetime
class primaryBeamToggles:

    # --- Which Regions to fit ---
    wFit_times = [0]
    invertedV_times = [
                        # ACESII HIGH FLYER
                        [datetime(2022, 11, 20, 17, 25,  1, 000000), datetime(2022, 11, 20, 17, 25, 3, 000000)], # 0 Dispersive Region
                        [datetime(2022, 11, 20, 17, 24, 12, 162000), datetime(2022, 11, 20, 17, 24, 18, 812000)], # 1 Very First ,Inverted-V, the high energy one
                        [datetime(2022, 11, 20, 17, 24, 45, 862000), datetime(2022, 11, 20, 17, 24, 49, 312000)], # 2 small inverted-V, after the High energy One
                        [datetime(2022, 11, 20, 17, 25, 23, 762000), datetime(2022, 11, 20, 17, 26, 8, 212000)],  # 3 Primary inverted-V
                        [datetime(2022, 11, 20, 17, 26, 11, 412000), datetime(2022, 11, 20, 17, 26, 19, 912000)],  # 4 Inverted-V right after the Primary-V, has STEBs on either sides of it
                        [datetime(2022, 11, 20, 17, 26, 35, 112000), datetime(2022, 11, 20, 17, 26, 40, 712000)], # 5 Inverted-V two after the Primary-V
                        [datetime(2022, 11, 20, 17, 28, 17, 112000), datetime(2022, 11, 20, 17, 28, 34, 612000)] # 6 Faint inverted-V on the most northside of the flight

                        # ACESII LOW FLYER
                        [datetime(2009, 1, 29, 9, 54, 4, 0), datetime(2009, 1, 29, 9, 54, 29, 000)] # 7 Very First ,Inverted-V, the high energy one
                       ]


    # --- controlling the noise floor ---
    countNoiseLevel = 2

    # --- accelerating potential toggles ---
    engy_Thresh = 120  # minimum allowable energy of the inverted-V potential
    maxfev = int(1E4) # number of iterations the LM fit is allowed
    useNoGuess = True # use an initial guess?

    # --- Levenberg-Marquart Fit toggles ---
    wPitchsToFit = [10, 20, 30] # give pitch angles in degrees
    wDistributionToFit = 'Maxwellian' # 'Maxwellian' or 'Kappa'
    numToAverageOver = 1 # HOW many datapoints are averaged together when fitting

    # Determine guesses for the fitted data
    V0_deviation = 0.18
    n_bounds = [0.001, 10]  # n [cm^-3]
    Te_bounds = [10, 300]
    kappa_bounds = [1.5, 101]

    if wDistributionToFit == 'Maxwellian':
        n_guess = 1
        T_guess = 100
        # can't do V0 guess, that's generated in the code itself
    elif wDistributionToFit == 'Kappa':
        n_guess = 1
        T_guess = 100
        kappa_guess = 20

    # --- fit refinement ---
    useFitRefinement = False
    beta_guess = 6 # altitude of the inverted-V
    n0guess_deviation = 0.99

class primaryBeamPlottingToggles:

    # -- Fit Statistics Toggles ---
    chiSquare_ThreshRange = [0.1, 100]  # range that the ChiSquare must fall into in order to be counted