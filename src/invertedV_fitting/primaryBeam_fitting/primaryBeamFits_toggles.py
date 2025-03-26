
##########################
# --- PRIMARY BEAM FIT ---
##########################
class primaryBeamToggles:

    # denote where the data should be taken from
    inputDataPath = GenToggles.input_diffNFiles[GenToggles.wFlyerFit]
    outputFolder = 'C:\Data\physicsModels\invertedV\primaryBeam_Fitting'

    # --- controlling the noise floor ---
    countNoiseLevel = 2

    # --- accelerating potential toggles ---
    engy_Thresh = 120  # minimum allowable energy of the inverted-V potential
    maxfev = int(1E4) # number of iterations the LM fit is allowed
    useNoGuess = True # use an initial guess?

    # --- Levenberg-Marquart Fit toggles ---
    wPitchsToFit = [10, 20, 30] # give pitch angles in degrees
    wDistributionToFit = 'Maxwellian' # 'Maxwellian' or 'Kappa'
    numToAverageOver = 3 # HOW many datapoints are averaged together when fitting

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