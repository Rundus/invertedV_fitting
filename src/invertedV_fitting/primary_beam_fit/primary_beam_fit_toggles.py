

class PrimaryBeamToggles:

    # --- FIT DATA COLLECTION ---
    pitch_angles_to_fit = [5, 15, 25]  # in degrees
    energy_thesh = 100  # in [eV]. The energy to start searching above for  electrostatic potentials

    #########################
    # --- FITTING TOGGLES ---
    #########################

    # General Fit
    fit_dist = 'maxwellian' # options: 'kappa', 'maxwellian'
    countNoiseLevel = 2
    maxfev = int(1E3)  # number of iterations the LM fit is allowed

    # guesses
    use_guess_bool = True # use an initial guess for the fit
    n0_guess = 1 # plasma density [cm^-3]
    T0_guess = 100 # Plasma temperature [eV]
    kappa0_guess = 20, # kappa parameter

    phi0_deviation = 0.18 # to match the detector's resolution
    n_bounds = [0.001, 10]  # n [cm^-3]
    Te_bounds = [10, 300]
    kappa_bounds = [1.5, 101]

    # Refine Fit
    use_kaeppler_fit_refinement_bool = False
    beta_guess = 6 # altitude of the inverted-V
    n0guess_deviation = 0.99
