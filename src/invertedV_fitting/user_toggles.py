# --- Imports ---
import numpy as np
import datetime as dt


class FileToggles:

    # --- User I/O Info ---
    mission_name = 'TRACERS'
    instr_name = 'ACE'
    data_year = '2026'
    data_month = '03'
    data_day = '1'
    data_level = 'l3'
    payload_designator = '2' #

    # --- Program File I/O ---
    PATH_TO_FOLDER = '/home/connor/Data/MODELS/invertedV_fitting/' # Path to folder where everything is saved
    RUN_PATH = f'{PATH_TO_FOLDER}/{mission_name}/{instr_name}/{data_year}/{data_month}/{data_day}/ts{payload_designator}/'


class FitDataToggles:

    # data keys for this mission - TRACERS
    differential_number_flux_key = None
    differential_energy_flux_key = 'ts2_l3_ace_pitch_def'
    counts_key = 'ts2_l3_ace_pitch_background_counts'
    epoch_key = 'Epoch'
    pitch_angle_key = 'ts2_l3_ace_pitch_angle'
    energy_key = 'ts2_l3_ace_energy'
    energy_bin_direction = 0  # 0 - energy bins decrease in value, 1 - energy bin increases in value
    dependency_structure = [epoch_key, energy_key, pitch_angle_key]
    geoFactor = [1.067e-04,
                 1.264e-04,
                 1.248e-04,
                 1.221e-04,
                 1.200e-04,
                 1.276e-04,
                 7.457e-05,
                 7.963e-05,
                 1.218e-04,
                 1.267e-04,
                 1.299e-04,
                 1.299e-04,
                 1.265e-04,
                 1.252e-04,
                 1.202e-04,
                 7.571e-05,
                 1.103e-04,
                 1.276e-04,
                 1.533e-04,
                 1.920e-04,
                 1.335e-04]
    deadtime = 80E-9  # in seconds
    integration_time = 0.9E-3  # in seconds.

class PrimaryBeamFitToggles:

    # --- Fit Data Region Definition ---
    # window of data which is fitted. Anything outside of this region is not fit
    start_hour = 16 # [UTC]
    start_minute = 0
    end_hour = 16# [UTC]
    end_minute = 20

    N_time_avg = 5 # number of time-slices to average together. Done sequentially.

    # --- FIT DATA COLLECTION ---
    pitch_angles_to_fit = [45, 55, 65]  # in degrees
    energy_thesh = 100  # in [eV]. The energy to start searching above for  electrostatic potentials

    # --- Inverted-V Fit Parameters (Marquart-Levenburg)
    # General Fit
    fit_dist = 'maxwellian' # options: 'kappa', 'maxwellian'
    countNoiseLevel = 2
    maxfev = int(1E3)  # number of iterations the LM fit is allowed

    # guesses
    use_guess_bool = True # use an initial guess for the fit
    n0_guess = 1 # [cm^-3] plasma density
    T0_guess = 100 # [eV] Plasma temperature
    kappa0_guess = 20, # kappa parameter. Only matters if fit_dist == 'kappa'

    phi0_deviation = 0.18 # to match the detector's resolution
    n_bounds = [0.001, 10]  # n [cm^-3]
    Te_bounds = [10, 300]
    kappa_bounds = [1.5, 101]

    # Refine Fit
    use_kaeppler_fit_refinement_bool = False
    beta_guess = 6 # altitude of the inverted-V
    n0guess_deviation = 0.99


class BackScatterToggles:

    # --- ENERGY GRID ---
    N_energyGrid = 500
    model_energyGrid = np.logspace(1, np.log10(2000), N_energyGrid)

    # --- model parameters ---
    modelParametersPitchAngle = 10#[degrees] - which pitch angle to use for the "primary beam"

    # --- Calculating backscatter ---
    betaChoice = 20 # which beta value to pick i.e. the height above the rocket of the invertedV
    niterations_backscatter = 6  # number of iterations for the secondaries calculations. >19 iterations is TOO many