# --- Imports ---
import numpy as np
import datetime as dt


class FileToggles:

    # --- User I/O Info ---
    mission_name = 'TRACERS'  # mission name. Used to build RUN_PATH
    instr_name = 'ACE'  # instrument name. Used to build RUN_PATH
    data_year = '2026'  # [YYYY] year of the data
    data_month = '04'  # [MM] month of the data
    data_day = '01'  # [DD] day of the month of the data
    data_level = 'l3'  # [l#] data processing level of the input file
    payload_designator = '2'  # TRACERS spacecraft number (1 or 2), i.e. ts1 or ts2

    # --- Program File I/O ---
    PATH_TO_FOLDER = '/home/connor/Data/MODELS/invertedV_fitting/'  # Path to folder where everything is saved
    # run-specific folder under PATH_TO_FOLDER, organized by mission/instrument/date/spacecraft
    RUN_PATH = f'{PATH_TO_FOLDER}/{mission_name}/{instr_name}/{data_year}/{data_month}/{data_day}/ts{payload_designator}/'


class FitDataToggles:

    # --- Data Keys (TRACERS) ---
    differential_number_flux_key = None  # key for the differential number flux. None if not present in the data file
    differential_energy_flux_key = 'ts2_l3_ace_pitch_def'  # key for the differential energy flux
    counts_key = 'ts2_l3_ace_pitch_background_counts'  # key for the background counts
    epoch_key = 'Epoch'  # key for the measurement timestamps
    pitch_angle_key = 'ts2_l3_ace_pitch_angle'  # key for the pitch-angle bin centers
    energy_key = 'ts2_l3_ace_energy'  # key for the energy bin centers
    energy_bin_direction = 0  # 0 - energy bins decrease in value, 1 - energy bins increase in value
    dependency_structure = [epoch_key, energy_key, pitch_angle_key]  # order of the data array dimensions: [time, energy, pitch angle]

    # --- Instrument Characteristics ---
    # [cm^2 sr eV/eV] geometric factor for each pitch-angle bin, ordered to match pitch_angle_key
    geoFactor = np.array([1.067e-04,
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
                 1.335e-04])
    deadtime = 80E-9  # [s] detector dead time
    integration_time = 0.9E-3  # [s] accumulation time for each energy bin


class PrimaryBeamFitToggles:

    # --- Fit Data Region Definition ---
    # window of data which is fitted. Anything outside of this region is not fit
    start_hour = '03'  # [UTC] hour at which the fit window starts
    start_minute = '08'  # [UTC] minute at which the fit window starts
    start_second = '59'
    end_hour = '03'  # [UTC] hour at which the fit window ends
    end_minute = '09'  # [UTC] minute at which the fit window ends
    end_second = '01'
    N_time_avg = 1  # number of time-slices to average together. Done sequentially.

    # --- FIT DATA COLLECTION ---
    pitch_angles_to_fit = [[5],[15],[25]]  # [degrees], If multiple pitch angles are to be averaged together put them in brackets, e.g. [[pitch1, pitch2,pitch3], [pitch,4],...] will average pitches 1 to 3, then fit whereas [pitch4] will only fit pitch 4.
    energy_thesh = 100  # in [eV]. The energy to start searching above for electrostatic potentials

    # --- Inverted-V Fit Parameters (Marquardt-Levenberg) ---
    # General Fit
    fit_dist = 'maxwellian'  # options: 'kappa', 'maxwellian'. Distribution function fit to the primary beam
    countNoiseLevel = 2  # [counts] bins at or below this count level are treated as noise and excluded from the fit
    maxfev = int(1E3)  # maximum number of function evaluations the LM fit is allowed
    N_fit_min = 3 # MINIMUM number of points each fit must have

    # guesses
    use_guess_bool = True  # use an initial guess for the fit
    n0_guess = 1  # [cm^-3] plasma density
    T0_guess = 100  # [eV] Plasma temperature
    kappa0_guess = 20,  # kappa parameter. Only matters if fit_dist == 'kappa'

    # bounds
    phi0_deviation = 0.18  # fractional range the fitted potential may deviate from its initial guess. 0.18 matches the detector's dE/E resolution
    n_bounds = [0.001, 10]  # [cm^-3] Bounds on the plasma density
    Te_bounds = [10, 300]  # [eV] Bounds on the plasma temperature
    kappa_bounds = [1.5, 101]  # Bounds on the kappa parameter. Only matters if fit_dist == 'kappa'

    # Refine Fit
    use_kaeppler_fit_refinement_bool = False  # run a second fit (Kaeppler et al. method) that refines the General Fit result
    beta_guess = 6  # initial guess for beta, the magnetic mirror ratio that sets the altitude of the inverted-V
    n0guess_deviation = 0.99  # fractional range n0 may deviate from the General Fit result during refinement


class BackScatterToggles:

    # --- ENERGY GRID ---
    N_energyGrid = 500  # number of points in the model energy grid
    model_energyGrid = np.logspace(1, np.log10(2000), N_energyGrid)  # [eV] log-spaced model energy grid from 10 eV to 2 keV

    # --- model parameters ---
    modelParametersPitchAngle = 10  # [degrees] - which pitch angle to use for the "primary beam"

    # --- Calculating backscatter ---
    betaChoice = 20  # which beta value to pick i.e. the height above the rocket of the invertedV
    niterations_backscatter = 6  # number of iterations for the secondaries calculations. >19 iterations is TOO many