
class UserToggles:
    path_to_data_folder = '/home/connor/Data'

    # --- User File I/O ---
    mission_name = 'TRACERS'
    data_category = 'SATELLITES'
    data_year = '2025'
    data_month = '10'
    data_day = '18'
    payload_designator = 'ts2'
    path_to_eESA_data = f'{path_to_data_folder}/{data_category}/{mission_name}/{data_year}/{data_month}/{data_day}/{payload_designator}/ts2_l3_ace_pitch-angle-dist_20251018_v1.0.0.cdf'

    # data keys for this mission - TRACERS
    differential_number_flux_key = None
    differential_energy_flux_key = 'def'
    Epoch_key = 'Epoch'
    pitch_angle_key = 'ts2_l3_ace_pitch_angle'
    energy_key = 'ts2_l3_ace_energy'
    dependency_structure = [Epoch_key, energy_key, pitch_angle_key]

    # --- Fitting toggles ---
    pitch_angles_to_fit = [5,15,25] # in degrees

    # --- Program File I/O ---
    run_folder_path = f'{path_to_data_folder}/MODELS/invertedV_fitting/{mission_name}/{data_year}/{data_month}/{data_day}/{payload_designator}/'
