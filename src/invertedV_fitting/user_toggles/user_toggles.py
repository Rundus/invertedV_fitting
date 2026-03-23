
import datetime as dt

class UserToggles:
    path_to_data_folder = '/home/connor/Data'

    # --- User File I/O ---
    mission_name = 'TRACERS'
    instr_name = 'ACE'
    data_category = 'SATELLITES'
    data_year = '2025'
    data_month = '10'
    data_day = '18'
    payload_designator = 'ts2'
    path_to_eESA_flux_data = f'{path_to_data_folder}/{data_category}/{mission_name}/{instr_name}/{data_year}/{data_month}/{data_day}/{payload_designator}/ts2_l3_ace_pitch-angle-dist_20251018_v1.0.0.cdf'
    path_to_eESA_counts_data = f'{path_to_data_folder}/{data_category}/{mission_name}/{instr_name}/{data_year}/{data_month}/{data_day}/{payload_designator}/ts2_l3_ace_pitch-angle-dist_20251018_v1.0.0.cdf'

    # --- Fitted Region ---
    datetime_low = dt.datetime(2025, 10, 18, 1, 00, 00)
    datetime_high = dt.datetime(2025, 10, 18, 1, 12, 00)

    # --- DATA STRUCTURE ---

    # data keys for this mission - TRACERS
    differential_number_flux_key = None
    differential_energy_flux_key = 'ts2_l3_ace_pitch_def'
    counts_key = 'ts2_l3_ace_pitch_background_counts'
    Epoch_key = 'Epoch'
    pitch_angle_key = 'ts2_l3_ace_pitch_angle'
    energy_key = 'ts2_l3_ace_energy'
    energy_bin_direction = 0 # 0 - energy bins decrease in value, 1 - energy bin increases in value
    dependency_structure = [Epoch_key, energy_key, pitch_angle_key]
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
    deadtime = 80E-9 # in seconds
    integration_time = 0.9E-3 # in seconds.

    # --- Program File I/O ---
    run_folder_path = f'{path_to_data_folder}/MODELS/invertedV_fitting/{mission_name}/{instr_name}/{data_year}/{data_month}/{data_day}/{payload_designator}/'
