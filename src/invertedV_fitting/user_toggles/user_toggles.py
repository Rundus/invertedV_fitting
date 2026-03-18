
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
    path_to_eESA_data = f'{path_to_data_folder}/{data_category}/{mission_name}/{instr_name}/{data_year}/{data_month}/{data_day}/{payload_designator}/ts2_l3_ace_pitch-angle-dist_20251018_v1.0.0.cdf'

    # data keys for this mission - TRACERS
    differential_number_flux_key = None
    differential_energy_flux_key = 'ts2_l3_ace_pitch_def'
    Epoch_key = 'Epoch'
    pitch_angle_key = 'ts2_l3_ace_pitch_angle'
    energy_key = 'ts2_l3_ace_energy'
    energy_bin_direction = 0 # 0 - energy bins decrease in value, 1 - energy bin increases in value
    dependency_structure = [Epoch_key, energy_key, pitch_angle_key]
    datetime_low = dt.datetime(2025,10,18,1,00,00)
    datetime_high = dt.datetime(2025,10, 18, 1, 12,00)

    # --- Fitting toggles ---
    pitch_angles_to_fit = [5,15,25] # in degrees
    energy_thesh = 100 # in [eV]. The energy to start searching above for  electrostatic potentials

    # --- Program File I/O ---
    run_folder_path = f'{path_to_data_folder}/MODELS/invertedV_fitting/{mission_name}/{instr_name}/{data_year}/{data_month}/{data_day}/{payload_designator}/'
