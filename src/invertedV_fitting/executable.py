# --- executable.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: performs fits on differential_number_flux data.


#################
# --- IMPORTS ---
#################
from src.invertedV_fitting.executable_classes import ExecutableClasses
from src.invertedV_fitting.executable_toggles import dict_executable

# Generate the Configuration File for this run
ExecutableClasses().generate_run_JSON()

if dict_executable['regen_EVERYTHING']==1:
    for key in dict_executable.keys():
        dict_executable[key] = 1

if dict_executable['fit_primary_beam']==1:
    print('\n--- Fitting Primary Inverted-V Beam ---', end='\n')
    from src.invertedV_fitting.fit_primary_beam.primary_beam_fit_generator import primary_beam_fit_generator
    primary_beam_fit_generator()


if dict_executable['calculate_backscatter']==1:
    print('\n--- Calculating Beam Backscatter ---', end='\n')
    from src.invertedV_fitting.backscatter.backscatter_generator import generateSecondaryBackScatter
    generateSecondaryBackScatter()
