def run_invertedV_fitting(dict_executable):

    #################
    # --- IMPORTS ---
    #################
    from src.invertedV_fitting.runners.executable_classes import ExecutableClasses

    # run overhead
    exec_obj = ExecutableClasses()
    exec_obj.generate_run_directories(dict_executable) # Verify the directories for this run exist
    exec_obj.load_fit_data()
    exec_obj.generate_run_JSON() # Generate the Configuration File for this run

    # Run the subroutines
    if dict_executable['regen_EVERYTHING']==1:
        for key in dict_executable.keys():
            dict_executable[key] = 1

    if dict_executable['primary_beam_fit']==1:
        print('\n--- Fitting Primary Inverted-V Beam ---', end='\n')
        from src.invertedV_fitting.primary_beam_fit.primary_beam_fit_generator import primary_beam_fit_generator
        primary_beam_fit_generator()

    if dict_executable['plot_primary_beam_fits']==1:
        print('\n--- PLOTTING Primary Inverted-V Beam ---', end='\n')
        from src.invertedV_fitting.plotting.primary_beam_plots.plot_primary_beam_fits_generator import plot_primary_beam_fits_generator
        plot_primary_beam_fits_generator()

    if dict_executable['calc_backscatter']==1:
        print('\n--- Calculating Beam Backscatter ---', end='\n')
        from src.invertedV_fitting.backscatter.backscatter_generator import generateSecondaryBackScatter
        generateSecondaryBackScatter()
