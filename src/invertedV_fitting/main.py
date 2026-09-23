dict_executable = {
    'regen_EVERYTHING': 0,
    'primary_beam_fit': 1,
    'plot_primary_beam_fits': 0,
    'calc_backscatter': 0,
}


if __name__ == "__main__":
    from src.invertedV_fitting.runners.execute_subroutines import run_invertedV_fitting
    run_invertedV_fitting(dict_executable)