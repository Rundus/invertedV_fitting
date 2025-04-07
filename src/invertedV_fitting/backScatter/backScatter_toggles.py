import numpy as np

class backScatterToggles:

    outputFolder = r'C:\Data\physicsModels\invertedV\backScatter'

    # --- ENERGY GRID ---
    N_energyGrid = 500
    model_energyGrid = np.logspace(1, np.log10(2000), N_energyGrid)

    # --- model parameters ---
    modelParametersPitchAngle = 10#[degrees] - which pitch angle to use for the "primary beam"

    # --- Calculating backScatter ---
    betaChoice = 20 # which beta value to pick i.e. the height above the rocket of the invertedV
    niterations_backscatter = 6  # number of iterations for the secondaries calculations. >19 iterations is TOO many