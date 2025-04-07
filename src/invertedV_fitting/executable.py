# --- executable.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: performs fits on diffNFlux data for ACES-II, but can be modified for other fluxes


#################
# --- IMPORTS ---
#################
import time
from fileIO_toggles import *
import spaceToolsLib as stl
start_time = time.time()


#################
# --- TOGGLES ---
#################
primaryBeam_fitting = False
primaryBeam_individualPlots = False
primaryBeam_fitParamPlots = False
backScatter_Calc = True
backScatter_Plotting = False


################################
# --- --- --- --- --- --- --- --
# --- ENVIRONMENT GENERATORS ---
# --- --- --- --- --- --- --- --
################################

if primaryBeam_fitting:
    stl.prgMsg('Generating Primary Beam Fit Parameters\n')
    from src.invertedV_fitting.primaryBeam_fitting.primaryBeamFits_generator import generatePrimaryBeamFit
    generatePrimaryBeamFit()
    stl.Done(start_time)

if primaryBeam_individualPlots:
    stl.prgMsg('Plotting Primary Beam Fits\n')
    from src.physicsModels.invertedV_fitting.primaryBeam_fitting.primaryBeamFits_Plotting import generatePrimaryBeamFitPlots
    generatePrimaryBeamFitPlots()
    stl.Done(start_time)

if primaryBeam_fitParamPlots:
    stl.prgMsg('Plotting Primary Beam Parameters\n')
    from src.physicsModels.invertedV_fitting.primaryBeam_fitting.primaryBeamFits_Plotting import generatePrimaryBeamFitPlots
    generatePrimaryBeamFitPlots()
    stl.Done(start_time)

if backScatter_Calc:
    stl.prgMsg('Generating Secondary/Backscatter Data\n')
    from src.physicsModels.invertedV_fitting.backScatter.backScatter_Generator import generateSecondaryBackScatter
    generateSecondaryBackScatter()
    stl.Done(start_time)

if backScatter_Plotting:
    stl.prgMsg('Plotting backScatter Beam Fits')
    from src.physicsModels.invertedV_fitting.backScatter.plotting.backScatter_Plotting import generateBackScatterPlots
    generateBackScatterPlots()
    stl.Done(start_time)
