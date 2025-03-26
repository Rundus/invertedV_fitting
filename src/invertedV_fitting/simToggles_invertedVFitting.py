# --- imports ---
import numpy as np
from spaceToolsLib.variables import m_to_km
from datetime import datetime

######################
# ---GENERAL SETUP ---
######################
class GenToggles:
    wFlyerFit = 0
    input_diffNFiles = ['C:\Data\ACESII\L2\high\ACESII_36359_l2_eepaa_fullCal.cdf',
                        'C:\Data\ACESII\L2\low\ACESII_36364_l2_eepaa_fullCal.cdf']
    input_attitudeFiles = [r'C:\Data\ACESII\attitude\high\ACESII_36359_Attitude_Solution.cdf',
                           r'C:\Data\ACESII\attitude\low\ACESII_36364_Attitude_Solution.cdf']

    if wFlyerFit == 0: # ACES-II High Flyer Data
        invertedV_times = [
                            [datetime(2022, 11, 20, 17, 25,  1, 000000), datetime(2022, 11, 20, 17, 25, 3, 000000)], # Dispersive Region
                            [datetime(2022, 11, 20, 17, 24, 12, 162000), datetime(2022, 11, 20, 17, 24, 18, 812000)], # Very First ,Inverted-V, the high energy one
                            [datetime(2022, 11, 20, 17, 24, 45, 862000), datetime(2022, 11, 20, 17, 24, 49, 312000)], # small inverted-V, after the High energy One
                            [datetime(2022, 11, 20, 17, 25, 23, 762000), datetime(2022, 11, 20, 17, 26, 8, 212000)],  # Primary inverted-V
                            [datetime(2022, 11, 20, 17, 26, 11, 412000), datetime(2022, 11, 20, 17, 26, 19, 912000)],  # Inverted-V right after the Primary-V, has STEBs on either sides of it
                            [datetime(2022, 11, 20, 17, 26, 35, 112000), datetime(2022, 11, 20, 17, 26, 40, 712000)], # Inverted-V two after the Primary-V
                            [datetime(2022, 11, 20, 17, 28, 17, 112000), datetime(2022, 11, 20, 17, 28, 34, 612000)] # Faint inverted-V on the most northside of the flight
                           ]
    elif wFlyerFit == 1: # ACES-II Low Flyer Data
        invertedV_times = [
                            [datetime(2009, 1, 29, 9, 54, 4, 0), datetime(2009, 1, 29, 9, 54, 29, 000)] # Very First ,Inverted-V, the high energy one
                            ]
















