from src.invertedV_fitting.user_toggles.user_toggles import UserToggles
from src.invertedV_fitting.primary_beam_fit.primary_beam_fit_toggles import PrimaryBeamToggles
import numpy as np

class PlotPrimaryBeamClasses:

    def calc_noise_level(self, counts_level, energy_value, pitch_angles):

        # (1) average the geofactors over the pitch range
        avg_indicies = [i for i in range(len(pitch_angles)) if pitch_angles[i] in PrimaryBeamToggles.pitch_angles_to_fit]
        geo_factor_avg = np.mean(np.array(UserToggles.geoFactor)[avg_indicies])

        return counts_level/(geo_factor_avg*UserToggles.integration_time*energy_value)