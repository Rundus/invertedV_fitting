# --- model_primaryBeam_classes --
import spaceToolsLib as stl
from scipy.special import gamma
import numpy as np
from src.invertedV_fitting.user_toggles import PrimaryBeamFitToggles

class PrimaryBeamClasses:

    def calc_noise_level(self, counts_level, energy_value, pitch_angles):

        # (1) average the geofactors over the pitch range
        avg_indicies = [i for i in range(len(pitch_angles)) if pitch_angles[i] in PrimaryBeamFitToggles.pitch_angles_to_fit]
        geo_factor_avg = np.mean(np.array(PrimaryBeamFitToggles.geoFactor)[avg_indicies])

        return counts_level/(geo_factor_avg*(PrimaryBeamFitToggles.integration_time - PrimaryBeamFitToggles.deadtime*counts_level)*energy_value)

    def calc_jN_error(self, counts_val,energy_value, pitch_angles):
        # (1) average the geofactors over the pitch range
        avg_indicies = [i for i in range(len(pitch_angles)) if pitch_angles[i] in PrimaryBeamFitToggles.pitch_angles_to_fit]
        geo_factor_avg = np.mean(np.array(PrimaryBeamFitToggles.geoFactor)[avg_indicies])

        return np.sqrt(counts_val) / (geo_factor_avg * PrimaryBeamFitToggles.integration_time * energy_value)

    def form_fit_params(self, phi0_guess, **kwargs):

        # form the guesses list
        guesses = [PrimaryBeamFitToggles.n0_guess, PrimaryBeamFitToggles.T0_guess, phi0_guess]

        # form the boundaries dictionary
        fit_param_boundaries = [PrimaryBeamFitToggles.n_bounds, PrimaryBeamFitToggles.Te_bounds, [(1 - PrimaryBeamFitToggles.phi0_deviation) * phi0_guess, (1 + PrimaryBeamFitToggles.phi0_deviation) * phi0_guess]]

        if PrimaryBeamFitToggles.fit_dist == 'kappa':
            guesses += [PrimaryBeamFitToggles.kappa0_guess]
            fit_param_boundaries += [PrimaryBeamFitToggles.kappa_bounds]

        fit_param_boundaries = tuple(np.array(fit_param_boundaries).T)

        # determine the fitting function
        fit_func = self.diffNFlux_fitFunc_Kappa if PrimaryBeamFitToggles.fit_dist == 'kappa' else self.diffNFlux_fitFunc_Maxwellian

        # form the fitting parameters
        kwargs_dict = {
            'maxfev': PrimaryBeamFitToggles.maxfev,
            'bounds':fit_param_boundaries
        }

        if PrimaryBeamFitToggles.use_guess_bool:
            kwargs_dict['p0'] = guesses

        return fit_func, kwargs_dict

    # --- FUNCTION for fitting ---
    def diffNFlux_fitFunc_Maxwellian(self, x, n, T, V):  # Used in primary_beam_fit
        '''
        :param x: scalar energy on the BEAM energy grid [eV]
        :param n: plasma density [cm^-3]
        :param T: electron temperature [eV]
        :param V: inverted-V parallel potential [eV]
        :return:
        jN for maxwellian
        '''

        Energy = (x - V)

        # Create the Distribution function in m^-6s^3
        Dist = (1E6 * n) * np.power(stl.m_e / (2 * np.pi * stl.q0 * T), 3 / 2) * np.exp((-Energy / T))

        # convert to diffNFlux in m^-2J^-1sr^-1s^1
        diffNFlux = (2*stl.q0*x/np.power(stl.m_e,2))*Dist

        # convert cm^-2eV^-1
        diffNFlux_converted = (stl.q0 / np.power(stl.cm_to_m, 2)) * diffNFlux

        return diffNFlux_converted

    def diffNFlux_fitFunc_Kappa(self, x, n, T, V, kappa):  # Used in primary_beam_fit
        '''
        :param x: scalar - energy on the BEAM energy grid [eV]
        :param n: scalar - plasma density [cm^-3]
        :param T: scalar - electron temperature [eV]
        :param V: scalar - inverted-V parallel potential [eV]
        :param kappa: scalar - kappa function value
        :return:
        jN for kappa
        '''
        # Input energy  (in eV)
        Energy = (x - V)

        # Kappa Ek
        Ek = T*(1 - 3/(2*kappa))

        # create the Distribution function in m^-6s^3
        Dist = ((1E6)*n * np.power(stl.m_e/(2*np.pi*kappa*stl.q0*Ek),3/2) * (gamma(kappa+1)/gamma(kappa-0.5)) * np.power(1 + Energy/(kappa*Ek),-(kappa +1)))

        # convert to diffNFlux in m^-2J^-1sr^-1s^1
        diffNFlux = (2*stl.q0*x/np.power(stl.m_e,2))*Dist

        # convert cm^-2 eV^-1
        diffNFlux_converted = (stl.q0/np.power(stl.cm_to_m,2))*diffNFlux
        return diffNFlux_converted



