# --- model_primaryBeam_classes --
import spaceToolsLib as stl
import numpy as np
from src.invertedV_fitting.fit_primary_beam.primary_beam_fit_toggles import PrimaryBeamToggles
from scipy.special import gamma


class PrimaryBeamClasses:

    def form_fit_params(self, phi0_guess):

        # form the guesses list
        guesses = [PrimaryBeamToggles.n0_guess, PrimaryBeamToggles.T0_guess, phi0_guess]

        # form the boundaries dictionary
        fit_param_boundaries = [PrimaryBeamToggles.n_bounds, PrimaryBeamToggles.Te_bounds, [(1 - PrimaryBeamToggles.phi0_deviation) * phi0_guess, (1 + PrimaryBeamToggles.phi0_deviation) * phi0_guess]]

        if PrimaryBeamToggles.fit_dist == 'kappa':
            guesses += [PrimaryBeamToggles.kappa0_guess]
            fit_param_boundaries += [PrimaryBeamToggles.kappa_bounds]

        fit_param_boundaries = tuple(np.array(fit_param_boundaries).T)

        # determine the fitting function
        fit_func = self.diffNFlux_fitFunc_Kappa if PrimaryBeamToggles.fit_dist == 'kappa' else self.diffNFlux_fitFunc_Maxwellian

        # form the fitting parameters
        kwargs_dict = {
            'maxfev': PrimaryBeamToggles.maxfev,
            'bounds':fit_param_boundaries
        }

        if PrimaryBeamToggles.use_guess_bool:
            kwargs_dict['p0'] = guesses

        return fit_func, kwargs_dict

    # --- FUNCTION for fitting ---
    def diffNFlux_fitFunc_Maxwellian(self, x, n, T, V):  # Used in fit_primary_beam
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

    def diffNFlux_fitFunc_Kappa(self, x, n, T, V, kappa):  # Used in fit_primary_beam
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



