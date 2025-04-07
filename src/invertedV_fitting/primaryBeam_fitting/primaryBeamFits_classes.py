# --- model_primaryBeam_classes --
import spaceToolsLib as stl
import numpy as np
from scipy.special import gamma
from copy import deepcopy

class helperFuncs:
    def distFunc_to_diffNFlux(self, Vperp, Vpara, dist, mass, charge):
        # Input: Velocities [m/s], distribution function [s^3m^-6]
        # output: diffNFlux [cm^-2 s^-1 eV^-1 str^-1]
        Emag = 0.5 * mass * (Vperp ** 2 + Vpara ** 2) / charge
        return (2 * Emag) * np.power(charge / (100 * mass), 2) * dist
    def diffNFlux_to_distFunc(self, Vperp, Vpara, diffNFlux, mass, charge):
        # Input: Vperp,Vpara in [m/s]. DiffNFlux in [cm^-2 s^-1 str^-1 eV^-1]
        # output: distribution function [s^3m^-6]
        Energy = 0.5 * mass * (Vperp ** 2 + Vpara ** 2) / charge
        return 0.5 * np.power((100 * mass / charge), 2) * diffNFlux / Energy
    def generateNoiseLevel(self, energyData, countNoiseLevel):
        count_interval = 0.8992E-3
        geo_factor = 8.63E-5
        deadtime = 324E-9

        # --- DEFINE THE NOISE LEVEL ---
        diffNFlux_NoiseCount = np.zeros(shape=(len(energyData)))

        for idx,engy in enumerate(energyData):
            deltaT = (count_interval) - (countNoiseLevel * deadtime)
            diffNFlux_NoiseCount[idx] = (countNoiseLevel) / (geo_factor * deltaT * engy)

        return diffNFlux_NoiseCount
    def groupAverageData(self, data_dict_diffFlux, targetTimes, N_avg, **kwargs):
        '''
        # Input:
        # data_dict_diffFlux - data_dict
            "Differential_Number_Flux", "Pitch_Angle" and "Energy"  variables
        # targetTimes - 1D array
            [T_min, T_max] where T values are datetimes of the min/max region of the data to average
        # N_avg - scalar
            the number of Epoch values to average over. Should be odd
        # fluxType (kwarg) - string.
            Option of either "diffNFlux" or "diffEFlux" to determine if number flux or energy flux should be used

        # Output:
        Epoch, Differential_Number_Flux and stdDevs averaged over N TIME-points specified in the primaryBeamToggles. Dimensions in Pitch and Energy are untouched
        '''

        if N_avg%2 == 0:
            raise Exception('Number of Points to Average over must be odd')

        # Determine if the output is differential energy or number flux
        if kwargs.get('fluxType',None) == 'diffEFlux':
            data = deepcopy(data_dict_diffFlux['Differential_Energy_Flux'][0])
            data_stdDev = np.multiply( deepcopy(data_dict_diffFlux['Differential_Number_Flux_stdDev'][0]), data_dict_diffFlux['Energy'][0] ,axis=2)
        else:
            data = deepcopy(data_dict_diffFlux['Differential_Number_Flux'][0])
            data_stdDev = deepcopy(data_dict_diffFlux['Differential_Number_Flux_stdDev'][0])



        ##############################
        # --- COLLECT THE FIT DATA ---
        ##############################
        # ensure the data is divided into chunks that can be sub-divided. If not, keep drop points from the end until it can be
        low_idx, high_idx = np.abs(data_dict_diffFlux['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_diffFlux['Epoch'][0] - targetTimes[1]).argmin()

        if (high_idx - low_idx) % N_avg != 0:
            high_idx -= (high_idx - low_idx) % N_avg

        # Handle the Epoch
        chunkedEpoch = np.split(data_dict_diffFlux['Epoch'][0][low_idx:high_idx], round(len(data_dict_diffFlux['Epoch'][0][low_idx:high_idx]) / N_avg))
        EpochFitData = np.array([chunkedEpoch[i][int((N_avg - 1) / 2)] for i in range(len(chunkedEpoch))])
        chunkedIlat = np.split(data_dict_diffFlux['ILat'][0][low_idx:high_idx], round(len(data_dict_diffFlux['ILat'][0][low_idx:high_idx]) / N_avg))
        ILatFitData = np.array([chunkedIlat[i][int((N_avg - 1) / 2)] for i in range(len(chunkedIlat))])
        chunkedAlt = np.split(data_dict_diffFlux['Alt'][0][low_idx:high_idx],round(len(data_dict_diffFlux['Alt'][0][low_idx:high_idx]) / N_avg))
        AltFitData = np.array([chunkedIlat[i][int((N_avg - 1) / 2)] for i in range(len(chunkedAlt))])

        # --- handle the multi-dimenional data ---

        # create the storage variable
        detectorPitchAngles = data_dict_diffFlux['Pitch_Angle'][0]
        diffFlux_avg = np.zeros(shape=(len(EpochFitData), len(detectorPitchAngles), len(data_dict_diffFlux['Energy'][0])))
        stdDevs_avg = np.zeros(shape=(len(EpochFitData), len(detectorPitchAngles), len(data_dict_diffFlux['Energy'][0])))

        for loopIdx, pitchValue in enumerate(detectorPitchAngles):

            chunkedyData = np.split(data[low_idx:high_idx, loopIdx, :], round(len(data[low_idx:high_idx, loopIdx,:]) / N_avg))
            chunkedStdDevs = np.split(data_stdDev[low_idx:high_idx, loopIdx, :], round(len(data_stdDev[low_idx:high_idx, loopIdx,:]) / N_avg))

            # --- Average the chunked data ---
            fitData = np.zeros(shape=(len(chunkedyData), len(data_dict_diffFlux['Energy'][0])))
            fitData_stdDev = np.zeros(shape=(len(chunkedStdDevs), len(data_dict_diffFlux['Energy'][0])))

            for i in range(len(chunkedEpoch)):
                # average the diffFlux data by only choosing data which is valid
                chunkedyData[i][chunkedyData[i] < 0] = np.NaN
                fitData[i] = np.nanmean(chunkedyData[i], axis=0)

                # average the diffFlux data by only choosing data which is valid
                chunkedStdDevs[i][chunkedStdDevs[i] < 0] = np.NaN
                fitData_stdDev[i] = np.nanmean(chunkedStdDevs[i], axis=0)

            diffFlux_avg[:, loopIdx, :] = fitData
            stdDevs_avg[:, loopIdx, :] = fitData_stdDev

        return EpochFitData, ILatFitData, AltFitData, diffFlux_avg, stdDevs_avg
    def removeDuplicates(self, a, b):
        from collections import defaultdict
        D = defaultdict(list)
        for i, item in enumerate(a):
            D[item].append(i)
        D = {k: v for k, v in D.items() if len(v) > 1}
        badIndicies = [D[key][1:] for key in D.keys()]
        badIndicies = [item for sublist in badIndicies for item in sublist]
        newA = np.delete(a, badIndicies, axis=0)
        newB = np.delete(b, badIndicies, axis=0)
        return newA, newB

class distributions_class:

    def generate_Maxwellian_Espace(self, n, T, energy_Grid):
        '''
        :param n: density [cm^-3]
        :param T: Temperature [eV]
        :param energy_Grid: Energy Grid [eV]
        :return: distribution function in SI units [s^3 m^-6]
        '''
        return (1E6 * n) * np.power(stl.m_e / (2 * np.pi * stl.q0*T), 3 / 2) * np.exp(-1 * energy_Grid / T)

    def calc_diffNFlux_Espace(self, dist, energy_Grid):
        '''
        :param dist: distribution function in SI units [s^3 m^-6]
        :param energy_Grid: Energy Grid [eV]
        :return: differential Number Flux [cm^-2s^-1sr^-1eV^-1]
        '''

        # in SI units
        diffNFlux = (2*stl.q0*energy_Grid/np.power(stl.m_e,2))*dist

        # in cm^-2 eV^-1
        diffNFlux_converted = (stl.q0/np.power(stl.cm_to_m,2)) * diffNFlux

        return diffNFlux_converted

    # --- Generate Distributions from Velocity Space ---
    def generate_Maxwellian_Vspace(self, mass, charge, n, T, Vperp, Vpara):
        # Input: density [cm^-3], Temperature [eV], Velocities [m/s]
        # output: the distribution function in SI units [s^3 m^-6]
        Emag = (0.5 * mass * (Vperp ** 2 + Vpara ** 2)) / charge
        return (1E6 * n) * np.power(mass / (2 * np.pi * charge * T), 3 / 2) * np.exp(-1 * Emag / T)

    def generate_kappa_Vspace(self, mass, charge, n, T, kappa, Vperp, Vpara):
        # Input: density [cm^-3], Temperature [eV], Velocities [m/s]
        # output: the distribution function in SI units [s^3 m^-6]
        Emag = (0.5 * mass * (Vperp ** 2 + Vpara ** 2)) / charge
        Ek = T*(1 - 3/(2*kappa))
        return (1E6)*n * np.power(mass/(2*np.pi*kappa*stl.q0*Ek),3/2) * (gamma(kappa+1)/gamma(kappa-0.5)) * np.power(1 + Emag/(kappa*Ek),-(kappa +1))
class primaryBeam_class:

    # --- FUNCTION for fitting ---
    def diffNFlux_fitFunc_Maxwellian(self, x, n, T, V):  # Used in primaryBeam_fitting
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

    def diffNFlux_fitFunc_Kappa(self, x, n, T, V, kappa):  # Used in primaryBeam_fitting
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



