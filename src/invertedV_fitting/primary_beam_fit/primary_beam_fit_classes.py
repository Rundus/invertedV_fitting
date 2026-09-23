# --- model_primaryBeam_classes --
import spaceToolsLib as stl
from scipy.special import gamma
import numpy as np
from src.invertedV_fitting.user_toggles import PrimaryBeamFitToggles
from src.invertedV_fitting.user_toggles import FitDataToggles
import warnings

class PrimaryBeamClasses:


    def __init__(self,epoch,energy,pitch_angle):
        self.epoch = epoch
        self.energy = energy
        self.pitch_angle = pitch_angle

    def standardize_flux_array(self, spectra, axis_hint=None):
        """
        Reorder a 3-D flux array to (epoch, energy, pitch) with standard bin directions.

        The axes of ``spectra`` are matched to the dependency arrays by length, and
        ``spectra`` is transposed into the canonical order (epoch, energy, pitch). The
        directions of the energy and pitch-angle arrays are then detected and
        standardized so that energy decreases (high to low) and pitch angle
        increases (low to high). Whenever a dependency array is reversed, the
        corresponding axis of ``spectra`` is reversed with it so the two remain aligned.

        Parameters
        ----------
        spectra : array_like
            3-D data array (e.g. differential flux) whose axes correspond to
            ``epoch``, ``energy`` and ``pitch`` in some unknown order.
        epoch : array_like
            1-D array of timestamps.
        energy : array_like
            1-D array of energy bin values. Must be strictly monotonic, in
            either direction.
        pitch : array_like
            1-D array of pitch-angle bin values. Must be strictly monotonic, in
            either direction.
        axis_hint : sequence of str, optional
            The axis order of ``spectra`` given as names, e.g.
            ``('epoch', 'pitch', 'energy')``. Only needed when two dependency
            arrays have the same length, in which case the order cannot be
            determined from shapes alone. If given, it is checked against the
            shapes.

        Returns
        -------
        spectra_out : numpy.ndarray
            ``spectra`` reordered to shape ``(len(epoch), len(energy), len(pitch))``,
            with the energy axis running high to low and the pitch axis running
            low to high. This is a view of ``spectra``, not a copy.
        energy_out : numpy.ndarray
            ``energy`` sorted high to low, aligned with axis 1 of ``spectra_out``.
        pitch_out : numpy.ndarray
            ``pitch`` sorted low to high, aligned with axis 2 of ``spectra_out``.
        info : dict
            ``'axis_order'``: tuple of names giving the detected axis order of
            the input ``spectra``.
            ``'energy_was_increasing'``: bool giving the detected direction of
            the input ``energy``.
            ``'pitch_was_increasing'``: bool giving the detected direction of
            the input ``pitch``.

        Raises
        ------
        ValueError
            If ``spectra`` is not 3-D, a dependency array is not 1-D, the shape of
            ``spectra`` does not match the dependency lengths, the axis order is
            ambiguous and no ``axis_hint`` is given, or ``energy`` or ``pitch``
            is not strictly monotonic.

        """

        import itertools
        import numpy as np

        CANONICAL_ORDER = ('epoch', 'energy', 'pitch')

        def _is_increasing(arr, name):
            """Return True if ``arr`` is strictly increasing, False if strictly decreasing."""
            steps = np.diff(arr)
            if np.all(steps > 0):
                return True
            if np.all(steps < 0):
                return False
            raise ValueError(f"'{name}' must be strictly increasing or strictly decreasing")

        spectra = np.asarray(spectra)
        deps = {'epoch': np.asarray(self.epoch),
                'energy': np.asarray(self.energy),
                'pitch': np.asarray(self.pitch_angle)}

        if spectra.ndim != 3:
            raise ValueError(f"spectra must be 3-D, got shape {spectra.shape}")
        for name, arr in deps.items():
            if arr.ndim != 1:
                raise ValueError(f"'{name}' must be 1-D, got shape {arr.shape}")

        # --- identify the axis order of spectra by matching lengths ---
        lengths = {name: arr.size for name, arr in deps.items()}
        lengths_str = ', '.join(f'{name}={n}' for name, n in lengths.items())
        candidates = [order for order in itertools.permutations(CANONICAL_ORDER)
                      if all(lengths[name] == spectra.shape[i] for i, name in enumerate(order))]

        if not candidates:
            raise ValueError(f"spectra has shape {spectra.shape}, which does not match the "
                             f"dependency lengths ({lengths_str})")

        if axis_hint is not None:
            axis_hint = tuple(axis_hint)
            if axis_hint not in candidates:
                raise ValueError(f"axis_hint {axis_hint} is inconsistent with spectra.shape "
                                 f"{spectra.shape} and dependency lengths ({lengths_str})")
            axis_order = axis_hint
        elif len(candidates) > 1:
            raise ValueError(f"Axis order is ambiguous for spectra.shape {spectra.shape} with "
                             f"dependency lengths ({lengths_str}). Possible orders: "
                             f"{candidates}. Pass axis_hint to choose one.")
        else:
            axis_order = candidates[0]

        # --- transpose spectra to (epoch, energy, pitch) ---
        spectra_out = spectra.transpose([axis_order.index(name) for name in CANONICAL_ORDER])

        # --- energy: standardize to high -> low ---
        energy_out = deps['energy']
        energy_was_increasing = _is_increasing(energy_out, 'energy')
        if energy_was_increasing:
            spectra_out = spectra_out[:, ::-1, :]
            energy_out = energy_out[::-1]

        # --- pitch angle: standardize to low -> high ---
        pitch_out = deps['pitch']
        pitch_was_increasing = _is_increasing(pitch_out, 'pitch')
        if not pitch_was_increasing:
            spectra_out = spectra_out[:, :, ::-1]
            pitch_out = pitch_out[::-1]

        info = {'axis_order': axis_order,
                'energy_was_increasing': energy_was_increasing,
                'pitch_was_increasing': pitch_was_increasing}
        return spectra_out, energy_out, pitch_out, info

    def pitch_average(self, spectra, groups, atol=1e-3):
        """
        Nan-average a 3-D array over groups of pitch angles.

        ``spectra`` is assumed to be in the standardized (epoch, energy, pitch) layout.
        Each entry of ``groups`` lists the pitch angles that are combined into one
        output bin: groups with several pitch angles are averaged with
        ``numpy.nanmean``, and groups with a single pitch angle are copied through
        unchanged. The pitch axis of the output has length ``len(groups)``, in the
        same order as ``groups``.

        Parameters
        ----------
        spectra : array_like
            3-D array of shape (N, M, Q), with pitch angle along axis 2.
        groups : sequence of sequences
            Pitch angles to combine, e.g. ``[[5, 15, 25], [35], [75]]``. Every
            pitch angle must match a value in ``pitch`` to within ``atol``.
        atol : float, optional
            Absolute tolerance, in the same units as ``pitch``, used to match the
            requested pitch angles to ``pitch``. Default is 1e-3.

        Returns
        -------
        numpy.ndarray
            Array of shape (N, M, len(groups)). Output bins where every averaged
            value is NaN are NaN.

        Raises
        ------
        ValueError
            If ``spectra`` is not 3-D, ``pitch`` does not match the length of axis 2 of
            ``spectra``, ``groups`` or one of its entries is empty, or a requested pitch
            angle is not found in ``pitch``.
        """
        spectra = np.asarray(spectra)
        pitch = np.asarray(self.pitch_angle)

        if spectra.ndim != 3:
            raise ValueError(f"spectra must be 3-D, got shape {spectra.shape}")
        if pitch.ndim != 1 or pitch.size != spectra.shape[2]:
            raise ValueError(f"pitch must be 1-D with length spectra.shape[2]={spectra.shape[2]}, "
                             f"got shape {pitch.shape}")
        if len(groups) == 0:
            raise ValueError("groups must contain at least one group of pitch angles")

        slices = []
        for group in groups:
            group = np.atleast_1d(group)
            if group.size == 0:
                raise ValueError("each group must contain at least one pitch angle")

            # --- find the index of each requested pitch angle ---
            idx = []
            for angle in group:
                i = np.argmin(np.abs(pitch - angle))
                if not np.isclose(pitch[i], angle, rtol=0, atol=atol):
                    raise ValueError(f"pitch angle {angle} not found (atol={atol}). "
                                     f"Available pitch angles: {pitch.tolist()}")
                idx.append(i)

            # --- single pitch angle: copy through; multiple: nan-average ---
            if len(idx) == 1:
                slices.append(spectra[:, :, idx[0]])
            else:
                with warnings.catch_warnings():
                    # all-NaN bins return NaN; suppress the "Mean of empty slice" warning
                    warnings.simplefilter('ignore', category=RuntimeWarning)
                    slices.append(np.nanmean(spectra[:, :, idx], axis=2))

        return np.stack(slices, axis=2)

    def time_average(self, spectra,n_avg, drop_remainder=False):
        """
        Nan-average consecutive time slices of an array along with their timestamps.

        Axis 0 of ``spectra`` (time) is divided into consecutive chunks of ``n_avg``
        slices, and each chunk is replaced by its NaN-ignoring mean. All other
        axes are unchanged, so a standardized (epoch, energy, pitch) array of
        shape (N, M, P) becomes (N', M, P), where N' is the number of chunks.
        Each chunk is assigned the midpoint between its first and last
        timestamps.

        Parameters
        ----------
        spectra : array_like
            Array with time along axis 0, e.g. the output of ``pitch_average``.
        epoch : array_like
            1-D array of the N timestamps of ``spectra``. May be numeric,
            ``numpy.datetime64``, or ``datetime.datetime`` objects.
        n_avg : int
            Number of consecutive time slices per chunk. 1 leaves the time
            resolution unchanged.
        drop_remainder : bool, optional
            How to handle trailing time slices when N is not divisible by
            ``n_avg``. If False (default), they form a shorter final chunk,
            giving N' = ceil(N / n_avg). If True, they are discarded, giving
            N' = N // n_avg.

        Returns
        -------
        A_spectra : numpy.ndarray
            Float64 array of shape ``(N', *spectra.shape[1:])``. Entries where every
            value in the chunk is NaN are NaN.
        epoch_avg : numpy.ndarray
            1-D array of the N' chunk midpoint times, aligned with axis 0 of
            ``A_spectra``. datetime64 input is returned as datetime64[ns],
            ``datetime.datetime`` objects as ``datetime.datetime`` objects, and
            numeric input as float.

        Raises
        ------
        ValueError
            If ``n_avg`` is not a positive integer, ``epoch`` is not 1-D with
            length ``spectra.shape[0]``, or ``drop_remainder`` is True and there are
            fewer than ``n_avg`` time slices.

        Notes
        -----
        Chunks are formed from consecutive indices, not fixed time windows, so a
        chunk that spans a data gap averages slices that are far apart in time.
        """
        spectra = np.asarray(spectra)
        epoch = np.asarray(self.epoch)

        if not (isinstance(n_avg, (int, np.integer)) and n_avg >= 1):
            raise ValueError(f"n_avg must be a positive integer, got {n_avg!r}")
        if epoch.ndim != 1 or epoch.size != spectra.shape[0]:
            raise ValueError(f"epoch must be 1-D with length spectra.shape[0]={spectra.shape[0]}, "
                             f"got shape {epoch.shape}")
        if np.issubdtype(epoch.dtype, np.datetime64):
            epoch = epoch.astype('datetime64[ns]')  # avoid truncating midpoints at coarse units

        n = spectra.shape[0]
        if drop_remainder:
            n = (n // n_avg) * n_avg
            if n == 0:
                raise ValueError(f"drop_remainder=True needs at least n_avg={n_avg} "
                                 f"time slices, got {spectra.shape[0]}")
            spectra, epoch = spectra[:n], epoch[:n]

        starts = np.arange(0, n, n_avg)
        ends = np.minimum(starts + n_avg, n)  # exclusive end index of each chunk

        # --- nan-mean of each chunk ---
        valid = ~np.isnan(spectra)
        sums = np.add.reduceat(np.where(valid, spectra, 0.0), starts, axis=0, dtype=np.float64)
        counts = np.add.reduceat(valid, starts, axis=0, dtype=np.int64)
        with np.errstate(invalid='ignore'):
            A_spectra = sums / counts  # 0 / 0 -> NaN for all-NaN chunks

        # --- midpoint time of each chunk ---
        first, last = epoch[starts], epoch[ends - 1]
        epoch_avg = first + (last - first) / 2

        return A_spectra, epoch_avg


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



