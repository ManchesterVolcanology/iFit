"""Contains functions for reading and fitting spectra."""
import os
import logging
import warnings
import numpy as np
import xarray as xr
import pandas as pd
from scipy.special import gamma
from scipy.optimize import curve_fit
from scipy.interpolate import griddata
from scipy.signal import savgol_filter


logger = logging.getLogger(__name__)

warnings.filterwarnings(
    "ignore",
    message="Covariance of the parameters could not be estimated"
)


# =============================================================================
# =============================================================================
# # Spectral Analyser
# =============================================================================
# =============================================================================

class Analyser(object):
    """
    Analyse UV spectra to retrieve trace gas information.
    """

    def __init__(
            self, params, fit_window, frs_path, model_padding=1.0,
            model_spacing=0.01, flat_flag=False, flat_path=None,
            stray_flag=False, stray_window=[280, 290], dark_flag=False,
            ils_type='Manual', ils_path=None, despike_flag=False,
            spike_limit=None, bad_pixels=None, update_params_flag=False,
            residual_limit=None, residual_type='Percentage',
            intensity_limit=None, pre_process_flag=True, interp_method='cubic',
            prefit_shift=0.0, spectrum_format=None, wavelength_calibration=None
        ):
        """Initialise the analyser."""
        self.ingester = Ingester(spectrum_format, wavelength_calibration)

        # Set the initial estimate for the fit parameters
        self.params = params.make_copy()
        self.p0 = self.params.fittedvalueslist()

        self.dark_spectrum = 0

        # ---------------------------------------------------------------------
        # Model Grid
        # ---------------------------------------------------------------------

        # Build model grid, a high res grid on which the forward model is build
        start = fit_window[0] - model_padding
        stop = fit_window[1] + model_padding + model_spacing
        self.init_grid = np.arange(start, stop, step=model_spacing)
        self.model_grid = self.init_grid.copy()

        # ---------------------------------------------------------------------
        # Flat Spectrum
        # ---------------------------------------------------------------------

        # Try importing flat spectrum
        if flat_flag:

            logger.info('Importing flat spectrum...')

            try:
                # Import the flat spectrum
                self.flat = np.loadtxt(flat_path, unpack=True)

            except OSError:
                # If no flat spectrum then report and turn off the flat flag
                logger.warning('No flat spectrum found!')
                self.flat_flag = False

        # ---------------------------------------------------------------------
        # Spectrometer ILS
        # ---------------------------------------------------------------------

        # Import measured ILS
        if ils_type == 'File':
            logger.info('Importing ILS')

            try:
                # Read in measured ILS shape
                x_ils, y_ils = np.loadtxt(ils_path, unpack=True)

                # Interpolate the measured ILS onto the model grid spacing
                grid_ils = np.arange(x_ils[0], x_ils[-1], model_spacing)
                ils = griddata(x_ils, y_ils, grid_ils, 'cubic')
                self.ils = ils / np.sum(ils)
                self.generate_ils = False

            except OSError:
                logger.error(f'{ils_path} file not found!')

        # Import ILS params
        if ils_type == 'Params':
            logger.info('Importing ILS parameters')
            try:
                # Import ils parameters
                ils_params = np.loadtxt(ils_path, unpack=True)

                # Build ILS
                self.ils = make_ils(model_spacing, *ils_params)

                self.generate_ils = False
                logger.info(f'ILS parameters imported: {ils_path}')

            except OSError:
                logger.error(f'{ils_path} file not found!')

        # Manually set the ILS params
        if ils_type == 'Manual':

            # Check if they're all fixed
            keys = ['fwem', 'k', 'a_w', 'a_k']
            vary_check = np.array([params[k].vary for k in keys])
            if vary_check.any():
                self.generate_ils = True
            else:
                ils_params = [params[k].value for k in keys]
                self.ils = make_ils(model_spacing, *ils_params)
                self.generate_ils = False

        # ---------------------------------------------------------------------
        # Import solar spectrum
        # ---------------------------------------------------------------------

        # Import solar reference spectrum
        logger.info('Importing solar reference spectrum...')
        sol_x, sol_y = np.loadtxt(frs_path, unpack=True)

        # Interpolate onto model_grid
        self.init_frs = griddata(sol_x, sol_y, self.model_grid, method='cubic')
        self.frs = self.init_frs.copy()

        # ---------------------------------------------------------------------
        # Import Gas spectra
        # ---------------------------------------------------------------------

        # Create an empty dictionary to hold the gas cross-sections
        self.init_xsecs = {}

        # Cycle through the parameters
        for name, param in self.params.items():

            # If a parameter has a xpath defined, read it in
            if param.xpath is not None:
                logger.info(f'Importing {name} reference spectrum...')

                # Read in the cross-section
                x, xsec = np.loadtxt(param.xpath, unpack=True)

                # Interpolate onto the model grid
                self.init_xsecs[name] = griddata(
                    x, xsec, self.model_grid, method='cubic'
                )

        logger.info('Analyser setup complete')

        # Create a copy of the cross-sections
        self.xsecs = self.init_xsecs.copy()

        # ---------------------------------------------------------------------
        # Other model settings
        # ---------------------------------------------------------------------

        self.init_fit_window = fit_window
        self.fit_window = fit_window
        self.stray_window = stray_window
        self.stray_flag = stray_flag
        self.flat_flag = flat_flag
        self.dark_flag = dark_flag
        self.model_padding = model_padding
        self.model_spacing = model_spacing
        self.despike_flag = despike_flag
        self.spike_limit = spike_limit
        self.bad_pixels = bad_pixels
        self.update_params_flag = update_params_flag
        self.residual_limit = residual_limit
        self.residual_type = residual_type
        self.intensity_limit = intensity_limit
        self.pre_process_flag = pre_process_flag
        self.interp_method = interp_method
        self.prefit_shift = prefit_shift

    def fit_spectrum(self, spectrum, calc_od=[], fit_window=None):
        """Fit the supplied spectrum.

        Parameters
        ----------
        spectrum : str or xarray DataArray
            The spectrum to analyse. If a string, then this is the path to the
            spectrum file. Otherwise a DataArray can be passed directly. Note
            that the coords must be labelled "wavelength"
        calc_od : list, optional, default=[]
            List of parameters to calculate the optical depth for.
        fit_window : tuple, optional
            Upper and lower limits of the fit window. This superceeds the main
            fit_window of the Analyser but must be contained within the window
            used to initialise the Analyser. Default is None.

        Returns
        -------
        fit_result : xarray Dataset
            Contains the fit results
        """
        # If a new fit window is given, trim the cross-sections down
        if fit_window is not None:

            # Check the new window is within the old one
            a = fit_window[0] < self.init_fit_window[0]
            b = fit_window[1] > self.init_fit_window[1]
            if a or b:
                logger.error(
                    'New fit window must be within initial fit window!'
                )
                raise ValueError

            # Pad the fit window
            pad_window = [
                fit_window[0] - self.model_padding,
                fit_window[1] + self.model_padding
            ]

            # Trim the model grid to the new fit window
            mod_idx = np.where(np.logical_and(
                self.init_grid >= pad_window[0],
                self.init_grid <= pad_window[1]
            ))
            self.model_grid = self.init_grid[mod_idx]

            # Trim the FRS to the new fit window
            self.frs = self.init_frs[mod_idx]

            # Trim the gas cross-sections to the new fit window
            for key in self.init_xsecs.keys():
                self.xsecs[key] = self.init_xsecs[key][mod_idx]

            # Update the fit window attribute
            self.fit_window = fit_window

        # If spectrum is a path, then ingest the file
        if isinstance(spectrum, str):
            spectrum = self.ingester.ingest(spectrum)

        # Check is spectrum requires preprocessing
        if self.pre_process_flag:
            spectrum = self.pre_process(spectrum, self.prefit_shift)

        # Initialise holder for fit data
        coords = spectrum.coords
        fit_data = {
            'spectrum': xr.DataArray(
                data=spectrum.data,
                coords=coords,
                attrs={'units': 'counts'}
            )
        }

        # Set residual units
        if self.residual_type == 'Absolute':
            resid_units = 'counts'
        elif self.residual_type == 'Percentage':
            resid_units = '%'

        # Fit the spectrum
        try:
            popt, pcov = curve_fit(
                self.fwd_model,
                spectrum.wavelength,
                spectrum.data,
                self.p0
            )

            # Calculate the parameter error
            perr = np.sqrt(np.diag(pcov))

             # Add the fit results to each parameter
            n = 0
            for par in self.params.values():

                if par.vary:
                    par.set(fit_val=popt[n], fit_err=perr[n])
                    n += 1
                else:
                    par.set(fit_val=par.value)
                    par.set(fit_err=0)

            # Set the success flag
            fit_quality = 1

            # Generate the fit
            fit = self.fwd_model(spectrum.wavelength, *popt)

            # Calculate the residual
            if self.residual_type == 'Absolute':
                resid = spectrum.data - fit
            elif self.residual_type == 'Percentage':
                resid = (spectrum.data - fit)/spectrum.data * 100

            # Check the fit quality
            if self.residual_limit is not None \
                    and max(abs(resid)) > self.residual_limit:
                logger.info('High residual detected')
                fit_quality = 2

            # Check for spectrum light levels
            if self.intensity_limit is not None:

                # Check for low intensity
                if min(spectrum.data) <= self.intensity_limit[0]:
                    logger.info('Low intensity detected')
                    fit_quality = 2

                # Check for high intensity
                elif max(spectrum.data) >= self.intensity_limit[1]:
                    logger.info('High intensity detected')
                    fit_quality = 2

            # Collate the fit results
            fit_data['fit'] = xr.DataArray(
                data=fit, coords=coords, attrs={'units': 'counts'}
            )
            fit_data['residual'] = xr.DataArray(
                data=resid, coords=coords, attrs={'units': resid_units}
            )

            # Calculate optical depth spectra
            for par_name in calc_od:
                if par_name in self.params:
                    meas_od, synth_od = self.calc_od(spectrum, par_name)
                    fit_data[f'{par_name}_meas_od'] = xr.DataArray(
                        data=meas_od, coords=coords, attrs={'units': 'arb'}
                    )
                    fit_data[f'{par_name}_synth_od'] = xr.DataArray(
                        data=synth_od, coords=coords, attrs={'units': 'arb'}
                    )

        # If the fit fails return nans
        except RuntimeError:
            logger.warn('Fit failed!')

            # Set fit quality to 0
            fit_quality = 0

            # Set the ouput parameter values
            n = 0
            for par in self.params.values():
                if par.vary:
                    par.set(fit_val=np.nan, fit_err=np.nan)
                    n += 1
                else:
                    par.set(fit_val=par.value)
                    par.set(fit_err=0)

            # Set output spectra to nans
            nan_arr = np.full(len(spectrum), np.nan)
            fit_data['fit'] = xr.DataArray(
                data=nan_arr, coords=coords, attrs={'units': 'counts'}
            )
            fit_data['residual'] = xr.DataArray(
                data=nan_arr, coords=coords, attrs={'units': resid_units}
            )

            # Calculate optical depth spectra
            for par_name in calc_od:
                if par_name in self.params:
                    fit_data[f'{par_name}_meas_od'] = xr.DataArray(
                        data=nan_arr, coords=coords, attrs={'units': 'arb'}
                    )
                    fit_data[f'{par_name}_synth_od'] = xr.DataArray(
                        data=nan_arr, coords=coords, attrs={'units': 'arb'}
                    )

        fit_result = xr.Dataset(
            data_vars=fit_data,
            attrs={'fit_quality': fit_quality, **spectrum.attrs}
        )

        # If the fit was good then update the initial parameters
        if self.update_params_flag:
            if fit_result.fit_quality == 1:
                self.p0 = popt
            else:
                logger.info('Resetting initial guess parameters')
                self.p0 = self.params.fittedvalueslist()

        return fit_result

# =============================================================================
#   Set dark spectrum
# =============================================================================

    def set_dark_spectrum(self, spectrum):

        # If spectrum is a path, then ingest the file
        if isinstance(spectrum, str):
            spectrum = self.ingester.ingest(spectrum)

        # Assign
        self.dark_spectrum = spectrum

# =============================================================================
#   Spectrum Pre-processing
# =============================================================================

    def pre_process(self, spectrum, prefit_shift=0.0):
        """Prepare spectrum for fit.

        Function to pre-process the measured spectrum to prepare it for the
        fit, correcting for the dark and flat spectrum, stray light, spiky
        pixels and extracting the fit wavelength window. Which corrections are
        applied depends on the settings of the Analyser object.

        Parameters
        ----------
        spectrum : 2D numpy array
            The spectrum as [wavelength, intensities].
        prefit_shift : float, optional
            Wavelength shift (in nm) applied to the spectrum wavelength
            calibration prior to the fit. Default is 0.0

        Returns
        -------
        processed_spec : 2D numpy array
            The processed spectrum.
        """
        # Unpack spectrum
        spectrum = spectrum.copy()
        x = spectrum.wavelength
        y = spectrum.data

        # Remove the dark spectrum from the measured spectrum
        if self.dark_flag:
            try:
                y = np.subtract(y, self.dark_spectrum)
            except ValueError:
                logger.exception(
                    'Error in dark correction. Is dark spectrum the same shape'
                    ' as the measurement?'
                )

        # Remove stray light
        if self.stray_flag:
            stray_idx = np.where(np.logical_and(
                x >= self.stray_window[0],
                x <= self.stray_window[1]
            ))

            if len(stray_idx[0]) == 0:
                logger.warn(
                    'No stray window outside spectrum, disabling '
                    'stray correction'
                )
                self.stray_flag = False

            else:
                y = np.subtract(y, np.average(y[stray_idx]))

        # Run de-spike
        if self.despike_flag:

            # Run a savgol filter on the spectrum
            sy = savgol_filter(y, 11, 3)

            # Calculate the difference
            dspec = np.abs(np.subtract(y, sy))

            # Find any points that are over the spike limit and replace with
            # smoothed values
            spike_idx = np.where(dspec > self.spike_limit)[0]
            for i in spike_idx:
                y[i] = sy[i]

        # Remove bad pixels
        if self.bad_pixels is not None:
            for i in self.bad_pixels:
                y[i] = np.average([y[i-1], y[i+1]])

        # Apply prefit shift
        x = np.add(x, prefit_shift)

        # Cut desired wavelength window
        fit_idx = np.where(np.logical_and(
            x >= self.fit_window[0],
            x <= self.fit_window[1]
        ))
        grid = x[fit_idx]
        spec = y[fit_idx]

        # Divide by flat spectrum
        if self.flat_flag:

            # Unpack the flat spectrum and trim to the fit window
            flat_x, flat_y = self.flat
            flat_idx = np.where(np.logical_and(
                flat_x >= self.fit_window[0],
                flat_x <= self.fit_window[1]
            ))
            flat = flat_y[flat_idx]

            # Divide the emasured spectrum by the flat spectrum
            try:
                spec = np.divide(spec, flat)
            except ValueError:
                logger.exception(
                    'Error in flat correction. Is flat spectrum'
                    ' the same shape as the measurement?'
                )

        out_spectrum = xr.DataArray(
            data = spec,
            coords={'wavelength': grid},
            attrs={
                'fit_window_lower_limit': self.fit_window[0],
                'fit_window_upper_limit': self.fit_window[1],
                **spectrum.attrs
            }
        )

        return out_spectrum

# =============================================================================
#   Forward Model
# =============================================================================

    def fwd_model(self, x, *p0):
        """Forward model for iFit to fit measured UV sky spectra.

        Parameters
        ----------
        x, array
            Measurement wavelength grid
        *p0, floats
            Forward model state vector. Should consist of:
                - bg_poly{n}: Background polynomial coefficients
                - offset{n}:  The intensity offset polynomial coefficients
                - shift{n}:   The wavelength shift polynomial
                - gases:      Any variable with an associated cross section,
                              including absorbing gases and Ring. Each "gas" is
                              converted to transmittance through:
                              gas_T = exp(-xsec . amt)

            For polynomial parameters n represents ascending intergers
            starting from 0 which correspond to the decreasing power of
            that coefficient

        Returns
        -------
        fit, array
            Fitted spectrum interpolated onto the spectrometer wavelength grid
        """
        # Get dictionary of fitted parameters
        params = self.params
        p = params.valuesdict()

        # Update the fitted parameter values with those supplied to the forward
        # model
        i = 0
        for par in params.values():
            if par.vary:
                p[par.name] = p0[i]
                i += 1
            else:
                p[par.name] = par.value

        # Unpack polynomial parameters
        bg_poly_coefs = [p[n] for n in p if 'bg_poly' in n]
        offset_coefs = [p[n] for n in p if 'offset' in n]
        shift_coefs = [p[n] for n in p if 'shift' in n]

        # Construct background polynomial
        bg_poly = np.polyval(bg_poly_coefs, self.model_grid)
        frs = np.multiply(self.frs, bg_poly)

        # Create empty arrays to hold optical depth spectra
        plm_gas_T = np.zeros((len(self.xsecs), len(self.model_grid)))
        sky_gas_T = np.zeros((len(self.xsecs), len(self.model_grid)))

        # Calculate the gas optical depth spectra
        for n, gas in enumerate(self.xsecs):
            if self.params[gas].plume_gas:
                plm_gas_T[n] = (np.multiply(self.xsecs[gas], p[gas]))
            else:
                sky_gas_T[n] = (np.multiply(self.xsecs[gas], p[gas]))

        # Sum the gas ODs
        sum_plm_T = np.sum(np.vstack([plm_gas_T, sky_gas_T]), axis=0)
        sky_plm_T = np.sum(sky_gas_T, axis=0)

        # Build the exponent term
        plm_exponent = np.exp(-sum_plm_T)
        sky_exponent = np.exp(-sky_plm_T)

        # Build the complete model
        sky_F = np.multiply(frs, sky_exponent)
        plume_F = np.multiply(frs, plm_exponent)

        # Add effects of light dilution
        if 'LDF' in p and p['LDF'] != 0:

            # Calculate constant light dilution
            ldf_const = - np.log(1-p['LDF'])*(310**4)

            # Add wavelength dependancy to light dilution factor
            rayleigh_scale = self.model_grid**-4
            ldf = 1-np.exp(-ldf_const * rayleigh_scale)

        else:
            ldf = 0

        # Construct the plume and diluting light spectra, scaling by the ldf
        dilut_F = np.multiply(sky_F, ldf)

        # Build the baseline offset polynomial
        offset = np.polyval(offset_coefs, self.model_grid)

        # Combine the undiluted light, diluted light and offset
        raw_F = np.add(dilut_F, plume_F) + offset

        # Generate the ILS
        if self.generate_ils:

            # Unpack ILS params
            ils = make_ils(
                self.model_spacing,
                p['fwem'], p['k'], p['a_w'], p['a_k']
            )
        else:
            ils = self.ils

        # Apply the ILS convolution
        F_conv = np.convolve(raw_F, ils, 'same')

        # Apply shift and stretch to the model_grid
        zero_grid = self.model_grid - min(self.model_grid)
        wl_shift = np.polyval(shift_coefs, zero_grid)
        shift_model_grid = np.add(self.model_grid, wl_shift)

        # Interpolate onto measurement wavelength grid
        fit = griddata(shift_model_grid, F_conv, x, method=self.interp_method)

        return fit

# =============================================================================
#   Calculate Optical Depths
# =============================================================================

    def calc_od(self, spectrum, par_name):
        """Calculate the optical depth for the given parameter.

        Parameters
        ----------
        par_name : str
            The key of the parameter to calculate optical depth for
        analyser : Analyser object
            The Analyser used to create the results

        Returns
        -------
        mead_od : numpy array
            The measured optical depth, calculated by removing the fitted gas
            from the measured spectrum
        synth_od : numpy array
            The synthetic optical depth, calculated by multiplying the
            parameter cross-section by the fitted amount
        """
        # Make a copy of the parameters to use in the OD calculation
        params = self.params.make_copy()

        # Set the parameter and any offset coefficients to zero
        params[par_name].set(fit_val=0)
        for par in params:
            if 'offset' in par:
                params[par].set(fit_val=0)

        # Calculate the fit without the parameter
        fit_params = params.popt_list()
        p = self.params.popt_dict()

        fit = self.fwd_model(spectrum.wavelength, *fit_params)

        # Calculate the shifted model grid
        shift_coefs = [p[n] for n in p if 'shift' in n]
        zero_grid = self.model_grid - min(self.model_grid)
        wl_shift = np.polyval(shift_coefs, zero_grid)
        shift_model_grid = np.add(self.model_grid, wl_shift)

        # Calculate the wavelength offset
        offset_coefs = [p[n] for n in p if 'offset' in n]
        offset = np.polyval(offset_coefs, self.model_grid)
        offset = griddata(
            shift_model_grid,
            offset,
            spectrum.wavelength,
            method='cubic'
        )

        # Calculate the parameter od
        par_od = np.multiply(self.xsecs[par_name], p[par_name])

        # Make the ILS
        if self.generate_ils:

            ils_params = []
            for name in ['fwem', 'k', 'a_w', 'a_k']:
                if params[name].vary:
                    ils_params.append(params[name].fit_val)
                else:
                    ils_params.append(params[name].value)

            # Unpack ILS params
            ils = make_ils(self.model_spacing, *ils_params)
        else:
            ils = self.ils

        # Convolve with the ILS and interpolate onto the measurement grid
        par_od = griddata(
            shift_model_grid,
            np.convolve(par_od, ils, mode='same'),
            spectrum.wavelength,
            method='cubic'
        )

        # Add to self
        meas_od = -np.log(np.divide(spectrum-offset, fit))
        synth_od = par_od

        return meas_od, synth_od


# =============================================================================
# =============================================================================
# # Spectrum Ingester
# =============================================================================
# =============================================================================

class Ingester(object):
    """
    Ingest UV spectra for analysis.
    """

    def __init__(self, spectrum_format=None, wavelength_calibration=None):
        """Initialise."""
        # Create a dictionary of available functions
        self.read_functions = {
            'netCDF': self._load_netcdf,
            'iFit': self._load_ifit,
            'iFit (old)': self._load_ifit_old,
            'Master.Scope': self._load_master_scope,
            'Spectrasuite': self._load_spectrasuite,
            'mobileDOAS': self._load_mobile_doas,
            'Basic': self._load_basic
        }

        self.spectrum_format = spectrum_format
        self.wavelength_calibration = wavelength_calibration

    def ingest(self, fname):
        """."""

        # Check if the spectrum format has been defined
        if self.spectrum_format is not None:
            da = self.read_functions[self.spectrum_format](fname)
            return da

        else:
            # Cycle through the formats until one works
            logger.debug('Auto-detecting spectrum format')
            for key, func in self.read_functions.items():
                try:
                    da = func(fname)
                    self.spectrum_format = key
                    logger.debug(f'{key} format found!')
                    return da
                except Exception:
                    continue

            # If nothing works, raise an exception
            raise ValueError(f'Unable to load file {fname}')

    def _load_netcdf(self, fname):
        """Load netCDF spectrum."""
        return xr.load_dataarray(fname)

    def _load_ifit(self, fname):
        """Load iFit text spectrum format."""
        # Load data into a numpy array
        grid, spec = np.loadtxt(fname, unpack=True)

        # Assert metadata data types
        int_types = ['spectrum_number', 'integration_time', 'coadds']
        float_types = ['lat', 'lon', 'alt']
        bool_types = ['elecdk_correction', 'nonlin_correction']
        time_types = ['timestamp']

        # Extract metadata from the header
        with open(fname, 'r') as r:

            # Initialise the metadata dictionary
            metadata = {}

            # Read the first line, this isn't used
            line = r.readline()

            # Cycle through until we hit the data header
            line = r.readline()
            while 'Wavelength' not in line:

                # Split the key from the value
                key, value = line[2:].strip().split('; ')

                # For each, record according to the perscribed datatype
                if key in int_types and value != 'None':
                    metadata[key] = int(value)
                elif key in float_types and value != 'None':
                    metadata[key] = float(value)
                elif key in bool_types and value != 'None':
                    if value == 'True':
                        value = True
                    elif value == 'False':
                        value = False
                    metadata[key] = value
                elif key in time_types and value != 'None':
                    ts = pd.Timestamp(value)
                    metadata[key] = ts
                else:
                    metadata[key] = value

                # Read the next line
                line = r.readline()

        # Colate into an output DataArray
        da = xr.DataArray(
            data = spec,
            coords={'wavelength': grid},
            attrs={'units': 'counts', 'wavelength units': 'nm', **metadata}
        )
        return da

    def _load_ifit_old(self, fname):
        """Read older iFit format files."""
        # Load data into a numpy array
        grid, spec = np.loadtxt(fname, unpack=True)

        # Extract metadata from the header
        with open(fname, 'r') as r:
            lines = r.readlines()
            serial_number = lines[1].strip().split(': ')[-1]
            try:
                integration_time = int(lines[2].strip().split(': ')[-1][:-4])
            except ValueError:
                integration_time = int(lines[2].strip().split(': ')[-1])
            coadds = int(lines[3].strip().split(': ')[-1])
            time_string = lines[4].strip().split(': ')[-1]
            timestamp = pd.Timestamp(time_string)
            elecdk_correction = bool(lines[5].strip().split(': ')[-1])
            nonlin_correction = bool(lines[6].strip().split(': ')[-1])

        # Get spectrum number
        try:
            spec_no = int(fname[-9:-4])
        except ValueError:
            spec_no = 0

        metadata = {
            'spectrum_number': spec_no,
            'serial_number': serial_number,
            'integration_time': integration_time,
            'coadds': coadds,
            'timestamp': timestamp,
            'elecdk_correction': elecdk_correction,
            'nonlin_correction': nonlin_correction
        }

        # Colate into an output DataArray
        da = xr.DataArray(
            data = spec,
            coords={'wavelength': grid},
            attrs={'units': 'counts', 'wavelength units': 'nm', **metadata}
        )
        return da

    def _load_master_scope(self, fname):
        """Load Master.Scope file."""
        # Load data into a numpy array, skipping the header data
        grid, spec = np.genfromtxt(
            fname, unpack=True, skip_header=19, skip_footer=1
        )

        # Extract metadata from the header
        with open(fname, 'r') as r:
            lines = r.readlines()
            serial_number = lines[4].strip().split(': ')[-1]
            integration_time = int(lines[6].strip().split(': ')[-1])
            coadds = int(lines[7].strip().split(': ')[-1])
            time_string = lines[2].strip().split(': ')[-1]
            timestamp = pd.Timestamp(time_string)
            elecdk_correction = lines[9].strip().split(': ')[-1]
            if elecdk_correction == 'Disabled':
                elecdk_correction = False
            else:
                elecdk_correction = True

        # Get spectrum number
        try:
            spec_no = int(fname[-18:-13])
        except ValueError:
            spec_no = 0

        metadata = {
            'spectrum_number': spec_no,
            'serial_number': serial_number,
            'integration_time': integration_time,
            'coadds': coadds,
            'timestamp': timestamp,
            'elecdk_correction': elecdk_correction
        }

        # Colate into an output DataArray
        da = xr.DataArray(
            data = spec,
            coords={'wavelength': grid},
            attrs={'units': 'counts', 'wavelength units': 'nm', **metadata}
        )
        return da

    def _load_spectrasuite(self, fname):
        """Load Spectrasuite file."""
        # Load data into a numpy array, skipping the header data
        grid, spec = np.genfromtxt(
            fname, unpack=True, skip_header=17, skip_footer=2
        )

        # Extract metadata from the header
        with open(fname, 'r') as r:
            lines = r.readlines()
            serial_number = lines[7].strip().split(': ')[-1]
            integration_time = int(lines[6].strip().split(': ')[-1])
            coadds = int(lines[7].strip().split(': ')[-1])
            time_string = lines[2].strip().split(': ')[-1]
            timestamp = pd.Timestamp(time_string)
            elecdk_correction = lines[9].strip().split(': ')[-1]
            if elecdk_correction == 'Disabled':
                elecdk_correction = False
            else:
                elecdk_correction = True

        # Get spectrum number
        try:
            spec_no = int(fname[-9:-4])
        except ValueError:
            spec_no = 0

        metadata = {
            'spectrum_number': spec_no,
            'serial_number': serial_number,
            'integration_time': integration_time,
            'coadds': coadds,
            'timestamp': timestamp,
            'elecdk_correction': elecdk_correction
        }

        # Colate into an output DataArray
        da = xr.DataArray(
            data = spec,
            coords={'wavelength': grid},
            attrs={'units': 'counts', 'wavelength units': 'nm', **metadata}
        )
        return da

    def _load_mobile_doas(self, fname):
        """Load mobileDOAS file."""

        # Ensure the waveength calibration has been defined
        if self.wavelength_calibration is None:
            raise ValueError('Wavelength calibration not defined!')

        # Get the spectrum number
        _, tail = os.path.split(fname)
        try:
            spec_no = int(tail[:5])
        except ValueError:
            spec_no = 0

        # Format created by the mobileDOAS program
        with open(fname, 'r') as r:

            # Read in the file
            lines = [line.strip() for line in r.readlines()]

            # Get the number of pixels
            npixels = int(lines[2])

            # Pull out the spectral information
            spec = np.array([float(y) for y in lines[3:3+npixels]])

            # Pull out the metadata
            mlines = lines[npixels+3:]

            # Get the timestamp
            timestamp = pd.Timestamp(f'{mlines[3]} {mlines[4]}')

            metadata = {
                'spectrum_number': spec_no,
                'filename': mlines[0],
                'serial_number': mlines[1],
                'timestamp': timestamp,
                'unknown1': float(mlines[6]),
                'unknown2': float(mlines[7]),
                'coadds': int(mlines[8].split(' ')[-1]),
                'integration_time': int(mlines[9].split(' ')[-1]),
                'site': mlines[10].split(' ')[-1],
                'lon': float(mlines[11].split(' ')[-1]),
                'lat': float(mlines[12].split(' ')[-1]),
                'alt': float(mlines[12].split(' ')[-1])
            }

        # Colate into an output DataArray
        da = xr.DataArray(
            data = spec,
            coords={'wavelength': np.loadtxt(self.wavelength_calibration)},
            attrs={'units': 'counts', 'wavelength units': 'nm', **metadata}
        )
        return da

    def _load_basic(self, fname):
        """Load basic text column file."""
        # Load the file
        data = np.loadtxt(fname)

        # Check the dimensions
        if data.ndim == 1:
            # For one dimension, look for the wavelength calibration
            if self.wavelength_calibration is None:
                raise ValueError('Wavelength calibration not defined!')
            grid = np.loadtxt(self.wavelength_calibration)
            spec = data
        else:
            grid = data[:, 0]
            spec = data[:, 1]

        # Colate into an output DataArray
        da = xr.DataArray(
            data = spec,
            coords={'wavelength': grid},
            attrs={'units': 'counts', 'wavelength units': 'nm'}
        )
        return da


# =============================================================================
# =============================================================================
# # Spectrometer ILS
# =============================================================================
# =============================================================================

# =============================================================================
# Super Gaussian function
# =============================================================================

def super_gaussian(grid, w, k, a_w, a_k, shift=0, amp=1, offset=0):
    """Return a super-Gaussian line shape."""
    # Compute A
    A = k / (2 * w * gamma(1/k))

    # Form empty array
    ils = np.zeros(len(grid))

    # Iterate over x grid. If negative do one thing, if positive do the other
    ils = np.array([
        left_func(x, w, k, a_w, a_k) if x <= 0
        else right_func(x, w, k, a_w, a_k)
        for x in grid
    ])

    # Shift the lineshape
    if shift != 0:
        mod_grid = grid + shift

        # Interpolate onto the measurement grid
        ils = griddata(mod_grid, ils, grid, method='cubic', fill_value=0.0)

    return ils * A * amp + offset


def left_func(x, w, k, a_w, a_k):
    """Left function for asymetric Gaussian."""
    return np.exp(-np.power(np.abs((x) / (w - a_w)), k - a_k))


def right_func(x, w, k, a_w, a_k):
    """Right function for asymetric Gaussian."""
    return np.exp(-np.power(np.abs((x) / (w + a_w)), k + a_k))


# =============================================================================
# make_ils
# =============================================================================

def make_ils(interval, FWEM, k=2, a_w=0, a_k=0):
    """Generate a synthetic instrument line shape.

    Generates a lineshape based on the super-Gaussian function:

    .                { exp(-| x / (w-a_w) | ^ (k-a_k)) for x <= 0
    G(x) = A(w, k) * {
    .                { exp(-| x / (w+a_w) | ^ (k+a_k)) for x > 0

    where A(w, k) = k / (2 * w * Gamma(1/k)).

    See Beirle et al (2017) for more details: doi:10.5194/amt-10-581-2017

    Parameters
    ----------
    interval : int
        The spacing of the wavelength grid on which the ILS is built
    FWEM : float
        The Full Width eth Maximum of the lineshape, defined as FWEM = 2*w
    k : float, optional
        Controls the shape of the lineshape (default = 2):
            - k < 2 -> sharp point and wide tails
            - k = 2 -> normal Gaussian
            - k > 2 -> flat top, approaches boxcar at k -> inf
    a_w and a_k : float, optional
        Controls the asymetry of the lineshape. Defaults are 0

    Returns
    -------
    ils : numpy array
        The calculated ILS function on a wavelength grid of the given spacing
        and 5 times the width of the supplied FWEM
    """
    # Create a 4 nm grid
    grid = np.arange(-2, 2, interval)

    # Calculate w as half of the FWEM
    w = 0.5 * FWEM

    # Make the line shape
    ils = super_gaussian(grid, w, k, a_w, a_k)

    ils = np.divide(ils, sum(ils))

    return ils
