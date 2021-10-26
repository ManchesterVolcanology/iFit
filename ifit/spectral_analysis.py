"""Contains functions for fitting spectra."""
import logging
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import griddata
from scipy.signal import savgol_filter

from ifit.make_ils import make_ils


logger = logging.getLogger(__name__)


# =============================================================================
# =============================================================================
# # Spectral Analyser
# =============================================================================
# =============================================================================

class Analyser():
    """Object to perform the spectral analysis.

    Parameters
    ----------
    params : Parameters object
        The model parameters to include in the fit
    fit_window : tuple
        The lower and upper wavelength limits of the fit window (nm)
    frs_path : str
        Path to the Fraunhofer Reference Spectrum. This is read in using
        numpy.loadtxt and must consist of two columns [wavelength, intensity]
    model_padding : float, optional
        The padding in nm around the fit window to account for wavelength
        shifts and edge effects from the ils convolution. The default is
        1.0.
    model_spacing : float, optional
        The wavelength spacing of the model grid (nm). The default is 0.01.
    flat_flag : bool, optional
        If True then a flat correction is applied. The default is False.
    flat_path : str, optional
        Path to the flat spectrum used if flat_flag is True. The default is
        None.
    stray_flag : bool, optional
        If True then a stray correction is applied. The default is True.
    stray_window : tuple, optional
        The lower and upper wavelength limits of the stray light window.
        The default is [280, 290].
    dark_flag : bool, optional
        If True then a dark correction is applied. The default is True.
    ils_type : str, optional
        Controls how the ILS is generated. The default is "Manual"
            - 'File': Import a measured ils which is interpolated onto
                      the model grid
            - 'Params': Read in the ILS parameters from a file
            - 'Manual': Define the ILS parameters manually
    ils_path : str, optional
        Path to the ILS file, either a measured ILS is ils_type is set to
        "File", or the ILS parameters if ils_type is set to "Params".
        The default is None.
    despike_flag : bool, optional
        Controls whether bad pixels are removed. Requires spike_limit to be
        defined. The default is False
    spike_limit : int, optional
        The limit of the spike size to remove. Defined as the difference
        between the raw spectrum and a savgol filtered spectrum. The
        default is None

    Attributes
    ----------
    params : Parameters object
        A copy of the Parameters object supplied
    p0 : list
        The first guess fit parameters
    model_grid : numpy array
        The wavelength grid on which the calculations take place
    flat : numpy array
        The flat field spectrum
    flat_flag : bool
        Controls whether the flat spectrum is corrected in pre-processing
    ils : numpy array
        The ILS to use to smooth the model spectrum to the instrument
        resolution
    generate_ils : bool
        Controls whether the stored ILS is used (self.ils) or if the ILS is
        calculated at each model iteration (such as when being fit)
    frs : numpy array
        The Fraunhofer Reference Spectrum to be used by the forward model,
        interpolated onto the model grid
    xsecs : dict
        The cross-sections to use in the forward model, interpolated onto the
        model grid
    """

    def __init__(self, params, fit_window, frs_path, model_padding=1.0,
                 model_spacing=0.01, flat_flag=False, flat_path=None,
                 stray_flag=False, stray_window=[280, 290], dark_flag=False,
                 ils_type='Manual', ils_path=None, despike_flag=False,
                 spike_limit=None, bad_pixels=None):
        """Initialise the model for the analyser."""
        # Set the initial estimate for the fit parameters
        self.params = params.make_copy()
        self.p0 = self.params.fittedvalueslist()

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

            logger.info('Importing flat spectrum')

            try:
                # Import the flat spectrum
                self.flat = np.loadtxt(flat_path, unpack=True)
                logger.info('Flat spectrum imported')

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

        logger.info('Solar reference spectrum imported')

        # ---------------------------------------------------------------------
        # Import Gas spectra
        # ---------------------------------------------------------------------

        logger.info('Importing gas cross-sections...')

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
                self.init_xsecs[name] = griddata(x, xsec, self.model_grid,
                                                 method='cubic')

                logger.info(f'{name} cross-section imported')

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
        x, y = spectrum

        # Remove the dark spectrum from the measured spectrum
        if self.dark_flag:
            try:
                y = np.subtract(y, self.dark_spec)
            except ValueError:
                logger.exception('Error in dark correction. Is dark spectrum'
                                 + ' the same shape as the measurement?')

        # Remove stray light
        if self.stray_flag:
            stray_idx = np.where(np.logical_and(x >= self.stray_window[0],
                                                x <= self.stray_window[1]))

            if len(stray_idx[0]) == 0:
                logger.warn('No stray window outside spectrum, disabling '
                            + 'stray correction')
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
        fit_idx = np.where(np.logical_and(x >= self.fit_window[0],
                                          x <= self.fit_window[1]))
        grid = x[fit_idx]
        spec = y[fit_idx]

        # Divide by flat spectrum
        if self.flat_flag:

            # Unpack the flat spectrum and trim to the fit window
            flat_x, flat_y = self.flat
            flat_idx = np.where(np.logical_and(flat_x >= self.fit_window[0],
                                               flat_x <= self.fit_window[1]))
            flat = flat_y[flat_idx]

            # Divide the emasured spectrum by the flat spectrum
            try:
                spec = np.divide(spec, flat)
            except ValueError:
                logger.exception('Error in flat correction. Is flat spectrum'
                                 + ' the same shape as the measurement?')

        return np.row_stack([grid, spec])

# =============================================================================
#   fit_spectrum
# =============================================================================

    def fit_spectrum(self, spectrum, update_params=False, resid_limit=None,
                     resid_type='Percentage', int_limit=None, calc_od=[],
                     pre_process=True, prefit_shift=0.0, interp_method='cubic',
                     fit_window=None):
        """Fit the supplied spectrum.

        Parameters
        ----------
        spectrum : tuple
            The measured spectrum in the form [wavelengths, intensities]
        update_params : bool, optional
            Flag whether to update the initial fit parameters upon a sucessful
            fit. Can speed up fitting for large datasets but requires robust
            quality checks to avoid getting stuck
        resid_limit : float, optional, default=None
            Sets the maximum residual value to allow the fit initial parameters
            to be updated. Ignored if None.
        resid_type : str, optional, default='Percentage'
            Controls how the residual is calculated:
                - 'Percentage': calculated as (spec - fit) / spec * 100
                - 'Absolute':   calculated as spec - fit
        int_limit : tuple of floats, optional, default=None
            Sets the minimu and maximum intensity value for the fit window to
            allow the fit initial parameters to be updated. Ignored if None.
        calc_od : list, optional, default=[]
            List of parameters to calculate the optical depth for.
        pre_process : bool, optional, default=True
            Flag to control whether the supplied spectrum requires
            pre-prossessing or not
        prefit_shift : float, optional
            Wavelength shift (in nm) applied to the spectrum wavelength
            calibration prior to the fit. Default is 0.0
        interp_method : str, optional, default="cubic"
            Controls whether the interpolation at the end of the forward model
            is cubic or linear. Must be either "cubic", "linear" or "nearest".
            See scipy.interpolate.griddata for details.
        fit_window : tuple, optional
            Upper and lower limits of the fit window. This superceeds the main
            fit_window of the Analyser but must be contained within the window
            used to initialise the Analyser. Default is None.

        Returns
        -------
        fit_result : ifit.spectral_analysis.FitResult object
            An object that contains the fit results
        """
        # If a new fit window is given, trim the cross-sections down
        if fit_window is not None:

            # Check the new window is within the old one
            a = fit_window[0] < self.init_fit_window[0]
            b = fit_window[1] > self.init_fit_window[1]
            if a or b:
                logger.error('New fit window must be within initial fit'
                             + ' window!')
                raise ValueError

            # Pad the fit window
            pad_window = [fit_window[0] - self.model_padding,
                          fit_window[1] + self.model_padding]

            # Trim the model grid to the new fit window
            mod_idx = np.where(np.logical_and(self.init_grid >= pad_window[0],
                                              self.init_grid <= pad_window[1]))
            self.model_grid = self.init_grid[mod_idx]

            # Trim the FRS to the new fit window
            self.frs = self.init_frs[mod_idx]

            # Trim the gas cross-sections to the new fit window
            for key in self.init_xsecs.keys():
                self.xsecs[key] = self.init_xsecs[key][mod_idx]

            # Update the fit window attribute
            self.fit_window = fit_window

        # Check is spectrum requires preprocessing
        if pre_process:
            spectrum = self.pre_process(spectrum, prefit_shift)

        # Define the interpolation mode
        self.interp_method = interp_method

        # Unpack the spectrum
        grid, spec = spectrum

        # Fit the spectrum
        try:
            popt, pcov = curve_fit(self.fwd_model, grid, spec, self.p0)

            # Calculate the parameter error
            perr = np.sqrt(np.diag(pcov))

            # Set the success flag
            nerr = 1

        # If the fit fails return nans
        except RuntimeError:
            popt = np.full(len(self.p0), np.nan)
            perr = np.full(len(self.p0), np.nan)
            nerr = 0

        # Put the results into a FitResult object
        fit_result = FitResult(self, spectrum, popt, perr, nerr,
                               self.fwd_model, self.params, resid_type,
                               resid_limit, int_limit, calc_od)

        # If the fit was good then update the initial parameters
        if update_params and fit_result.nerr == 1:
            self.p0 = popt
        else:
            logger.info('Resetting initial guess parameters')
            self.p0 = self.params.fittedvalueslist()

        return fit_result

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
        plm_F = np.multiply(frs, plm_exponent)

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
        plume_F = np.multiply(plm_F, 1-ldf)

        # Build the baseline offset polynomial
        offset = np.polyval(offset_coefs, self.model_grid)

        # Combine the undiluted light, diluted light and offset
        raw_F = np.add(dilut_F, plume_F) + offset

        # Generate the ILS
        if self.generate_ils:

            # Unpack ILS params
            ils = make_ils(self.model_spacing,
                           p['fwem'],
                           p['k'],
                           p['a_w'],
                           p['a_k'])
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
# =============================================================================
# # Fit Result
# =============================================================================
# =============================================================================

class FitResult():
    """Contains the fit results.

    Parameters
    ----------
    analyser : Analyser
        Analyser object used to generate the results
    spectrum : tuple
        The fitted spectrum in the form [wavelengths, intensities]
    popt : list
        The optimised fit parameter values
    perr : list
        The calculated parameter fit errors
    nerr : int
        Signals success of the fit. 1 = successful, 0 = failed
    fwd_model : function
        The forward model used to perform the fit
    params : Parameters object
        The fit parameters
    resid_type : str
        Controls how the fit residual is calculated:
            - 'Percentage': calculated as (spec - fit) / spec * 100
            - 'Absolute':   calculated as spec - fit
    resid_limit : float
        Limit on the maximum residual for a good fit
    int_limit : float
        Limit on the spectrum intensity for a good fit
    calc_od : list of strings
        The parameters to calculate the optical depth spectrum for

    Attributes
    ----------
    params : Parameters object
        The Parameters object with the fitted values
    grid : numpy array
        The selected wavelength window
    spec : numpy array
        The selected intensity spectrum
    popt : numpy array
        The optimised fit results
    perr : numpy array
        The fit errors
    nerr : int
        The error code of the fit. 0 if the fit failed, 1 if the fit was
        successful and 2 if the fit was sucessful but failed a quality check
    fwd_model : function
        The forward model used to complete the fit
    meas_od : dict
        Dictionary of the measured optical depths
    synth_od : dict
        Dictionary of the synthetic optical depths
    int_lo, int_av, int_hi : float
        The min, average and max intensity in the fit window
    fit : numpy array
        The final fitted spectrum
    resid : numpy array
        The fit residual
    """

    def __init__(self, analyser, spectrum, popt, perr, nerr, fwd_model,
                 params, resid_type, resid_limit, int_limit, calc_od):
        """Initialize."""
        # Make a copy of the parameters
        self.params = params.make_copy()

        # Assign the variables
        self.grid, self.spec = spectrum
        self.popt = popt
        self.perr = perr
        self.nerr = nerr
        self.fwd_model = fwd_model
        self.meas_od = {}
        self.synth_od = {}

        # Calculate the spectrum intensity values
        self.int_lo = np.min(self.spec)
        self.int_hi = np.max(self.spec)
        self.int_av = np.average(self.spec)

        # Add the fit results to each parameter
        n = 0
        for par in self.params.values():

            if par.vary:
                par.set(fit_val=self.popt[n],
                        fit_err=self.perr[n])
                n += 1
            else:
                par.set(fit_val=par.value)
                par.set(fit_err=0)

        # If fit was successful then calculate the fit and residual
        if self.nerr:

            # Generate the fit
            self.fit = self.fwd_model(self.grid, *self.popt)

            # Calculate the residual
            if resid_type == 'Absolute':
                self.resid = self.spec - self.fit
            elif resid_type == 'Percentage':
                self.resid = (self.spec - self.fit)/self.spec * 100

            # Calculate residual values
            self.resid_max = np.nanmax(np.abs(self.resid))
            self.resid_std = np.nanstd(self.resid)

            # Check the fit quality
            if resid_limit is not None and max(abs(self.resid)) > resid_limit:
                logger.info('High residual detected')
                self.nerr = 2

            # Check for spectrum light levels
            if int_limit is not None:

                # Check for low intensity
                if min(self.spec) <= int_limit[0]:
                    logger.info('Low intensity detected')
                    self.nerr = 2

                # Check for high intensity
                elif max(self.spec) >= int_limit[1]:
                    logger.info('High intensity detected')
                    self.nerr = 2

            # Calculate optical depth spectra
            for par in calc_od:
                if par in self.params:
                    self.calc_od(par, analyser)

        # If not then return nans
        else:
            logger.warn('Fit failed!')
            self.fit = np.full(len(self.spec), np.nan)
            self.resid = np.full(len(self.spec), np.nan)
            for par in calc_od:
                if par in self.params:
                    self.meas_od[par] = np.full(len(self.spec), np.nan)
                    self.synth_od[par] = np.full(len(self.spec), np.nan)

# =============================================================================
#   Calculate Optical Depths
# =============================================================================

    def calc_od(self, par_name, analyser):
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

        fit = self.fwd_model(self.grid, *fit_params)

        # Calculate the shifted model grid
        shift_coefs = [p[n] for n in p if 'shift' in n]
        zero_grid = analyser.model_grid - min(analyser.model_grid)
        wl_shift = np.polyval(shift_coefs, zero_grid)
        shift_model_grid = np.add(analyser.model_grid, wl_shift)

        # Calculate the wavelength offset
        offset_coefs = [p[n] for n in p if 'offset' in n]
        offset = np.polyval(offset_coefs, analyser.model_grid)
        offset = griddata(shift_model_grid,
                          offset,
                          self.grid,
                          method='cubic')

        # Calculate the parameter od
        par_od = np.multiply(analyser.xsecs[par_name], p[par_name])

        # Make the ILS
        if analyser.generate_ils:

            ils_params = []
            for name in ['fwem', 'k', 'a_w', 'a_k']:
                if params[name].vary:
                    ils_params.append(params[name].fit_val)
                else:
                    ils_params.append(params[name].value)

            # Unpack ILS params
            ils = make_ils(analyser.model_spacing, *ils_params)
        else:
            ils = analyser.ils

        # Convolve with the ILS and interpolate onto the measurement grid
        par_od = griddata(shift_model_grid,
                          np.convolve(par_od, ils, mode='same'),
                          self.grid,
                          method='cubic')

        # Add to self
        self.meas_od[par_name] = -np.log(np.divide(self.spec-offset, fit))
        self.synth_od[par_name] = par_od

        return self.meas_od[par_name], self.synth_od[par_name]
