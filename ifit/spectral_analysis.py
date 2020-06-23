#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 11:28:07 2020
@author: mqbpwbe2
"""

import logging
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import griddata

from ifit.make_ils import make_ils

# =============================================================================
# =============================================================================
# # Spectral Analyser
# =============================================================================
# =============================================================================

class Analyser():

    '''
    Object to perform the spectral analysis
    Parameters
    ----------
    common : dict
        Common dictionary of parameters and variables passed from the main
        program to subroutines
    '''

    def __init__(self, params, fit_window, frs_path, model_padding=1.0, 
                 model_spacing=0.01, flat_flag=False, flat_path=None,
                 stray_flag=False, stray_window=[280, 290], dark_flag=False,
                 ils_type='Manual', ils_path=None):
        '''
        Initialise the model for the analyser

        Parameters
        ----------
        params : Parameters object
            The model parameters to include in the fit.
        fit_window : tuple
            The lower and upper wavelength limits of the fit window (nm).
        frs_path : str
            Path to the Fraunhofer Reference Spectrum. This is read in using
            numpy.loadtxt and must consist of two columns: 
            [wavelength, intensity]
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
            Controls how the ILS is generated:
                - "File"   : Import a measured ils which is interpolated onto 
                             the model grid
                - "Params" : Read in the ILS parameters from a file
                - "Manual" : Define the ILS parameters manually.
            The default is 'Manual'.
        ils_path : str, optional
            Path to the ILS file, either a measured ILS is ils_type is set to
            "File", or the ILS parameters if ils_type is set to "Params". 
            The default is None.

        Returns
        -------
        None.

        '''
        
        # Set the initial estimate for the fit parameters
        self.params = params.make_copy()
        self.p0 = self.params.fittedvalueslist()
        
        # ---------------------------------------------------------------------
        # Model Grid
        # ---------------------------------------------------------------------
        
        # Build model grid, a high res grid on which the forward model is build
        start = fit_window[0] - model_padding
        stop = fit_window[1] + model_padding + model_spacing
        model_grid = np.arange(start, stop, step = model_spacing)
    
        # Add the model grid to the common dict
        self.model_grid = model_grid
        
        # ---------------------------------------------------------------------
        # Flat Spectrum
        # ---------------------------------------------------------------------
        
        # Try importing flat spectrum
        if flat_flag == True:
    
            logging.info('Importing flat spectrum')
    
            try:
                # Import flat spectrum and extract the fit window
                flat_x, flat_y = np.loadtxt(flat_path, unpack=True)
                f_idx = np.where(np.logical_and(flat_x >= fit_window[0],
                                                flat_x <= fit_window[1]))
                self.flat = flat_y[f_idx]
    
                logging.info('Flat spectrum imported')
    
            except OSError:
                # If no flat spectrum then report and turn off the flat flag
                logging.warning('No flat spectrum found!')
                self.flat_flag = False
        
        # ---------------------------------------------------------------------
        # Spectrometer ILS
        # ---------------------------------------------------------------------
                
        # Import measured ILS
        if ils_type == 'File':
            logging.info('Importing ILS')
    
            try:
                # Read in measured ILS shape
                x_ils, y_ils = np.loadtxt(ils_path, unpack=True)
    
                # Interpolate the measured ILS onto the model grid spacing
                grid_ils = np.arange(x_ils[0], x_ils[-1], model_spacing)
                ils = griddata(x_ils, y_ils, grid_ils, 'cubic')
                self.ils = ils / np.sum(ils)
                self.generate_ils = False
    
            except OSError:
                logging.error(f'{ils_path} file not found!')
                
        # Import ILS params
        if ils_type == 'Params':
            logging.info('Importing ILS parameters')
            try:
                # Import ils parameters
                ils_params = np.loadtxt(ils_path, unpack=True)
    
                # Build ILS
                self.ils = make_ils(model_spacing, *ils_params)
                
                self.generate_ils = False
                logging.info('ILS imported')
    
            except OSError:
                logging.error(f'{ils_path} file not found!')
                
        # Manually set the ILS params
        if ils_type == 'Manual':
            self.generate_ils = True
        
        # ---------------------------------------------------------------------
        # Import solar spectrum
        # ---------------------------------------------------------------------
            
        # Import solar reference spectrum
        logging.info('Importing solar reference spectrum...')
        sol_x, sol_y = np.loadtxt(frs_path, unpack=True)
    
        # Interpolate onto model_grid
        self.frs = griddata(sol_x, sol_y, model_grid, method='cubic')
    
        logging.info('Solar reference spectrum imported')
        
        # ---------------------------------------------------------------------
        # Import Gas spectra
        # ---------------------------------------------------------------------

        logging.info('Importing gas cross-sections...')
    
        # Create an empty dictionary to hold the gas cross-sections
        self.xsecs = {}
        
        # Cycle through the parameters
        for name, param in self.params.items():
            
            # If a parameter has a xpath defined, read it in
            if param.xpath != None:                
                logging.info(f'Importing {name} reference spectrum...')
    
                # Read in the cross-section
                x, xsec = np.loadtxt(param.xpath, unpack=True)
    
                # Interpolate onto the model grid
                self.xsecs[name] = griddata(x, xsec, model_grid, 
                                            method='cubic')
        
                logging.info(f'{name} cross-section imported')
        
        # ---------------------------------------------------------------------
        # Other model settings
        # ---------------------------------------------------------------------
            
        self.fit_window = fit_window
        self.stray_window = stray_window
        self.stray_flag = stray_flag
        self.flat_flag = flat_flag
        self.dark_flag = dark_flag
        self.model_spacing = model_spacing    

# =============================================================================
#   Spectrum Pre-processing
# =============================================================================

    def pre_process(self, spectrum):

        '''
        Function to pre-process the measured spectrum to prepare it for the
        fit, correcting for the dark and flat spectrum, stray light and
        extracting the fit wavelength window
        Parameters
        ----------
        spectrum : 2D numpy array
            The spectrum in [wavelength, intensities]
        common : dict
            Common dictionary of parameters and variables passed from the main
            program to subroutines
        Returns
        -------
        processed_spec : 2D numpy array
            The processed spectrum, corrected for dark, flat, stray light and
            cut to the desired wavelength window
        '''

        # Unpack spectrum
        x, y = spectrum

        # Remove the dark spectrum from the measured spectrum
        if self.dark_flag == True:
            try:
                y = np.subtract(y, self.dark_spec)
            except ValueError:
                logging.exception('Error in dark correction. Is dark spectrum'+
                                  ' the same shape as the measurement?')

        # Remove stray light
        if self.stray_flag == True:
            stray_idx = np.where(np.logical_and(x >= self.stray_window[0],
                                                x <= self.stray_window[1]))

            if len(stray_idx[0]) == 0:
                logging.warn('No stray window outside spectrum, disabling ' + 
                             'stray correction')
                self.stray_flag = False
            
            else: y = np.subtract(y, np.average(y[stray_idx]))

        # Cut desired wavelength window
        fit_idx = np.where(np.logical_and(x >= self.fit_window[0],
                                          x <= self.fit_window[1]))
        grid = x[fit_idx]
        spec = y[fit_idx]

        # Divide by flat spectrum
        if self.flat_flag == True:
            try:
                spec = np.divide(spec, self.flat)
            except ValueError:
                logging.exception('Error in flat correction. Is flat spectrum'+
                                  ' the same shape as the measurement?')
                

        return np.row_stack([grid, spec])
        
# =============================================================================
#   fit_spectrum
# =============================================================================

    def fit_spectrum(self, spectrum, update_params=False, resid_limit=None,
                     resid_type='Percentage', int_limit=None, calc_od=[],
                     pre_process=True, interp_method='cubic'):

        '''
        Fit the supplied spectrum using a non-linear least squares minimisation
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
        interp_method : str, optional, default="cubic"
            Controls whether the interpolation at the end of the forward model
            is cubic or linear. Must be either "cubic", "linear" or "nearest".
            See scipy.interpolate.griddata for details.
        Returns
        -------
        fit_result : ifit.spectral_analysis.FitResult object
            An object that contains the fit results
        '''

        # Check is spectrum requires preprocessing
        if pre_process:
            spectrum = self.pre_process(spectrum)

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
        fit_result = FitResult(spectrum, popt, perr, nerr, self.fwd_model,
                               self.params, resid_type)

        if nerr:

            # Update the initial guess
            if update_params:
                reset_flag = False

                # Check the residual level
                if resid_limit != None:
                    if fit_result.resid.max()>resid_limit:
                        logging.info('High residual detected')
                        reset_flag = True
                        fit_result.nerr = 2

                # Check for spectrum light levels
                if int_limit != None:
                    
                    # Check for low intensity
                    if min(spec) <= int_limit[0]:
                        logging.info('Low intensity detected')
                        reset_flag = True
                        fit_result.nerr = 2

                    # Check for high intensity
                    elif max(spec) >= int_limit[1]:
                        logging.info('High intensity detected')
                        reset_flag = True
                        fit_result.nerr = 2

                # Reset the parameters if the fit was ok
                if reset_flag:
                    logging.info('Resetting initial guess parameters')
                    self.p0 = self.params.fittedvalueslist()
                else:
                    self.p0 = popt

            # Calculate the Optical depth spectra
            for par in calc_od:
                if par in self.params:
                    fit_result.calc_od(par, self)

        else:
            logging.warn(f'Fit failed!')
            if update_params:
                self.p0 = self.params.fittedvalueslist()
            for par in calc_od:
                if par in self.params:
                    fit_result.meas_od[par] = np.full(len(spec), np.nan)
                    fit_result.synth_od[par] = np.full(len(spec), np.nan)


        return fit_result

# =============================================================================
#   Forward Model
# =============================================================================

    def fwd_model(self, x, *p0):

        '''
        iFit forward model to fit measured UV sky spectra:
        I(w) = ILS *conv* {I_off(w) + I*(w) x P(w) x exp( SUM[-xsec(w) . amt])}
        where w is the wavelength.
        Requires the following to be defined in the common dictionary:
            - params:       Parameters object holding the fit parameters
            - model_grid:   The wavelength grid on which the forward model is
                            built
            - frs:          The Fraunhofer reference spectrum interpolated onto
                            the model_grid
            - xsecs:        Dictionary of the absorber cross sections that have
                            been pre-interpolated onto the model grid.
                            Typically includes all gas spectra and the Ring
                            spectrum
            - generate_ils: Boolian flag telling the function whether to
                            build the ILS or not. If False then the ILS
                            must be predefined in the common
            - ils           The instrument line shape of the spectrometer. Only
                            used if generate ILS is False.
        Parameters
        ----------
        grid, array
            Measurement wavelength grid
        *x0, list
            Forward model state vector. Should consist of:
                - bg_polyx: Background polynomial coefficients
                - offsetx:  The intensity offset polynomial coefficients
                - shiftx:   The wavelength shift polynomial
                - gases:    Any variable with an associated cross section,
                            including absorbing gases and Ring. Each "gas" is
                            converted to transmittance through:
                                      gas_T = exp(-xsec . amt)
                For polynomial parameters x represents ascending intergers
                starting from 0 which correspond to the decreasing power of
                that coefficient
        Returns
        -------
        fit, array
            Fitted spectrum interpolated onto the spectrometer wavelength grid
        '''

        # Get dictionary of fitted parameters
        params = self.params
        p = params.valuesdict()

        # Update the fitted parameter values with those supplied to the forward
        # model
        i = 0
        for par in params.values():
            if par.vary:
                p[par.name] = p0[i]
                i+=1
            else:
                p[par.name] = par.value

        # Unpack polynomial parameters
        bg_poly_coefs = [p[n] for n in p if 'bg_poly' in n]
        offset_coefs  = [p[n] for n in p if 'offset' in n]
        shift_coefs   = [p[n] for n in p if 'shift' in n]

        # Construct background polynomial
        bg_poly = np.polyval(bg_poly_coefs, self.model_grid)
        frs = np.multiply(self.frs, bg_poly)

        # Create empty array to hold optical depth spectra
        gas_T = np.zeros((len(self.xsecs),
                          len(self.model_grid)))

        # Calculate the gas optical depth spectra
        for n, gas in enumerate(self.xsecs):
            gas_T[n] = (np.multiply(self.xsecs[gas], p[gas]))

        # Sum the gas ODs
        sum_gas_T = np.sum(gas_T, axis=0)

        # Build the exponent term
        exponent = np.exp(-sum_gas_T)

        # Build the baseline polynomial
        offset = np.polyval(offset_coefs, self.model_grid)

        # Build the complete model
        raw_F = np.multiply(frs, exponent) + offset

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
        wl_shift = np.polyval(shift_coefs, self.model_grid)
        shift_model_grid = np.add(self.model_grid, wl_shift)

        # Interpolate onto measurement wavelength grid
        fit = griddata(shift_model_grid, F_conv, x,
                       method=self.interp_method)

        return fit

# =============================================================================
# =============================================================================
# # Fit Result
# =============================================================================
# =============================================================================

class FitResult():

    '''
    Contains the fit results including:
        - params:    the Parameters object with the fitted values
        - grid:      the cut wavelength window
        - spec:      the cut intensity spectrum
        - popt:      the optimised fit results (as a list)
        - perr:      the fit errors (as a list)
        - nerr:      the error code of the fit. 1 if successful, 0 if failed
        - fwd_model: the forward model used to complete the fit
        - meas_od:   dictionary of the measured optical depths
        - synth_od:  dictionary of the synthetic optical depths
        - fit:       the final fitted spectrum
        - resid:     the fit residual
    Parameters
    ----------
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
    common : dictionary
        Common dictionary of program variables
    resid_type : str
        Controls how the fit residual is calculated:
            - 'Percentage': calculated as (spec - fit) / spec * 100
            - 'Absolute':   calculated as spec - fit
    '''

    def __init__(self, spectrum, popt, perr, nerr, fwd_model, params,
                 resid_type):

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
                par.set(fit_val = self.popt[n],
                        fit_err = self.perr[n])
                n += 1
            else:
                par.set(fit_val = par.value)
                par.set(fit_err = 0)

        # If fit was successful then calculate the fit and residual
        if self.nerr:

            # Generate the fit
            self.fit = self.fwd_model(self.grid, *self.popt)

            # Calculate the residual
            if resid_type == 'Absolute':
                self.resid = self.spec - self.fit
            elif resid_type == 'Percentage':
                self.resid = (self.spec - self.fit)/self.spec * 100

        # If not then return nans
        else:
            self.fit = np.full(len(self.spec), np.nan)
            self.resid = np.full(len(self.spec), np.nan)

# =============================================================================
#   Calculate Optical Depths
# =============================================================================

    def calc_od(self, par_name, analyser):
        '''Calculates the optical depth for the given parameter'''

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
        wl_shift = np.polyval(shift_coefs, analyser.model_grid)
        shift_model_grid = np.add(analyser.model_grid, wl_shift)

        # Calculate the wavelength offset
        offset_coefs = [p[n] for n in p if 'offset' in n]
        offset = np.polyval(offset_coefs, analyser.model_grid)
        offset = griddata(shift_model_grid,
                          offset,
                          self.grid,
                          method = 'cubic')

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
                          method = 'cubic')

        # Add to self
        self.meas_od[par_name] = -np.log(np.divide(self.spec-offset, fit))
        self.synth_od[par_name] = par_od