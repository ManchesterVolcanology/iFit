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

class Analyser():

    '''
    Object to perform the spectral analysis
    Parameters
    ----------
    common : dict
        Common dictionary of parameters and variables passed from the main
        program to subroutines
    '''

    def __init__(self, common):

        # Set the common dictionary
        self.common = common

        # Set the initial estimate for the fit parameters
        self.params = common['params'].make_copy()
        self.p0 = self.params.fittedvalueslist()

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
        if self.common['dark_flag'] == True:
            y = np.subtract(y, self.common['dark'])

        # Remove stray light
        if self.common['stray_flag'] == True:
            y = np.subtract(y, np.average(y[self.common['stray_idx']]))

        # Cut desired wavelength window
        grid = x[self.common['fit_idx']]
        spec = y[self.common['fit_idx']]

        # Divide by flat spectrum
        if self.common['flat_flag'] == True:
            spec = np.divide(spec, self.common['flat'])

        if 'solar_resid' in self.common.keys():
            spec = np.subtract(spec, self.common['solar_resid'])

        return np.row_stack([grid, spec])

#==============================================================================
#================================ fit_spectrum ================================
#==============================================================================

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
        int_limit : tuple of float, optional, default=None
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
        self.common['interp_method'] = interp_method

        # Unpack the spectrum
        grid, spec = spectrum

        # Fit the spectrum
        try:
            popt, pcov = curve_fit(self.fwd_model, grid, spec, self.p0,
                                   method='lm')

            # Calculate the parameter error
            perr = np.sqrt(np.diag(pcov))

            # Set the success flag
            nerr = 1

        except RuntimeError:
            popt = np.full(len(self.p0), np.nan)
            perr = np.full(len(self.p0), np.nan)
            nerr = 0

        # Put the results into a FitResult object
        fit_result = FitResult(spectrum, popt, perr, nerr, self.fwd_model,
                               self.params, self.common, resid_type)

        if nerr:

            # Update the initial guess
            if update_params:

                reset_flag = False

                # Check the residual level
                if resid_limit != None and fit_result.resid.max() > resid_limit:
                    logging.info('High residual detected')
                    reset_flag = True
                    fit_result.nerr = 2

                # Check for spectrum low light
                if int_limit != None and min(spec) <= int_limit[0]:
                    logging.info('Low intensity detected')
                    reset_flag = True

                # Check for spectrum saturation
                if int_limit != None and max(spec) >= int_limit[1]:
                    logging.info('High intensity detected')
                    reset_flag = True

                # Reset the parameters if the fit was ok
                if reset_flag:
                    logging.info('Resetting initial guess parameters')
                    self.p0 = self.params.fittedvalueslist()
                else:
                    self.p0 = popt

            # Calculate the Optical depth spectra
            for par in calc_od:
                if par in self.params:
                    fit_result.calc_od(par, self.common)

        else:
            logging.warn(f'Fit failed!')
            if update_params:
                self.p0 = self.params.fittedvalueslist()
            for par in calc_od:
                if par in self.params:
                    fit_result.meas_od[par] = np.full(len(spec), np.nan)
                    fit_result.synth_od[par] = np.full(len(spec), np.nan)


        return fit_result

#==============================================================================
#================================= fwd_model ==================================
#==============================================================================

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
        bg_poly = np.polyval(bg_poly_coefs, self.common['model_grid'])
        frs = np.multiply(self.common['frs'], bg_poly)

        # Create empty array to hold transmission spectra
        gas_T = np.zeros((len(self.common['xsecs']),
                          len(self.common['model_grid'])))

        # Calculate the gas transmission spectra
        for n, gas in enumerate(self.common['xsecs']):
            gas_T[n] = -(np.multiply(self.common['xsecs'][gas], p[gas]))

        # Sum the gas transmissions
        sum_gas_T = np.sum(gas_T, axis=0)

        # Build the exponent term
        exponent = np.exp(sum_gas_T)

        # Build the baseline polynomial
        line = np.linspace(0, 1, num = len(self.common['model_grid']))
        offset = np.polyval(offset_coefs, line)

        # Build the complete model
        raw_F = np.multiply(frs, exponent) + offset

        # Generate the ILS
        if self.common['generate_ils']:

            # Unpack ILS params
            ils = make_ils(self.common['model_spacing'],
                           p['fwem'],
                           p['k'],
                           p['a_w'],
                           p['a_k'])
        else:
            ils = self.common['ils']

        # Apply the ILS convolution
        F_conv = np.convolve(raw_F, ils, 'same')

        # Apply shift and stretch to the model_grid
        wl_shift = np.polyval(shift_coefs, line)
        shift_model_grid = np.add(self.common['model_grid'], wl_shift)

        # Interpolate onto measurement wavelength grid
        fit = griddata(shift_model_grid, F_conv, x,
                       method=self.common['interp_method'])

        return fit


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

    def __init__(self, spectrum, popt, perr, nerr, fwd_model, params, common,
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

    def calc_od(self, par_name, common):
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
        line = np.linspace(0, 1, num = len(common['model_grid']))
        wl_shift = np.polyval(shift_coefs, line)
        shift_model_grid = np.add(common['model_grid'], wl_shift)

        # Calculate the wavelength offset
        offset_coefs = [p[n] for n in p if 'offset' in n]
        offset = np.polyval(offset_coefs, line)
        offset = griddata(shift_model_grid,
                          offset,
                          self.grid,
                          method = 'cubic')

        # Calculate the parameter od
        par_od = np.multiply(common['xsecs'][par_name], p[par_name])

        # Make the ILS
        if common['generate_ils']:

            ils_params = []
            for name in ['fwem', 'k', 'a_w', 'a_k']:
                if params[name].vary:
                    ils_params.append(params[name].fit_val)
                else:
                    ils_params.append(params[name].value)

            # Unpack ILS params
            ils = make_ils(common['model_spacing'], *ils_params)
        else:
            ils = common['ils']

        # Convolve with the ILS and interpolate onto the measurement grid
        par_od = griddata(shift_model_grid,
                          np.convolve(par_od, ils, 'same'),
                          self.grid,
                          method = 'cubic')

        # Add to self
        self.meas_od[par_name] = -np.log(np.divide(self.spec-offset, fit))
        self.synth_od[par_name] = par_od