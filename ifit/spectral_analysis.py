#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 11:45:29 2019

@author: mqbpwbe2
"""

import copy
import logging
import numpy as np
from scipy.interpolate import griddata
from scipy.optimize import leastsq

#==============================================================================
#================================ fit_spectrum ================================
#==============================================================================

def fit_spectrum(spectrum, common, update_params=False, resid_limit=5,
                 calc_od=[],**kwargs):

    # Unpack the spectrum
    grid, spec = spectrum

    # Get the fit first guesses
    fit_params = common['params'].fittedvalueslist()

    # Perform the fit!
    logging.debug('Begin fit')
    fit_results = leastsq(residual, fit_params, args = (grid, spec, common),
                          full_output = True, **kwargs)
    logging.debug('Fit complete')

    # Format the fit results into a FitResult object
    fit_result = FitResult(fit_results, 
                           common, 
                           spectrum, 
                           update_params, 
                           resid_limit)

    # Calculate the Optical depth spectra
    for par in calc_od:
        fit_result.calc_od(par, common)

    return fit_result



#==============================================================================
#================================ FitSpectrum =================================
#==============================================================================

class FitResult():

    '''
    Object containing the results of the fit

    **Parameters**

    fit_results : list
        The output of scipy.optimize.leastsq

    common : dict
        The common dictionary containing information for the fit

    spectrum : 2D numpy array
        The measured spectrum to which the fit is done as:
        [wavelength, intensity]

    update_params : bool, optional, default = False
        If True then the first guess of the fit parameters is updated with the
        previous results

    **Attributes**

    params : ifit.parameters.Parameters object
        Contains the fit parameters as they were for this fit

    popt : list
        The optimised fit parameters given by the fitting function

    fit : numpy array
        The final output of the forward model for the fitted spectrum

    resid : numpy array
        The percentage difference between the fit and the measured spectrum.
        Calculated as (spec - fit) / spec * 100

    meas_od : dict
        Dictionary of the "measured" optical depths for each parameter.
        Calculated by dividing the measured spectrum by the forward model
        without the selected parameter and taking the negative natural log

    synth_od : dict
        Dictionary of the "synthetic" optical depths for each parameter.
        Calculated by multiplying the "gas" xsec by the fitted "gas" amt and
        interpolating onto the measurement grid

    pcov : array
        The covarience array given by the fitting function. Used to calculate
        the fit error for each parameter

    info : dictionary
        Contains the other outputs of the fitting function. See
        scipy.optimize.leastsq for more details

    mesg : str
        String message giving the cause of the fit failure.

    nerr : int
        Integer code showing the success or failure of the fit. If equal to
        1, 2, 3 or 4 the fit was successful. Otherwise the fit failled. More
        information is given in mesg


    '''

    def __init__(self, fit_results, common, spectrum, update_params=False,
                 resid_limit=5):

        # Create a copy of the parameters
        self.params = common['params'].make_copy()

        # Get the fit first guesses
        fit_params = self.params.fittedvalueslist()

        # Unpack the measured spectrum
        self.spectrum = spectrum
        grid, spec = spectrum

        # Unpack the fit results
        self.popt = fit_results[0]
        self.pcov = fit_results[1]
        self.info = fit_results[2]
        self.mesg = fit_results[3]
        self.nerr = fit_results[4]

        # Check for a fit error
        if self.nerr in [1,2,3,4]:

            # Calculate the error from the partial covarience matrix
            if (len(spec) > len(fit_params)) and self.pcov is not None:
                s_sq = (residual(self.popt, grid, spec, common)**2).sum()/ \
                       (len(spec)-len(fit_params))
                self.pcov = self.pcov * s_sq
            else:
                self.pcov = np.inf
            try:
                self.perr = np.sqrt(np.diag(self.pcov))
            except ValueError:
                logging.warning('Unable to calculate error matrix')
                self.perr = np.full(len(self.popt), np.nan)

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

            # Recalulate the model
            self.fit = ifit_fwd_model(grid, *self.popt, **common)

            # Calculate the residual
            self.resid = (spec-self.fit)/spec * 100
            
            # Check if the spectrum was saturated
            if max(spec) >= 35000:
                logging.info(f'Spectrum saturating, resetting fit parameters')
                common['params'].update_values(common['x0'])
                self.nerr = 5

            # Use the fitted parameters as the next first guess
            if update_params:
                if max(self.resid) < resid_limit:
                    common['params'].update_values(self.popt)
                else:
                    logging.info(f'Low fit quality, resetting fit parameters')
                    common['params'].update_values(common['x0'])
                    self.nerr = 6
            

        else:
            logging.warning(f'Fit failed: {self.mesg}')
            self.fit = np.full(len(grid), np.nan)
            self.resid = np.full(len(grid), np.nan)


        self.meas_od = {}
        self.synth_od = {}


    def print_result(self):
        '''Prints the result of the fit in a readable way'''

        msg = self.params.pretty_print()

        msg += f'\nNumber of function calls: {self.info["nfev"]}'

        msg += f'\n\nAverage of Residual: {np.average(self.resid):.3g}'
        msg += f'\nStdev of Residual:   {np.std(self.resid):.3g}'

        return msg


    def calc_od(self, par_name, common):
        '''Calculates the optical depth for the given parameter'''

        # Unpack the measured spectrum
        grid, spec = self.spectrum

        # Create a copy of the common
        c = copy.deepcopy(common)

        # Make a copy of the parameters
        params = self.params.make_copy()

        # Set the parameter to zero
        params[par_name].set(fit_val=0)

        # Calculate the fit without the parameter
        fit_params = params.popt_list()
        p = params.valuesdict()

        fit = ifit_fwd_model(grid, *fit_params, **c)

        # Calculate a shifted model grid
        #shift_model_grid = np.add(c['model_grid'], p['shift'])
        #line = np.linspace(0, 1, num = len(shift_model_grid))
        #shift_model_grid = np.add(shift_model_grid,
        #                          np.multiply(line, p['stretch']))

        # Calculate the parameter od
        par_od = np.multiply(c['xsecs'][par_name], p[par_name])

        # Convolve with the ILS and interpolate onto the measurement grid
        par_od = griddata(c['model_grid'],
                          np.convolve(par_od, c['ils'], 'same'),
                          grid,
                          method = 'cubic')

        # Add to self
        self.meas_od[par_name] = -np.log(np.divide(spec, fit))
        self.synth_od[par_name] = par_od

#==============================================================================
#================================ pre_process =================================
#==============================================================================

def pre_process(spectrum, common):
    '''
    Function to pre-process the measured spectrum to prepare it for the fit,
    correcting for the dark and flat spectrum, stray light and extracting the
    fit wavelength window

    **Parameters**

    spectrum : 2D numpy array
        The spectrum in [wavelength, intensities]

    common : dict
        Common dictionary of parameters and variables passed from the main
        program to subroutines

    **Returns**

    processed_spec : 2D numpy array
        The processed spectrum, corrected for dark, flat, stray light and cut
        to the desired wavelength window
    '''

    # Unpack spectrum
    x, y = spectrum

    # Remove the dark spectrum from the measured spectrum
    if common['dark_flag'] == True:
        y = np.subtract(y, common['dark'])

    # Remove stray light
    if common['stray_flag'] == True:
        y = np.subtract(y, np.average(y[common['stray_idx']]))

    # Cut desired wavelength window
    grid = x[common['fit_idx']]
    spec = y[common['fit_idx']]

    # Divide by flat spectrum
    if common['flat_flag'] == True:
        spec = np.divide(spec, common['flat'])

    return np.row_stack([grid, spec])

#==============================================================================
#================================== residual ==================================
#==============================================================================

def residual(fit_params, *args):
    '''
    Function to calculate the residual between the forward model and measured
    data
    '''

    grid, data, common = args

    fit = ifit_fwd_model(grid, *fit_params, **common)

    return np.subtract(fit, data)

#==============================================================================
#=============================== ifit_fwd_model ===============================
#==============================================================================

def ifit_fwd_model(meas_grid, *x0, **com):

    '''
    iFit forward model to fit measured UV sky spectra:

    I(w) = ILS *conv* { I_offset(w) + I*(w) x P(w) x exp( SUM[-xsec(w) . amt])}

    where w is the wavelength.

    **Parameters**

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

            For polynomial parameters x represents ascending intergers starting
            from 0 which correspond to the decreasing power of that coefficient

    **com, dict
        The common dictionary containing the information required for the fit.
        Must contain:
            - params: Parameters object holding the fit parameters
            - model_grid: The wavelength grid on which the forward model is
                          built
            - frs:        The Fraunhofer reference spectrum interpolated onto
                          model_grid
            - xsecs:      Dictionary of the absorber cross sections that have
                          been pre-interpolated onto the moddel grid. Typically
                          includes all gas spectra and the Ring spectrum
            - ils         The instrument line shape of the spectrometer


    **Returns**

    fit, array
        Fitted spectrum interpolated onto the spectrometer wavelength grid
    '''

    # Check for nans
    if np.isnan(x0).any():
        return np.zeros(len(meas_grid))

    # Get dictionary of fitted parameters
    p = com['params'].fittedvaluesdict()

    # Update the fitted parameter values with those supplied to the forward
    # model
    for i, par in enumerate(p):
        p[par] = x0[i]

    # Unpack polynomial parameters
    bg_poly_coefs = [p[n] for n in p if 'bg_poly' in n]
    offset_coefs  = [p[n] for n in p if 'offset' in n]
    shift_coefs   = [p[n] for n in p if 'shift' in n]

    # Construct background polynomial
    bg_poly = np.polyval(bg_poly_coefs, com['model_grid'])
    frs = np.multiply(com['frs'], bg_poly)

    # Create empty array to hold transmission spectra
    gas_T = np.zeros((len(com['xsecs']), len(com['model_grid'])))

    # Calculate the gas transmission spectra
    for n, gas in enumerate(com['xsecs']):
        gas_T[n] = -(np.multiply(com['xsecs'][gas], p[gas]))

    # Sum the gas transmissions
    sum_gas_T = np.sum(gas_T, axis=0)

    # Build the exponent term
    exponent = np.exp(sum_gas_T)

    # Build the baseline polynomial
    line = np.linspace(0, 1, num = len(com['model_grid']))
    bl_poly = np.polyval(offset_coefs, line)

    # Build the complete model
    raw_F = np.add(np.multiply(frs, exponent), bl_poly)

    # Apply the ILS convolution
    ils = com['ils']

    F_conv = np.convolve(raw_F, ils, 'same')

    # Apply shift and stretch to the model_grid
    wl_shift = np.polyval(shift_coefs, line)
    shift_model_grid = np.add(com['model_grid'], wl_shift)

    # Interpolate onto measurement wavelength grid
    fit = griddata(shift_model_grid, F_conv, meas_grid, method = 'cubic')

    return fit

def ifit_od_fwd_model(meas_grid, *x0, **com):

    # Get dictionary of fitted parameters
    p = com['params'].fittedvaluesdict()

    # Update the fitted parameter values with those supplied to the forward
    # model
    for i, par in enumerate(p):
        p[par] = x0[i]

    # Unpack polynomial parameters
    offset_coefs  = [p[n] for n in p if 'offset' in n]
    shift_coefs   = [p[n] for n in p if 'shift' in n]














