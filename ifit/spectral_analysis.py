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
from collections import OrderedDict

#==============================================================================
#================================ fit_spectrum ================================
#==============================================================================

def fit_spectrum(spectrum, common, update_params=False, **kwargs):

    # Make a copy of the parameters
    params = common['params'].make_copy()

    # Unpack the spectrum
    grid, spec = spectrum

    # Get the fit first guesses
    fit_params = params.fittedvalueslist()

    # Perform the fit!
    logging.debug('Begin fit')
    results = leastsq(residual, fit_params, args = (grid, spec, common),
                      full_output = True, **kwargs)
    logging.debug('Fit complete')

    # Unpack the fit results
    popt = results[0]
    pcov = results[1]
    info = results[2]
    mesg = results[3]
    nerr = results[4]

    # Check for a fit error
    if nerr in [1,2,3,4]:

        # Calculate the error from the partial covarience matrix
        if (len(spec) > len(fit_params)) and pcov is not None:
            s_sq = (residual(popt, grid, spec, common)**2).sum()/ \
                   (len(spec)-len(fit_params))
            pcov = pcov * s_sq
        else:
            pcov = np.inf
        try:
            perr = np.sqrt(np.diag(pcov))
        except ValueError:
            logging.warning('Unable to calculate error matrix')
            perr = np.full(len(popt), np.nan)

        # Add the fit results to each parameter
        n = 0
        for par in params.values():

            if par.vary:
                par.set(fit_val = popt[n],
                        fit_err = perr[n])
                n += 1
            else:
                par.set(fit_val = par.value)
                par.set(fit_err = 0)

        # Recalulate the model
        fit = ifit_fwd_model(grid, *popt, **common)

        # Calculate the residual
        resid = (spec-fit)/spec * 100

        # Calculate the SO2 OD
        #self.calc_od('SO2_295', spectrum, common)

        # Use the fitted parameters as the next first guess
        if update_params:
            common['params'].update_values(popt)

    else:
        logging.warning(f'Fit failed: {mesg}')
        fit = np.full(len(grid), np.nan)
        resid = np.full(len(grid), np.nan)

    fit_result = FitResult(params, popt, pcov, info, mesg, nerr, perr, fit, resid)

    return fit_result



#==============================================================================
#================================ FitSpectrum =================================
#==============================================================================

class FitResult(OrderedDict):

    def __init__(self, params, popt, pcov, info, mesg, nerr, perr, fit, resid):

        self.params = params
        self.popt = popt
        self.pcov = pcov
        self.info = info
        self.mesg = mesg
        self.nerr = nerr
        self.perr = perr
        self.fit = fit
        self.resid = resid

    def print_result(self):
        '''Prints the result of the fit in a readable way'''
        msg = self.params.pretty_print()

        msg += f'\nNumber of function calls: {self.info["nfev"]}'

        msg += f'\n\nAverage of Residual: {np.average(self.resid):.3g}'
        msg += f'\nStdev of Residual:   {np.std(self.resid):.3g}'

        return msg

class FitSpectrum(OrderedDict):
    '''
    An object to hold the fit results in a concise and useful way

    **Parameters**

    spectrum : 2D numpy array
        The spectrum to be fitted as [wavelength, intensity]

    common : dict
        The common dictionary holding the model information and settings

    update_params : bool, optional, default = False

    kwargs : tuple, optional
        Keywork arguments to be passed to the leastsq
    '''

    def __init__(self, spectrum, common, update_params = False, **kwargs):

        self.params = common['params'].make_copy()

        # Unpack the spectrum
        grid, spec = spectrum

        # Get the fit first guesses
        fit_params = self.params.fittedvalueslist()

        # Perform the fit!
        logging.debug('Begin fit')
        results = leastsq(residual, fit_params, args = (grid, spec, common),
                          full_output = True, **kwargs)
        logging.debug('Fit complete')

        # Unpack the fit results
        self.popt = results[0]
        self.pcov = results[1]
        self.info = results[2]
        self.mesg = results[3]
        self.nerr = results[4]

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

            # Calculate the SO2 OD
            #self.calc_od('SO2_295', spectrum, common)

            # Use the fitted parameters as the next first guess
            if update_params:
                common['params'].update_values(self.popt)

        else:
            logging.warning(f'Fit failed: {self.mesg}')
            self.fit = np.full(len(grid), np.nan)
            self.resid = np.full(len(grid), np.nan)


    def calc_od(self, par, spectrum, common):
        '''Calculate the measured and synthetic OD of a give gas'''

        # Create a copy of the common
        c = copy.deepcopy(common)

        # Unpack the spectrum
        grid, spec = spectrum

        # Set the desired parameter to zero and fix it
        c['params'][par].set(value=0, vary=False)

        # Calculate the fit without the parameter
        sim_params = c['params'].fittedvalueslist()
        p = c['params'].valuesdict()
        fit_no_par = ifit_fwd_model(grid, *sim_params, **c)

        # Calculate a shifted model grid
        shift_model_grid = np.add(c['model_grid'], p['shift'])
        line = np.linspace(0, 1, num = len(shift_model_grid))
        shift_model_grid = np.add(shift_model_grid,
                                  np.multiply(line, p['stretch']))

        # Calculate the parameter od
        par_od = np.exp(-(np.multiply(c['gas_xsec'][par], p[par])))

        # Convolve with the ILS and interpolate onto the measurement grid
        par_od = griddata(shift_model_grid,
                          np.convolve(par_od, c['ils'], 'same'),
                          grid,
                          method = 'cubic')

        # Add to self
        self.meas_od = fit_no_par#-np.log(np.divide(spec, fit_no_par))
        self.synth_od = spec#-np.log(par_od)


    def print_result(self):
        '''Prints the result of the fit in a readable way'''
        msg = self.params.pretty_print()

        msg += f'\nNumber of function calls: {self.info["nfev"]}'

        msg += f'\n\nAverage of Residual: {np.average(self.resid):.3g}'
        msg += f'\nStdev of Residual:   {np.std(self.resid):.3g}'

        return msg

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

    I = ILS *conv* { I* x P x exp( SUM[-xsec.gas_SCD] + Ring)}

    **Parameters**

    grid, array
        Measurement wavelength grid

    *x0, list
        Forward model state vector

    **com, dict
        The common dictionary containing the informations required for the fit

    **Returns**

    fit, array
        Fitted spectrum interpolated onto the spectrometer wavelength grid
    '''

    # Check for nans
    if np.isnan(x0).any():
        return np.zeros(len(meas_grid))

    # Create empty dictionary to hold the parameters and a counter
    p = {}
    i = 0

    # Unpack the fit paramters
    for par in com['params'].values():

        # If the parameter is varying pull it from fit_params
        if par.vary:
            p[par.name] = x0[i]
            i += 1
        else:
            # If it is fixed pull its value from the Parameters object
            p[par.name] = par.value

    # Unpack polynomial parameters
    bg_poly_coefs = []
    offset_coefs = []
    shift_coefs = []
    for n, par in enumerate(com['params'].values()):
        if 'bg_poly' in par.name:
            bg_poly_coefs.append(x0[n])
        elif 'offset' in par.name:
            offset_coefs.append(x0[n])
        elif 'shift' in par.name:
            shift_coefs.append(x0[n])

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
    raw_F = np.multiply(frs, exponent) + bl_poly

    # Apply the ILS convolution
    ils = com['ils']

    F_conv = np.convolve(raw_F, ils, 'same')

    # Apply shift and stretch to the model_grid
    wl_shift = np.polyval(shift_coefs, line)
    shift_model_grid = np.add(com['model_grid'], wl_shift)

    # Interpolate onto measurement wavelength grid
    fit = griddata(shift_model_grid, F_conv, meas_grid, method = 'cubic')

    return fit

