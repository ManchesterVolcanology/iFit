# -*- coding: utf-8 -*-
"""
Created on Mon May  8 09:07:32 2017

@author: mqbpwbe2
"""

import logging
import numpy as np
from scipy.interpolate import griddata

from ifit.make_ils import make_ils

#==============================================================================
#================================ model_setup =================================
#==============================================================================

def model_setup(settings):

    '''
    Function to set up the model, including setting the model grid, reading in
    reference spectra and building the spectrometer ILS

    **Parameters**

    settings : dict
        Dictionary of model settings containing:
            - w_lo:          float, the lower limit of the fitting window
            - w_hi:          float, the upper limit of the fitting window
            - model_spacing: float, the spacing of the model wavelength grid
            - model_padding: float, the padding of the model wavelength grid
            - spec_name:     str, the serial number of the spectrometer. Used
                               to identify the correct ILS and flat files
            - flat_flag:     bool, if true the the flat spectrum is imported
            - ils_type:      str, controls how the ILS is build:
                                 ~ "File": import the ils params form file
                                 ~ "Shape": import the actual measured ILS
                                 ~ "Fix": Use a user defined width
            - sol_path:      str, path to the solar spectrum
            - ring_path:     str, path to the Ring spectrum
            - gas_data:      dict, contains the gases to include in the fit and
                               the path to each reference cross-section

    **Returns**

    common, dict
        Common dictionary to be passed to the fitting function

    '''

    common = {}

#============================== Build model grid ==============================

    # Build model grid, a high res grid on which the forward model is build
    start = settings['w_lo'] - settings['model_padding']
    stop = settings['w_hi'] + settings['model_padding'] + settings['model_spacing']
    model_grid = np.arange(start, stop, step = settings['model_spacing'])

    # Add the model grid to the common dict
    common['model_grid'] = model_grid

#============================ Import flat spectrum ============================

    # Try importing flat spectrum
    if settings['flat_flag'] == True:

        logging.info('Importing flat spectrum')

        try:
            # Import flat spectrum and extract the fit window
            flat_x, flat_y = np.loadtxt(settings['flat_path'] , unpack = True)
            f_idx = np.where(np.logical_and(settings['w_lo'] <= flat_x,
                                            flat_x <= settings['w_hi']))
            common['flat'] = flat_y[f_idx]

            logging.info('Flat spectrum imported')

        except OSError:
            # If no flat spectrum then report and turn off the flat flag
            logging.warning('No flat spectrum found!')
            common['flat_flag'] = False

#============================== Import ILS ====================================

    # Try importing ils
    if settings['ils_type'] == 'File':
        logging.info('Importing ILS parameters')
        try:
            # Import ils width
            FWEM, k, a_w, a_k = np.loadtxt(settings['ils_path'], unpack = True)

            # Build ILS
            common['ils_params'] = [FWEM, k, a_w, a_k]
            common['ils'] = make_ils(settings['model_spacing'],
                                     FWEM,
                                     k,
                                     a_w,
                                     a_k)
            
            common['generate_ils'] = False
            logging.info('ILS imported')

        except OSError:
            logging.error(f'No ILS file found for {settings["spec_name"]}!')

    elif settings['ils_type'] == 'Shape':
        logging.info('Importing ILS')

        try:

            # Read in measured ILS shape
            x_ils, y_ils = np.loadtxt(settings['ils_path'], unpack = True)

            # Interpolate the measured ILS onto the model grid spacing
            grid_ils = np.arange(x_ils[0], x_ils[-1], settings['model_res'])
            common['ils'] = griddata(x_ils, y_ils, grid_ils, 'cubic')
            common['ils'] = common['ils'] / np.sum(common['ils'])
            common['generate_ils'] = False

        except OSError:
            logging.error(f'No ILS file found for {settings["spec_name"]}!')

    elif settings['ils_type'] == 'Fix':
        logging.info('Building ILS...')

        # Build the ILS
        common['ils'] = make_ils(settings['model_spacing'],
                                 settings['ils_width'])

        logging.info('ILS built')
        
    # Otherwise tell the program to build the ILS at the model.
    else: common['generate_ils'] = True

#=========================== Import solar spectrum ============================

    # Import solar reference spectrum
    logging.info('Importing solar reference spectrum...')
    sol_x, sol_y = np.loadtxt(settings['frs_path'], unpack = True)

    # Interpolate onto model_grid
    common['frs'] = griddata(sol_x, sol_y, model_grid, method = 'cubic')

    logging.info('Solar reference spectrum imported')

#============================= Import gas spectra =============================

    logging.info('Importing gas cross-sections...')

    # Create an empty dictionary to hold the gas cross-sections
    xsecs = {}

    for gas in settings['gas_data']:

        logging.info(f'Importing {gas} reference spectrum...')

        # Read in the cross-section
        x, xsec = np.loadtxt(settings['gas_data'][gas][0], unpack = True)

        # Interpolate onto the model grid
        xsecs[gas] = griddata(x, xsec, model_grid, method = 'cubic')

        logging.info(f'{gas} cross-section imported')

    # Add the cross-sections to the common dict
    common['xsecs'] = xsecs

#=============================== Other settings ===============================

    common['stray_flag'] = settings['stray_flag']
    common['flat_flag'] = settings['flat_flag']
    common['dark_flag'] = settings['dark_flag']
    common['model_spacing'] = settings['model_spacing']

    return common
