# -*- coding: utf-8 -*-
"""
Created on Mon May  8 09:07:32 2017

@author: mqbpwbe2
"""

import numpy as np
from scipy.interpolate import griddata

from ifit_lib.make_ils import make_ils

#========================================================================================
#==================================== build_fwd_data ====================================
#========================================================================================

def build_fwd_data(common, settings, self):

    '''
    Read in required spectra and cross-sections etc and place them on the model grid as
    required

    INPUTS
    ------
    common, dict
        Common dictionary including file paths to spectra and constants used in the
        program

    settings, dict
        Dictionary of GUI settings determined by the user

    OUTPUTS
    -------
    common, dict
        Inputted common array updated with required spectra etc

    '''

#=================================== Build model grid ===================================

    # Build model grid, a high res grid on which the forward model is build. It extends
    #  3 nm beyond the measurement grid and has a spacing controlled by the user
    start = settings['wave_start'] - settings['model_pad']
    stop = settings['wave_stop'] + settings['model_pad'] + float(settings['model_res'])
    model_grid = np.arange(start, stop, step = float(settings['model_res']))

    common['model_grid'] = model_grid

    # Make a Gaussian with the resolution of the model
    #gauss = make_gauss(0.01, settings['model_res'])

#================================= Import flat spectrum =================================

    # Try importing flat spectrum. If not found set to 1
    if common['flat_flag'] == True:
        self.print_output('Importing flat spectrum', add_line = False)
        try:
            # Import flat spectrum and extract window of interest
            flat_path = 'data_bases/Spectrometer/flat_' + common['spec_name'] + '.txt'
            flat_grid, flat = np.loadtxt(flat_path , unpack = True)
            flat_idx = np.where(np.logical_and(settings['wave_start'] <= flat_grid,
                                               flat_grid <= settings['wave_stop']))
            common['flat'] = flat[flat_idx]

            self.print_output('Flat spectrum imported', add_line = False)

        except OSError:
            self.print_output('No flat spectrum found', add_line = False)
            common['flat_flag'] = False

#=================================== Import ILS =========================================

    # Try importing ils
    if settings['Fit ILS'] == 'File':
        self.print_output('Importing ILS parameters', add_line = False)
        try:
            # Import ils width
            ils_path = 'data_bases/Spectrometer/ils_params_'+common['spec_name']+'.txt'
            FWEM, k, a_w, a_k = np.loadtxt(ils_path, unpack = True)
            common['params']['ils_width'][0] = FWEM

            # Build ILS
            common['ils'] = make_ils(settings['model_res'], FWEM, k, a_w, a_k)

            self.print_output('ILS built', add_line = False)

        except OSError:
            self.print_output('No ILS file found', add_line = False)

    elif settings['Fit ILS'] == 'Shape':
        self.print_output('Importing ILS', add_line = False)

        try:

            # Read in measured ILS shape
            ils_path = 'data_bases/Spectrometer/ils_' + common['spec_name'] + '.txt'
            x_ils, y_ils = np.loadtxt(ils_path, unpack = True)

            # Interpolate the measured ILS onto the model grid spacing
            grid_ils = np.arange(x_ils[0], x_ils[-1], settings['model_res'])
            common['ils'] = griddata(x_ils, y_ils, grid_ils, 'cubic')
            common['ils'] = common['ils'] / np.sum(common['ils'])

        except OSError:
            self.print_output('No ILS file found', add_line = False)

    elif settings['Fit ILS'] == 'Fix':
        self.print_output('Building ILS...', add_line = False)

        # Build the ILS
        common['ils'] = make_ils(settings['model_res'], settings['ils_width'])

        self.print_output('ILS built', add_line = False)

#================================ Import solar spectrum =================================

    # Import solar reference spectrum
    self.print_output('Importing solar reference spectrum...', add_line = False)
    sol_x, sol_y = np.loadtxt(settings['sol_path'], unpack = True)

    # Normalise to typical intensity for simplicity
    sol_y = np.divide(sol_y, (sol_y.max() / 70000))

    # Smooth and interpolate onto model_grid
    #sol_y = np.convolve(sol_y, gauss, mode = 'same')
    common['sol'] = griddata(sol_x, sol_y, model_grid, method = 'cubic')

    self.print_output('Solar reference spectrum imported', add_line = False)

#================================ Import solar residual =================================

    # Import solar residual spectrum
    if common['solar_resid_flag'] == 'Remove':
        try:
            x, common['solar_resid'] = np.loadtxt(settings['solar_resid_path'],
                                                  unpack=True)

        except FileNotFoundError:
            self.print_output('Solar Residual not found', add_line = False)
            common['solar_resid'] = 1

#================================= Import Ring spectrum =================================

    # Import ring spectrum and interpolate onto the model_grid
    self.print_output('Importing ring spectrum...', add_line = False)
    ring_x, ring_y = np.loadtxt(settings['ring_path'], unpack = True)

    #ring_y = np.convolve(ring_y, gauss, mode = 'same')
    common['ring'] = griddata(ring_x, ring_y, model_grid, method = 'cubic')
    self.print_output('Ring spectrum imported', add_line = False)

    self.print_output('Importing gas cross-sections...', add_line = False)

#================================= Import SO2 spectrum ==================================

    # Import SO2 data
    so2_grid, so2_xsec = np.loadtxt(settings['so2_path'], unpack = True)
    #so2_xsec = np.convolve(so2_xsec, gauss, mode = 'same')
    common['so2_xsec'] = griddata(so2_grid, so2_xsec, model_grid,
                                  method = 'cubic')
    self.print_output('SO2 cross-section imported', add_line = False)

#================================= Import BO2 spectrum ==================================

    # Import NO2 data
    no2 = np.loadtxt(settings['no2_path'], skiprows=43)
    no2_xsec = no2[:,2]#np.convolve(no2[:,2], gauss, mode = 'same')
    common['no2_xsec'] = griddata(no2[:,0], no2_xsec, model_grid,
                                  method = 'cubic')
    self.print_output('NO2 cross-section imported', add_line = False)

#================================= Import O3 spectrum ===================================

    # Import O3 data
    o3_xsec = np.loadtxt(settings['o31_path'])

    # Get column number from temperature
    temps = ['298K', '283K', '273K', '263K', '253K', '243K', '233K', '223K', '213K',
             '203K', '193K']
    col1_n = temps.index(settings['o31_temp']) + 1
    col2_n = temps.index(settings['o32_temp']) + 1

    # Set first ozone xsec
    o3 = o3_xsec[:,col1_n]#np.convolve(o3_xsec[:,col1_n], gauss, mode = 'same')
    common['o31_xsec'] = griddata(o3_xsec[:,0], o3, model_grid,
                                  method = 'cubic')

    # Set second ozone xsec
    common['o32_xsec'] = griddata(o3_xsec[:,0], o3_xsec[:,col2_n], model_grid,
                                  method = 'cubic')

    self.print_output('O3 cross-section imported', add_line = False)

#================================= Import BrO spectrum ==================================

    # Import BrO data
    bro = np.loadtxt(settings['bro_path'], skiprows=12)

    # As BrO is in terms of wavenumber(cm^-1), need to convert to nm
    bro[:,0] = np.divide(1e7, bro[:,0])

    # Reverse array so wavelength is accending
    bro = bro[::-1]

    # Interpolate onto the grid
    common['bro_xsec'] = griddata(bro[:,0], bro[:,1], model_grid, method = 'cubic')
    self.print_output('BrO cross-section imported')

    # Turn off flag to build forward model
    self.build_model_flag = False

    return common
