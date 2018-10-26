# -*- coding: utf-8 -*-
"""
Created on Mon May  8 09:07:32 2017

@author: mqbpwbe2
"""

import numpy as np
from scipy.interpolate import griddata

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

    # Build model grid, a high res grid on which the forward model is build. It extends
    #  3 nm beyond the measurement grid and has a spacing controlled by the user
    start = settings['wave_start'] - settings['model_pad']
    stop = settings['wave_stop'] + settings['model_pad'] + float(settings['model_res'])
    model_grid = np.arange(start, stop, step = float(settings['model_res']))

    common['model_grid'] = model_grid
    
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
    
    # Try importing ils
    if settings['Fit ILS'] == 'File':
        self.print_output('Importing ILS width', add_line = False)
        try:
            # Import ils width
            ils_path = 'data_bases/Spectrometer/ils_width_'+common['spec_name']+'.txt'
            common['params']['ils_width'][0] = np.loadtxt(ils_path)           
            
            self.print_output('ILS width imported', add_line = False)
            
        except OSError:
            self.print_output('No ILS file found', add_line = False)  
    
    if settings['Fit ILS'] == 'Shape':
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
            
    
    # Import solar reference spectrum
    self.print_output('Importing solar reference spectrum...', add_line = False)
    sol_x, sol_y = np.loadtxt(settings['sol_path'], unpack = True)
    
    # Normalise to typical intensity for simplicity
    sol_y = np.divide(sol_y, (sol_y.max() / 70000))
    
    # Smooth and interpolate onto model_grid
    common['sol'] = griddata(sol_x, sol_y, model_grid, method = 'cubic')
        
    self.print_output('Solar reference spectrum imported', add_line = False)
     
    # Import solar residual spectrum
    if common['solar_resid_flag'] == 'Remove':
        try:
            x, common['solar_resid'] = np.loadtxt(settings['solar_resid_path'], 
                                                  unpack=True)
            
        except FileNotFoundError:
            self.print_output('Solar Residual not found', add_line = False)
            common['solar_resid'] = 1
    
    # Import ring spectrum and interpolate onto the model_grid
    self.print_output('Importing ring spectrum...', add_line = False)
    ring_x, ring_y = np.loadtxt(settings['ring_path'], unpack = True)
    
    common['ring'] = griddata(ring_x, ring_y, model_grid, method = 'cubic')
    self.print_output('Ring spectrum imported', add_line = False)
    
    self.print_output('Importing gas cross-sections...', add_line = False)
    
    
    # Import SO2 data
    so2_xsec = np.loadtxt(settings['so2_path'], skiprows=1)
    common['so2_xsec'] = griddata(so2_xsec[:,0], so2_xsec[:,1], model_grid, 
                                  method = 'cubic')
    self.print_output('SO2 cross-section imported', add_line = False)
    

    # Import NO2 data
    no2_xsec = np.loadtxt(settings['no2_path'], skiprows=43)
    common['no2_xsec'] = griddata(no2_xsec[:,0], no2_xsec[:,2], model_grid,
                                  method = 'cubic')
    self.print_output('NO2 cross-section imported', add_line = False)  
    
    
    # Import O3 data
    o3_xsec = np.loadtxt(settings['o3_path'])
    
    # Get column number from temperature
    temps = ['298K', '283K', '273K', '263K', '253K', '243K', '233K', '223K', '213K', 
             '203K', '193K']
    col_n = temps.index(settings['o3_temp']) + 1
    common['o3_xsec'] = griddata(o3_xsec[:,0], o3_xsec[:,col_n], model_grid, 
                                 method = 'cubic')
    self.print_output('O3 cross-section imported', add_line = False)
    
    
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
