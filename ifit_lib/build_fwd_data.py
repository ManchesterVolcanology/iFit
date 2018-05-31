# -*- coding: utf-8 -*-
"""
Created on Mon May  8 09:07:32 2017

@author: mqbpwbe2
"""

import numpy as np
from scipy.interpolate import griddata
from ifit_lib.find_nearest import extract_window

#========================================================================================
#=====================================build_fwd_data=====================================
#========================================================================================

# Function to build the gas cross sections and place them on a model grid that extends 
#  2 nm beyond the measurement grid

# INPUTS: common: common dictionary of parameters and constants passed from the program 
#                   to subroutines. Also contains the filepaths to the required data

def build_fwd_data(common, settings, self):

    '''
    Read in required spectra and cross-sections etc and place them on the model grid as 
    required

    INPUTS
    ------
    self: requirement from tkinter
    common: common dictionary including file paths to spectra and constants used in the
            program
    settings: dictionary of settings determined by the user
    
    OUTPUTS
    -------
    common: inputted common array updated with required spectra etc
        
    '''
    
    # Build model grid, a high res grid on which the forward model is build. It extends
    #  2 nm beyond the measurement grid and has a spacing controlled by the user
    spacing = 1/float(settings['model_resolution'])
    npts = ((common['wave_stop'] + 3) - (common['wave_start'] - 3)) * spacing
    model_grid = np.linspace(common['wave_start'] - 3, common['wave_stop'] + 3, 
                             num = npts + 1)
    
    # Try importing flat spectrum. If not found set to 1
    self.print_output('Importing flat spectrum', add_line = False)
    try:
        # Import flat spectrum and extract window of interest
        flat_grid, flat = np.loadtxt(settings['flat_path'] , unpack = True)
        x, i1, i2 = extract_window(flat_grid, common['wave_start'], common['wave_stop'])
        flat = flat[i1:i2]
        self.print_output('Flat spectrum imported', add_line = False)
        
    except FileNotFoundError:
        self.print_output('No flat spectrum found', add_line = False)
        common['flat_flag'] = False    
    
    
    # Import solar reference spectrum
    self.print_output('Importing solar reference spectrum...', add_line = False)
    sol_x, sol_y = np.loadtxt(settings['sol_path'], unpack = True)
    self.print_output('Solar reference spectrum imported', add_line = False)
    
    # Interpolate onto model_grid
    sol = griddata(sol_x, sol_y, model_grid)
    
    
    # Import solar residual spectrum
    if common['solar_resid_flag'] == 'Remove':
        try:
            common['solar_resid'] = np.loadtxt(settings['solar_resid_path'])
            
        except FileNotFoundError:
            self.print_output('Solar Residual not found', add_line = False)
            common['solar_resid'] = 1
    
    # Import ring spectrum and interpolate onto the model_grid
    self.print_output('Importing ring spectrum...', add_line = False)
    ring_x, ring_y = np.loadtxt(settings['ring_path'], unpack = True)
    ring = griddata(ring_x, ring_y, model_grid)
    self.print_output('Ring spectrum imported', add_line = False)
    
    self.print_output('Importing gas cross-sections...', add_line = False)
    
    
    # Import SO2 data
    so2_xsec = np.loadtxt(settings['so2_path'], skiprows=1)
    so2_xsec = griddata(so2_xsec[:,0], so2_xsec[:,1], model_grid)
    self.print_output('SO2 cross-section imported', add_line = False)
    
    
    # Import NO2 data
    no2_xsec = np.loadtxt(settings['no2_path'], skiprows=43)
    no2_xsec = griddata(no2_xsec[:,0], no2_xsec[:,2], model_grid)
    self.print_output('NO2 cross-section imported', add_line = False)
    
    
    # Import O3 data
    o3_xsec = np.loadtxt(settings['o3_path'], skiprows=41)
    o3_xsec = griddata(o3_xsec[:,0], o3_xsec[:,1], model_grid)
    self.print_output('O3 cross-section imported', add_line = False)
    
    
    # Import BrO data
    bro = np.loadtxt(settings['bro_path'], skiprows=12)
    
    # As BrO is in terms of wavenumber(cm^-1), need to convert to nm
    bro[:,0] = 10000000/bro[:,0]
    
    # Reverse array so wavelength is accending
    bro = bro[::-1]
    
    # Interpolate onto the grid
    bro_xsec = griddata(bro[:,0], bro[:,1], model_grid)
    self.print_output('BrO cross-section imported')
    
    
    # Add the data to the common dictionary
    common['model_grid']  = model_grid
    common['flat']        = flat
    common['sol']         = sol
    common['ring']        = ring
    common['so2_xsec']    = so2_xsec
    common['no2_xsec']    = no2_xsec
    common['o3_xsec']     = o3_xsec
    common['bro_xsec']    = bro_xsec
    
    return common
