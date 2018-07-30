# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 09:13:44 2018

@author: mqbpwbe2
"""

import numpy as np
from scipy.interpolate import griddata

#========================================================================================
#=======================================Get Shift========================================
#========================================================================================

def get_shift(grid, spec, common):
    
    '''
    
    Function to pre-calculate the spectrometer shift using a cross-correlation of the 
    solar spectrum with the measured spectrum
    
    INPUTS
    ------
    grid:   Measurement wavelength grid
    spec:   Measured spectrum
    common: Common dictionary of parameters form the main program
    
    OUTPUTS
    -------
    shift: Calculated wavelength shift
    
    '''
    
    # Interpolate measurment onto high resolution model grid
    spec = griddata(grid, spec, common['model_grid'], method = 'cubic')
    
    # Center spectra on 0
    norm_spec = np.subtract(spec, np.average(spec))
    norm_sol = np.subtract(common['sol'], np.average(common['sol']))
    
    # Perform cross coreelation and calc the shift
    corr = np.correlate(norm_spec, norm_sol, 'full')
    shift = (corr.argmax() - (len(corr)/2))/100
    
    return shift
    
    
    
    
    
    
    
    
    
    
    