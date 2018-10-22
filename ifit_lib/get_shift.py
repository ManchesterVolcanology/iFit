# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 09:13:44 2018

@author: mqbpwbe2
"""

import numpy as np
from scipy.interpolate import griddata

from ifit_lib.make_ils import make_ils

#========================================================================================
#====================================== get_shift =======================================
#========================================================================================

def get_shift(grid, spec, common):
    
    '''
    
    Function to pre-calculate the spectrometer shift using a cross-correlation of the 
    solar spectrum with the measured spectrum
    
    INPUTS
    ------
    grid, array
        Wavelength grid (must be evenly spaced)
        
    spec, array
        Measured spectrum
        
    common, dict
        Common dictionary of parameters from the main program
    
    OUTPUTS
    -------
    shift, float
        Calculated wavelength shift
    
    err, bool
        Flag showing sensibility of the fit. If 1 there is a problem.
    '''
    
    # Find spectral spacing
    spacing = common['model_grid'][1] - common['model_grid'][0]

    # Interpolate measurment onto high resolution model grid
    spec = griddata(grid, spec, common['model_grid'], method = 'cubic')

    # Convolve high resolution spectrum with ils
    ils = make_ils(common['ils_width'], spacing, common['gauss_weight'])
    smooth_sol = np.convolve(common['sol'], ils, 'same')
    
    # Center spectra on 0
    norm_spec = np.subtract(spec, np.average(spec))
    norm_sol = np.subtract(smooth_sol, np.average(smooth_sol))

    # Perform cross coreelation and calc the shift
    corr = np.correlate(norm_spec, norm_sol, 'full')
    shift = (corr.argmax() - (len(corr)/2)) * spacing
    
    # Check sensibilty of calculated shift
    if np.abs(shift) > 3.0:
        err = 1
    else:
        err = 0
    
    return shift, err