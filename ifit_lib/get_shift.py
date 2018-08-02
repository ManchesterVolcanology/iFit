# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 09:13:44 2018

@author: mqbpwbe2
"""

import numpy as np
from scipy.interpolate import griddata

from ifit_lib.make_ils import make_ils

#========================================================================================
#=======================================Get Shift========================================
#========================================================================================

def get_shift(grid, spec, common):
    
    '''
    
    Function to pre-calculate the spectrometer shift using a cross-correlation of the 
    solar spectrum with the measured spectrum
    
    INPUTS
    ------
    grid:   Wavelength grid (must be evenly spaced)
    spec:   Spectrum
    common: Common dictionary of parameters form the main program
    
    OUTPUTS
    -------
    shift: Calculated wavelength shift
    
    '''
    
    # Find spectral spacing
    spacing = common['model_grid'][1] - common['model_grid'][0]

    # Interpolate measurment onto high resolution model grid
    spec = griddata(grid, spec, common['model_grid'], method = 'cubic')

    # Convolve high resolution spectrum with ils
    ils = make_ils(common['ils_width'], spacing, common['ils_gauss_weight'])
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



#########################################################################################
'''    
common = {}

common['model_grid'] = np.arange(start = 302, stop = 321 + 0.02, step = 0.02)
common['ils_width'] = 0.51
sol_path = 'C:/Users/mqbpwbe2/Dropbox/python_scripts/iFit/data_bases/gas data/sao2010.txt'
sol_x, sol_y = np.loadtxt(sol_path, unpack = True)    
conv_factor = sol_y.max() / 70000
sol_y = np.divide(sol_y, conv_factor)
common['sol'] = griddata(sol_x, sol_y, common['model_grid'], method = 'cubic')
common['ils_gauss_weight'] = 1.0

from ifit_lib.read_spectrum import read_spectrum
spec_path = 'C:/Users/mqbpwbe2/Dropbox/Fieldwork/2018 01 Central America/14-01-2018/ifrit_out/spectra/spectrum_00350.txt'
x, y, date, time, n, err = read_spectrum(spec_path)

shift = get_shift(x, y, common)
'''

    
    
    
    
#########################################################################################    
    
    