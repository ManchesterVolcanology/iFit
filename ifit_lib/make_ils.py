# -*- coding: utf-8 -*-
"""
Created on Wed May 10 14:16:04 2017

@author: mqbpwbe2
"""

import numpy as np
from scipy import signal

from ifit_lib.find_nearest import extract_window
#from find_nearest import extract_window

#========================================================================================
#=========================================make_ils=======================================
#========================================================================================

# Function to produce an instrument line shape for a spectrometer. Function is a gaussian 
#  defined by the equation:
#          
#                        f(x) = a*exp(-(x-b)^2 / (2c^2))
#
#  where x is the wavelength, a is the amplitude, b is the position, and c is the 
#  standard deviation. For this function the position is not important so b = 0, and a is
#  defined such that the total area is equal to 1

# INPUTS: grid;    wavelength array over which to make the ils
#         std_dev; standard deviation of the gaussian shape

# OUTPUTS: g; gaussian instrument line shape
 
def make_gauss(grid, stdev):

    # Make gaussian shape
    g = signal.gaussian(len(grid), stdev/(grid[1] - grid[0]))
    
    # Normalise such that the total area is one
    g = np.divide(g, np.sum(g))
    return g
    
#========================================================================================
#=========================================calc_ils=======================================
#========================================================================================    

def make_ils(grid, fit_res, ils_gauss_weight):
    '''
    # Avoid problems with silly resolutions
    if fit_res < 0.1 or fit_res > 2.0:
        fit_res = 0.3
        
    if np.isfinite(fit_res) == False:
        fit_res = 0.3
    ''' 
    # Calculate pure gaussian lineshape
    gauss = make_gauss(grid, fit_res)
    
    # Calculate pure boxcar lineshape
    box = np.zeros(len(grid))
    midpoint = grid[int(len(grid) / 2)]
    half_width = fit_res / 2
    window, ind1, ind2 = extract_window(grid, midpoint-half_width, midpoint+half_width)
    box[ind1:ind2] = 1.0 / len(window)
    
    # Calculate boxcar weighting
    ils_box_weight = 1.0 - ils_gauss_weight

    # Apply weightings
    gauss = np.multiply(gauss, ils_gauss_weight)
    box = np.multiply(box, ils_box_weight)
    
    # Combine
    ils = np.add(gauss, box)
    
    # Normalise
    ils = np.divide(ils, np.sum(ils))
    
    return ils



def make_gauss_ils(grid, fit_res):
    
    # Calculate pure gaussian lineshape
    gauss = make_gauss(grid, fit_res)
    
    # Normalise
    ils = np.divide(gauss, np.sum(gauss))
    
    return ils

'''
import matplotlib.pyplot as plt
grid = np.linspace(308,318, num = 10000)

fit_res = 0.22

gauss_wt = 0.5

ils1 = make_gauss(grid,fit_res)
ils2 = calc_ils(grid, fit_res, gauss_wt)

plt.plot(grid, ils1, grid, ils2)
plt.show()
'''
