# -*- coding: utf-8 -*-
"""
Created on Wed May 10 14:16:04 2017

@author: mqbpwbe2
"""

import numpy as np

#========================================================================================
#=========================================make_ils=======================================
#========================================================================================
 
def make_gauss(interval, fwhm):
    
    '''
    Gaussian function
    
    INPUTS
    ------
    interval: Model grid spacing
    fwhm: The full width half maximum of the Gaussain curve
    
    OUTPUTS
    -------
    gauss: Outputted normalised gaussian lineshape
    '''
    
    # Create grid 5 times bigger than FWHM
    n = (fwhm * 5) / interval
    x = np.divide(np.arange(n), (n-1))
    x = np.subtract(x, 0.5)
    x = np.multiply(x, (fwhm*5))
    
    # Define sigma
    sigma = fwhm / 2.36
    
    # Create gauss
    gauss = (1/(sigma*np.sqrt(2*np.pi))) * np.exp(-0.5*np.power((np.divide(x, sigma)),2))
    
    return gauss
    
#========================================================================================
#=========================================make_ils=======================================
#========================================================================================    

def make_ils(fit_res, interval, ils_gauss_weight):
    
    '''
    Function to calculate the intrument line shape (ILS) of the spectrometer
    
    INPUTS
    ------
    fit_res: Resolution of the fit (FWHM for a pure gaussian ILS)
    interval: Model grid spacing
    ils_gaus_weight: Weighting of gaussian contributio to ILS. Float from 0-1. Boxcar 
                     contribution is caclulated as 1 - ils_gauss_wieght
    
    OUTPUTS
    -------
    ils: the instrument line shape (ILS) of the spectrometer
    '''
    
    # Calculate pure Gaussian lineshape
    gauss_ils = make_gauss(interval, fit_res)
    
    # Calculate pure boxcar lineshape and normalise
    box_npts = int(fit_res/interval)+1
    box_ils = np.divide(np.ones(box_npts), box_npts)
    
    # Apply weightings
    gauss_ils = np.multiply(gauss_ils, ils_gauss_weight)
    box_ils = np.multiply(box_ils, (1 - ils_gauss_weight))
    
    # Combine
    idx0 = int(((len(gauss_ils)/2)+1) - (box_npts/2))
    idx1 = int((len(gauss_ils)/2) + (box_npts/2))
    idx0 = idx0 + ((idx1-idx0) - (box_npts))
    gauss_ils[idx0:idx1] = np.add(gauss_ils[idx0:idx1], box_ils)
    ils = gauss_ils
    
    # Normalise
    ils = np.divide(ils, np.sum(ils))
    
    return ils

