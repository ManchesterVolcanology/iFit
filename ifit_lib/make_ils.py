# -*- coding: utf-8 -*-
"""
Created on Wed May 10 14:16:04 2017

@author: mqbpwbe2
"""

import numpy as np
from scipy.signal import gaussian
from scipy.special import gamma
    
#========================================================================================
#======================================== make_ils ======================================
#========================================================================================  

def make_ils(interval, FWEM, k = 2, a_w = 0, a_k = 0):
    
    '''
    Function to generate a synthetic instrument line shape based on a super gaussian
    function:
        
                     { exp(-| x / (w-a_w) | ^ (k-a_k)) for x <= 0
    G(x) = A(w, k) * {
                     { exp(-| x / (w+a_w) | ^ (k+a_k)) for x > 0
    
    where A(w, k) = k / (2 * w * Gamma(1/k)).
    
    See Beirle et al (2017) for more details: doi:10.5194/amt-10-581-2017
    
    INPUTS
    ------
    interval, int
        The spacing of the wavelength grid on which the ILS is built
        
    FWEM, float
        The Full Width eth Maximum of the lineshape, defined as 2*w = FWEM
        
    k, float
        Controls the shape of the lineshape (default = 2): 
            - k < 2 -> sharp point and wide tails
            - k = 2 -> normal Gaussian
            - k > 2 -> flat top, approaches boxcar at k -> inf
            
    a_w and a_k, float,
        Controls the asymetry of the lineshape
        
    OUTPUTS
    -------
    ils, numpy array
        The calculated ILS function on a wavelength grid of the given spacing and 5 times
        the width of the supplied FWEM
    '''

    # Create a grid 6 times that of the width
    grid = np.arange(-FWEM * 2, FWEM * 2, interval)
    
    # Calculate w as half of the FWEM
    w = 0.5 * FWEM
    
    # Form empty array
    super_g = np.zeros(len(grid))
    
    # Calculate A
    A = k / (FWEM * gamma(1/k))
    
    # Split the x grid into =ve and -ve arrays
    neg_idx = np.where(grid <= 0)
    pos_idx = np.where(grid > 0)
    neg_grid = grid[neg_idx]
    pos_grid = grid[pos_idx]
    
    # Calculate the asymetric supergaussian function
    neg_g = np.multiply(A, np.exp(-np.power(np.abs((neg_grid) / (w - a_w) ), k - a_k)))
    pos_g = np.multiply(A, np.exp(-np.power(np.abs((pos_grid) / (w + a_w) ), k + a_k)))
    
    # Combine
    super_g = np.append(neg_g, pos_g)
    
    super_g = np.divide(super_g, sum(super_g))
    
    return super_g  

#========================================================================================
#======================================== make_ils ======================================
#========================================================================================
 
def make_gauss(fwhm, interval):
    
    '''
    Gaussian function
    
    INPUTS
    ------
    interval, float
        Model grid spacing
        
    fwhm, float
        The full width half maximum of the desired Gaussain
    
    OUTPUTS
    -------
    gauss, array
        Outputted normalised gaussian lineshape
    '''
    
    # Create grid 5 times bigger than FWHM
    n = int((fwhm * 5) / interval)
    
    # Define sigma
    sigma = fwhm / 2.355
    
    # Create gauss
    gauss = gaussian(n, sigma / interval)
    
    # Normalise
    gauss = np.divide(gauss, np.sum(gauss))
    
    return gauss
    
#========================================================================================
#======================================== make_ils ======================================
#========================================================================================    

def make_ils_old(fit_res, interval, ils_gauss_weight):
    
    '''
    Function to calculate the intrument line shape (ILS) of the spectrometer
    
    INPUTS
    ------
    fit_res, float
        Resolution of the fit (FWHM for a pure gaussian ILS)
        
    interval, float
        Model grid spacing
        
    ils_gaus_weight, float (0 - 1)
        Weighting of gaussian contributio to ILS. Boxcar contribution is caclulated as 
        1 - ils_gauss_wieght
    
    OUTPUTS
    -------
    ils, array
        The instrument line shape (ILS) of the spectrometer
    '''
    
    # Calculate pure Gaussian lineshape
    gauss_ils = make_gauss(fit_res, interval)
    
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

