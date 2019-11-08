# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 10:32:32 2019

@author: mqbpwbe2
"""

import numpy as np
from scipy.special import gamma

#==============================================================================
#================================== make_ils ==================================
#==============================================================================

def make_ils(interval, FWEM, k = 2, a_w = 0, a_k = 0):

    '''
    Function to generate a synthetic instrument line shape based on a super
    gaussian function:

                     { exp(-| x / (w-a_w) | ^ (k-a_k)) for x <= 0
    G(x) = A(w, k) * {
                     { exp(-| x / (w+a_w) | ^ (k+a_k)) for x > 0

    where A(w, k) = k / (2 * w * Gamma(1/k)).

    See Beirle et al (2017) for more details: doi:10.5194/amt-10-581-2017

    **Parameters**

    interval : int
        The spacing of the wavelength grid on which the ILS is built

    FWEM : float
        The Full Width eth Maximum of the lineshape, defined as 2*w = FWEM

    k : float, optional, default = 2
        Controls the shape of the lineshape (default = 2):
            - k < 2 -> sharp point and wide tails
            - k = 2 -> normal Gaussian
            - k > 2 -> flat top, approaches boxcar at k -> inf

    a_w and a_k : float, optional, default = 0
        Controls the asymetry of the lineshape

    **Returns**

    ils, numpy array
        The calculated ILS function on a wavelength grid of the given spacing
        and 5 times the width of the supplied FWEM
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
    neg_g = np.multiply(A, np.exp(-np.power(np.abs((neg_grid) / (w - a_w) ),
                                            k - a_k)))
    pos_g = np.multiply(A, np.exp(-np.power(np.abs((pos_grid) / (w + a_w) ),
                                            k + a_k)))

    # Combine
    super_g = np.append(neg_g, pos_g)

    ils = np.divide(super_g, sum(super_g))

    return ils

#==============================================================================
#=================================== smooth ===================================
#==============================================================================

def smooth(y, width):

    '''
    Function to smooth a 1D array using a boxcar of defined width This function
    will have edge effects! Apply smooth BEFORE cutting window

    **Parameters**

    y : array
        The array to smooth

    width : int
        The width of the boxcar (in array elements)

    **Returns**

    smooth_y : array
        The smoothed array
    '''

    # Create the boxcar window
    window = np.ones(int(width))/float(width)

    # Pad the array with values to avoid edge effects
    pre_array = np.ones(width-1) * y[0]
    post_array = np.ones(width-1) * y[-1]

    # Add padding to the origional array
    y = np.append(pre_array,y)
    y = np.append(y,post_array)

    # Convolve with boxcar to smooth
    smooth_y = np.convolve(y, window, 'same')

    # Cut array to origional size
    smooth_y = smooth_y[width-1:-(width-1)]

    return smooth_y