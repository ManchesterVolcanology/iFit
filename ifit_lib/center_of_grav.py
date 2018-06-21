import numpy as np
from ifit_lib.find_nearest import find_nearest

#========================================================================================
#=====================================centre_of_grav=====================================
#========================================================================================

def cog(array):
    
    '''
    Define function to find the centre of gravity of an array by finding the index of the 
    value equal to half the total cumulative frequency

    INPUTS
    ------
    array: array of values to find the centre of gravity of
    
    OUTPUTS
    -------
    index: index of the centre of gravity in the input array
        
    '''
    
    cul_freq = np.zeros(len(array))
    
    # Find cumulative frequency of the array
    for n, val in enumerate(array):
        if n == 0:
            cul_freq[n] = val
        else:
            cul_freq[n] = cul_freq[n-1] + val
    
    # Calulate half value
    half_total = cul_freq[-1] / 2.0
    
    # Find nearest point to the given value
    index, val = find_nearest(cul_freq, half_total)
    
    return index

