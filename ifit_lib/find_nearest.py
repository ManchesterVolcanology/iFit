import numpy as np

#========================================================================================
#===================================== find_nearest =====================================
#========================================================================================

def find_nearest(array, value):
    
    '''
    Function to find nearest value in an array
    
    INPUTS
    ------
    array, array of floats/integers
        Array to search
        
    value, float
        Value to find in array
    
    OUTPUTS
    -------
    idx, int
        Index of nearest value
        
    val, float
        Value of nearest value in array
    '''
    
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]

#========================================================================================
#==================================== extract window ====================================
#========================================================================================
 
# Function to extract a window from an ordered array
 
# INPUTS: array      ; array to window
#         start, stop; lower and upper values defining the window
 
# OUTPUTS: windowed_array; cut version of the input array
#          ind1/ind2     ; indicies of the window
 
def extract_window(array, start, stop, ordered = True):
    
    '''
    Function to extract a window from an array
    
    INPUTS
    ------
    array, array
        The array to extract from.
        
    start, float
        Lower window bound
    
    stop, float
        Upper window bound
        
    ordered, bool (optional)
        Flag to use if the input array is ordered or not (default = True)
        
    OUTPUT
    ------
    windowed_array, array
        The origional array cut to the desired window
        
    cond, tuple
        The upper and lower bound indices
    '''
    
    if ordered:
    
        # Find the nearest values in the array to the start and stop values
        idx0 = (np.abs(array - start)).argmin()
        idx1 = (np.abs(array - stop)).argmin()
         
        # Cut the array
        windowed_array = array[idx0:idx1]
        
        # Return the windowed array
        return windowed_array, idx0, idx1
        
    else:
        cond = np.zeros(len(array))
        
        for n, i in enumerate(array):
            if i > start and i < stop:
                cond[n] = 1
                
        windowed_array = np.extract(cond, array)
                
        return windowed_array, cond