import numpy as np

#====================================================================================
#===================================== smooth =======================================
#====================================================================================

# Define function smooth: smooths an array using a boxcar of a defined width

# NOTE! This function will have edge effects! Apply smooth BEFORE cutting window

# INPUTS: y; array to be smoothed

#         window; width of boxcar to smooth over  

def smooth(y, width, function = 'boxcar'):
    
    '''
    Function to smooth a 1D array using a boxcar of defined width
    
    NOTE! This function will have edge effects! Apply smooth BEFORE cutting window
    
    INPUTS
    ------
    y, array
        The array to smooth
        
    width, int
        The width of the boxcar (in array elements)
        
    function, str (optional)
        Defines the function used to smooth. Default is boxcar.
        (Currently only boxcar is suported)
        
    OUTPUTS
    -------
    smooth_y, array
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