import numpy as np

#========================================================================================
#======================================find_nearest======================================
#========================================================================================

# Function to find the nearest value in an array given a desired value

# INPUTS: array: array to search through
#         value: desired value to find in the array

# OUTPUTS: idx: index of the nearest array value
#          array[idx]: the true array value

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]

#========================================================================================
#=====================================extract window=====================================
#========================================================================================

# Function to extract a window from an ordered array

# INPUTS: array      ; array to window
#         start, stop; lower and upper values defining the window

# OUTPUTS: windowed_array; cut version of the input array
#          ind1/ind2     ; indicies of the window

def extract_window(array, start, stop):
    
    # Find the nearest values in the array to the start and stop values
    ind1 = (np.abs(array - start)).argmin()
    ind2 = (np.abs(array - stop)).argmin()
    
    # Cut the array
    windowed_array = array[ind1:ind2]
    
    # Return the windowed array
    return windowed_array, ind1, ind2