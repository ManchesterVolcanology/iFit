import numpy as np
from scipy import signal

#====================================================================================
#=====================================smoother=======================================
#====================================================================================

# Define function smooth: smooths an array using a boxcar of a defined width

# NOTE! This function will have edge effects! Apply smooth BEFORE cutting window

# INPUTS: y; array to be smoothed

#         window; width of boxcar to smooth over  

def smooth(y,width):

    # Create the boxcar window
    window = np.ones(int(width))/float(width)
    
    # Pad the array with values to avoid edge effects
    pre_array = np.ones(width-1) * y[0]
    post_array = np.ones(width-1) * y[-1]
    
    # Add padding to the origional array
    y = np.append(pre_array,y)
    y = np.append(y,post_array)
    
    # Convolve with boxcar to smooth
    y = np.convolve(y,window,'same')
    
    # Cut array to origional size    
    y = y[width-1:-(width-1)]
    
    return y
    
#====================================================================================
#====================================smooth ils======================================
#====================================================================================

# Define function smooth: smooths an array using a gaussian ils of a defined width

# NOTE! This function will have edge effects! Apply smooth BEFORE cutting window

# INPUTS: y; array to be smoothed

#         window; standard deviation of gaussian 

def smooth_ils(grid,y,width):

    # Create the boxcar window
    window = signal.gaussian(len(y),width/(grid[1]-grid[0]))
    
    # Pad the array with values to avoid edge effects
    pre_array = np.zeros(width-1) * y[0]
    
    post_array = np.zeros(width-1) * y[-1]
    
    # Add padding to the origional array
    #y = np.append(pre_array,y)
   # y = np.append(y,post_array)
    
    # Convolve with boxcar to smooth
    y = np.convolve(y,window,'same')
    
    # Cut array to origional size    
   # y = y[width-1:-(width-1)]
    
    return y













