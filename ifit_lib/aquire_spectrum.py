# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 13:24:15 2017

@author: mqbpwbe2
"""

import numpy as np
import datetime

#========================================================================================
#======================================read_spectro======================================
#========================================================================================

def aquire_spectrum(spec, integration_time_ms, coadds, dk_flag, nonlin_flag, q = None):
    
    '''
    Subprocess to read an ocean optics spectrometer

    INPUTS: 
    -------
    spec:                the spectrometer object
    integration_time_ms: integration time in milliseconds
    coadds:              number of spectra to average per reading. 
    dk_flag:             boolian flag to controll the correction of the electronic dark
    nonlin_flag:         boolian flag to controll the correction of non-linearity 
    q:                   Queue to which to add the output if threaded (default = None)
    
    OUTPUTS:
    --------
    x:      wavelength array
    y:      intensity array
    header: text string containing the information of the reading, formatted as a header 
              for a text file
    '''
    
    # Read wavelength grid
    x = spec.wavelengths()
    
    # Create empty array to hold intensities
    y = np.zeros(2048)
    
    # Loop over the coadds
    for i in range(coadds):
        y = np.add(y, spec.intensities(correct_dark_counts = dk_flag,
                                       correct_nonlinearity = nonlin_flag))
        
    # Divide by number of coadds
    y = np.divide(y, coadds)
    
    # Take time of reading
    read_time = str(datetime.datetime.now())
    
    #Create header and comments of file with spectrometer settings and 
    # date and time of the reading
    header  = 'Ocean optics spectrum file\n' + \
              'Spectrometer: ' + str(spec.serial_number) +'\n' + \
              'Integration time: ' + str(integration_time_ms) + ' (ms)\n' + \
              'Number of coadds: ' + str(coadds) + '\n' + \
              'Date/Time (end of read): ' + read_time + '\n' + \
              'Electronic dark correction: ' + str(dk_flag) + '\n' + \
              'Non-linearity correction: ' + str(nonlin_flag) + '\n' + \
              'Wavelength (nm),       Intensity (arb)'

    
    # Return the spectrum and file header
    if q == None:
        return x, y, header, read_time
    
    else:
       # Pack the spectrum and file header into the function output
        output = (x, y, header, read_time)
        
        # Add the output to the queue
        q.put(('spectrum', output)) 