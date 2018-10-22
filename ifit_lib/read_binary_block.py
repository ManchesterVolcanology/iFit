# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 16:43:49 2018

@author: mqbpwbe2
"""

import numpy as np

#========================================================================================
#================================= read_binary_block ====================================
#======================================================================================== 

def read_binary_block(fpath):
    
    '''
    Function to read spectra and header info for a data block created by SpectraLan
    
    INPUTS
    ------
    fpath, str
        File path to data block
    
    OUTPUTS
    -------
    error, bool
        An error code, 0 if all is OK, 1 if an error was produced
        
    wavelength, array
        Wavelength grid calculated from the calibration coefficients
        
    header, str
        Header parameters: Vbatt, Vpanel, IBatt, Temp.
        
    info_block, array
        Spectra info: spec no, hours, minutes, seconds, motor position
        
    spec_block, array
        Array of the measured spectra for the scan block
    
    Written by Ben Esse, June 2018
    '''
    
    try:
        # Get spec details
        scanner, spec_name, int_val, c1, c2, c3 = get_spec_details(fpath)
        
        # Calculate the pixel grid
        pixel_no = np.add(np.arange(2046), 1.0)
        
        # Calculate the wavelength grid
        if c3 == 0.0:
            wavelength = int_val + np.multiply(pixel_no, c1)    + \
                         np.multiply(np.power(pixel_no, 2), c2) 
            
        else:
            wavelength = int_val + np.multiply(pixel_no, c1)    + \
                         np.multiply(np.power(pixel_no, 2), c2) + \
                         np.multiply(np.power(pixel_no, 3), c3)
                         
        # Get spec details
        with open(fpath, 'rb') as rb:
            data = rb.read()
        
        # Strucure of file is as follows:
        # First have header as a string, separated by CR's:
        # Vbatt. Vpanel IBatt. Temp.
        # 12.9268 12.9116 -0.0323 29.7
        # nN_Acq Hour Min Sec MotorPos CoAdding Int_Time Ch_Num Scan_Numb 
        #    Scan_MemInt_Time_Count Pixel_Mode Pixel_Mode_Param ---> Spectral_Data...
        
        # Then spectra, each preceded by an info string
        # 00000 06 56 29 00000
        
        # Split header row using CR's (ASCII code 10)
        CR_ind = []
        for n, i in enumerate(data[:200]):
            if i == 10:
                CR_ind.append(n)
                
        # Unpack indices
        head_idx0 = CR_ind[0] + 1
        head_idx1 = CR_ind[1] - 1
        start_data_idx = CR_ind[2] + 1
        
        # Extract the Vbatt. Vpanel IBatt. Temp.
        header = data[head_idx0:head_idx1].decode('utf-8').split(' ')
        
        # Read in all scans
        all_scans = data[start_data_idx:]
    
        # Calc how many spectra based on file size
        n_spectra = int(len(all_scans)/4131.0)
        
        # Define arrays
        info_block = np.ndarray([5,n_spectra])
        spec_block = np.ndarray([n_spectra, 2046])
        
        # Loop through each scan and add data to array
        # Header is always 37 bytes and the spectrum is 4094 bytes consisting of 2047 two
        #   byte small-endian unsigned integers, miss out the last number
        
        for n in range(n_spectra):
            
            # Extract the entire scan from the data
    
            # Define indices for start and end
            idx0 = n * 4131 
            idx1 = idx0 + 4130
    
            scan = all_scans[idx0:idx1]
            
            # Separate header from the data
            scan_info = scan[0:36]
    
            # Extract jut string data from the header
            info_block[0,n] = int(scan_info[0:5].decode('utf-8'))     # n_acq
            info_block[1,n] = int(scan_info[6:9].decode('utf-8'))     # hour
            info_block[2,n] = int(scan_info[9:12].decode('utf-8'))    # minute
            info_block[3,n] = int(scan_info[12:15].decode('utf-8'))   # second
            info_block[4,n] = float(scan_info[15:20].decode('utf-8')) # motor pos
            
            # Extract the spectral data and copy into the spec_block
            scan_data = scan[37:]
            spec = np.arange(2046)
            
            for i in range(2046):
                
                spec[i] = int.from_bytes([scan_data[i*2], scan_data[i*2 + 1]], 
                                         byteorder='big')
            
            spec_block[n] = spec
        
        
        return 0, wavelength, header, info_block, spec_block
    
    except:
        return 1, 1, 1, 1, 1
    
#========================================================================================
#================================== get_spec_details ====================================
#========================================================================================    
    
def get_spec_details(fpath):
    
    '''
    Function to get spectrometer details given the FLAME network station name
    
    INPUTS
    ------
    fpath, str
        Filename of the spectra block, containing the station name
    
    OUTPUTS
    -------
    scanner, str
        Scanner station name
        
    spec_name, str
        Spectrometer serial number
        
    intercept, float
        Intercept
        
    c1, c2, c3, floats
        Calibration coeficents
    
    Written by Ben Esse, June 2018
    '''
    
    # Define staion parameters
    #         Station  Spec Name  Intercept    C1           C2            C3
    params = {'ecur': ['I2J5773', 296.388306,  0.048254,    -5.345750e-6, -1.687630e-11],
              'enic': ['I2J5769', 296.6937822, 0.047945548, -4.89043e-6,  -1.77072e-10 ],
              'eili': ['I2J5774', 296.2133877, 0.048868629, -5.55088e-6,  3.97945E-11  ],
              'emil': ['I2J5768', 295.9804618, 0.049231176, -5.52944e-6,  -6.98636E-12 ],
              'even': ['I2J5775', 296.2851694, 0.04864757,  -5.17264e-6,  -0.106506e-10],
              'etst': ['I2J5770', 295.1845975, 0.049603817, -5.54717e-6,  -3.531373-11 ]}
    
    # Get just filename
    fname = fpath.split('/')[-1]
    
    # Split fname into info pieces and extract station name
    scanner = fname.split('_')[2]
    
    # Get spectrometer parameters from the station name
    spec_name, intercept, c1, c2, c3 = params[scanner]
    
    return scanner, spec_name, intercept, c1, c2, c3