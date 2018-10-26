# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 10:21:57 2018

@author: mqbpwbe2
"""

import numpy as np

def read_csv(fpath, delim = ',', skip_rows = 0):
    
    '''
    Function to read in a .csv file and return the results as a dictionary.
    
    INPUTS
    ------
    fpath, str
        File path to the file to read in
        
    delim, str (optional)
        Delimiter used in the file. Default is ','
        
    skip_rows, int (optional)
        Number of rows to skip from the top of the document. Default is 0
        
    OUTPUTS
    -------
    data, dict
        Dictionary of the imported data, using the column headings as keys
    '''
    
    # Open file
    with open(fpath, 'r') as r:
        
        # Read in the data
        lines = r.readlines()
        
        # First line is the header. Unpack to get the keys
        keys = lines[skip_rows].strip().split(delim)

        # Unpack the data into an array
        raw_data = [line.strip().split(delim) for line in lines[1:]]
        
        # Initialise the results dictionary
        data = {}
        for k in keys:
            data[k] = [''] * len(raw_data)
        
        for n in range(len(raw_data)):
            
            # Add the data point in each line to the correspoding dictionary array
            for i, k in enumerate(keys):
                
                data[k][n] = raw_data[n][i]
                
        # Try to get the correct data format
        for k in keys:
            
            try:
                data[k] = np.array(data[k], dtype = float)
                
            except ValueError:
                pass
        
        
    return data
        
        
        
        
        
#data = read_csv('C:/Users/mqbpwbe2/Dropbox/python_scripts/iFit/Results/iFit/2018-01-14/iFit_out.csv')