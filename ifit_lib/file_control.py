# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 15:40:42 2017

@author: mqbpwbe2
"""

import os
from glob import glob

#========================================================================================
#=====================================make_directory=====================================
#========================================================================================

def make_directory(fpath, overwrite = False):
    
    '''
    Function used to build a new file directory, automatically changing the filepath to
    avoid overwriting existing data
    
    INPUTS
    ------
    fpath: desired file path to the directory
    overwrite: Boolian flag controlling if an existing directory is overwritten
    
    OUTPUTS
    -------
    new_fpath: actual filepath to the new directory
    '''
    
    # To avoid overwriting data, create a new directory if one already exists
    new_fpath = fpath

    if not os.path.exists(fpath):
        # Make the directory
        os.makedirs(fpath)
        return fpath 
    
    else:
        
        # Create counter
        i = 1
        
        # Cycle through possible filenames
        while True:
            
            # If it doesn't exist, create the folder
            if not os.path.exists(new_fpath):
                os.makedirs(new_fpath)
                break
            
            # If it exists and overwrite is on, check if it is empty (appart from notes)
            elif overwrite == True and len(glob(new_fpath + '*')) < 2:
                break
                
            # If it exists and isn't empty then create new folder name
            new_fpath = fpath[:-1] + '(' +str(i) + ')/'

            i += 1
    
    return new_fpath

#========================================================================================
#=====================================make_directory=====================================
#========================================================================================

def make_csv_file(fname, common):
    
    '''
    Function used to build a new file directory, automatically changing the filepath to
    avoid overwriting existing data
    
    INPUTS
    ------
    fpath: desired file path to the directory
    '''
    
    # Open excel file and write header line
    with open(fname, 'w') as writer:
        
        # Write header line
        #writer.write('File,Number,Date,Time,so2 (ppm.m),so2 error,')
        writer.write('File,Number,Date,Time,')
        
        for i in common['params'].keys():
            writer.write(i + ',' + i + '_e,')
            
        # Add column for fit info
        writer.write('Fit Quality\n')
            
        # Write other fit info
        #writer.write('Fit window: ' + str(common['wave_start']) + ' - ' + \
        #             str(common['wave_stop']) + ' nm\n')