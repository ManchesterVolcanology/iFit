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



def make_directory(fpath):
    
    '''
    Function used to build a new file directory, automatically changing the filepath to
    avoid overwriting existing data
    
    INPUTS
    ------
    fpath: desired file path to the directory
    
    OUTPUTS
    -------
    new_fpath: actual filepath to the new directory
    '''
    
    # To avoid overwriting data, create a new directory if one already exists
    new_fpath = fpath
    
    if not os.path.exists(fpath):
        # Make the directory
        os.makedirs(fpath)

    else:
        i = 1
        while os.path.exists(new_fpath):
            if fpath[-2] == ')':
                new_fpath = fpath[:-4] + '(' +str(i) + ')/'
                
            else:
                new_fpath = fpath[:-1] + '(' +str(i) + ')/'
            
            i += 1
            
        os.makedirs(new_fpath)
        
    return new_fpath