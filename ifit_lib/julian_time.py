# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 09:36:21 2018

@author: mqbpwbe2
"""

import numpy as np
from datetime import datetime as dt

def hms_to_julian(time_arr, str_format = None, out_format = 'decimal days'):
    
    '''
    Function to convert an array of datetime objects or strings to julian time
    
    INPUTS
    ------
    time_arr  : array; origional time objects
    str_format: string; time format if converting from strings. Default is None
    out_format: string; Format of output:
                   
    OUTPUTS
    -------
    jul_time: Array of julian time values calculated form the input array
    '''
    
    # Check if output format is correct
    if out_format not in ['decimal days', 'decimal hours']:
        msg = 'out_format value is not supported. Only supported values are:\n' + \
              '    decimal days\n    decimal hours'
        raise ValueError(msg)
    
    # Make julian time array object to populate
    jul_time = np.zeros(len(time_arr))
    
    # Iterate through time array and convert
    for n, time in enumerate(time_arr):
        
        # If format is given, use this to convert to datetime object
        if str_format != None:
            
            time = dt.strptime(str(time), format)
        
        # Extract individual hour, minute, second and microsecond info
        hours = time.hour
        mins  = time.minute
        secs  = time.second
        usecs = time.microsecond
        
        # Convert to either decimal days or decimal hours
        if out_format == 'decimal days':
            jul_time[n] = (hours/24) + (mins/1440) + (secs/86400) + (usecs/8.64e10)
        if out_format == 'decimal hours':
            jul_time[n] = (hours) + (mins/60) + (secs/3600) + (usecs/3.6e+9)
        
    return jul_time
            
        
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            