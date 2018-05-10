# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 14:30:41 2017

@author: mqbpwbe2
"""

import linecache
import numpy as np

#========================================================================================
#======================================read_spectrum=====================================
#========================================================================================

# Reads spectrum file from spectro_gui

def read_spectrum(fname, spec_type):

    if spec_type == 'IFRiT':
        
        # Load data into a numpy array
        x, y = np.loadtxt(fname, unpack = True, skiprows = 8)
        
        # Extract date and time string
        date_time = linecache.getline(fname, 5)[27:47]
        
        # Unpack date_time string
        date = date_time[0:10]
        time = date_time[11:19]
        
        # Get spectrum number
        spec_no = int(fname[-9:-4])
        
    if spec_type == 'Master.Scope':
        
        # Load data into a numpy array, skipping the header data
        x,y = np.genfromtxt(fname, unpack=True, skip_header = 19, skip_footer = 1)
        
        # Extract date and time string
        date_time = linecache.getline(fname, 3)[6:]
        date, time = [x.strip() for x in date_time.split(',')]
        
        # Get spectrum number
        spec_no = int(fname[-18:-13])
        
    if spec_type == 'Jai Spec':
    
        # Load data into a numpy array, skipping the header data
        x,y = np.loadtxt(fname, unpack=True, skiprows = 13, delimiter = ';')
        
        # Extract date and time string
        data_date_time = linecache.getline(fname, 2)
        
        # Extract day, month, year and time from the string
        date = data_date_time[0:10]
        time = data_date_time[11:19]
        
        # Get spectrum number
        spec_no = int(fname[-10:-4])
        
    if spec_type == 'Spectrasuite':
        
        # Load data into a numpy array, skipping the header data
        x,y = np.genfromtxt(fname, unpack=True, skip_header = 17, skip_footer = 2)
        
        # Extract date and time string
        line = linecache.getline(fname, 3)
        
        # Extract day, month, year and time from the string
        day = line[14:16]
        month = line[10:13]
        year = line[30:34]
        time = line[17:25]
        
        # Convert month to numerical format
        months = {'Jan': '01',
                  'Feb': '02',
                  'Mar': '03',
                  'Apr': '04',
                  'May': '05',
                  'Jun': '06',
                  'Jul': '07',
                  'Aug': '08',
                  'Sep': '09',
                  'Oct': '10',
                  'Nov': '11',
                  'Dec': '12'}
                  
        month = months[month]
        
        # Reform numerical date and time and return
        date = year + '-' + month + '-' + day
        
        # Get spectrum number
        spec_no = int(fname[-9:-4])
    
    return x, y, date, time, spec_no
  
#========================================================================================
#=====================================average_spectra====================================
#========================================================================================
  
def average_spectra(files, spec_type):
    
    for fname in files:
        
        # Check if file exists, if true then read
        try:
          
            # Load spectrum
            grid, y, spec_date, spec_time, spec_no = read_spectrum(fname, spec_type)
            
            # Sum spectra
            if fname == files[0]:
                spec = y
            else:
                spec += y
         
        except FileNotFoundError:
            spec = spec
            
    # Divide to get average spectrum
    spec = spec / len(files)
    
    return (grid, spec)