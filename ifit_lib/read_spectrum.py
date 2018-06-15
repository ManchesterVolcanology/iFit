# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 14:30:41 2017

@author: mqbpwbe2
"""

import linecache
import numpy as np
import datetime
from pandas import read_csv as read_csv

#========================================================================================
#======================================read_spectrum=====================================
#========================================================================================

# Reads spectrum file from spectro_gui

def read_spectrum(fname, spec_type='iFit'):
    
    try:
        
        if spec_type == 'iFit':
            
            # Load data into a numpy array
            x, y = np.loadtxt(fname, unpack = True, skiprows = 8)
            
            # Extract date and time string
            date_time = linecache.getline(fname, 5)[27:].strip()
            
            # Unpack date_time string
            try:
                date_time = datetime.datetime.strptime(date_time, '%Y-%m-%d %H:%M:%S.%f')
            except ValueError:
                date_time = datetime.datetime.strptime(date_time, '%Y-%m-%d %H:%M:%S')
            
            # Get spectrum number
            spec_no = int(fname[-9:-4])
            
        if spec_type == 'Master.Scope':
            
            # Load data into a numpy array, skipping the header data
            x,y = np.genfromtxt(fname, unpack=True, skip_header = 19, skip_footer = 1)
            
            # Extract date and time string
            date_time = linecache.getline(fname, 3)[6:].strip()        
            
            # Unpack date_time string
            date_time = datetime.datetime.strptime(date_time, '%m-%d-%Y, %H:%M:%S')
            
            # Get spectrum number
            spec_no = int(fname[-18:-13])
            
        if spec_type == 'Jai Spec':
        
            # Load data into a numpy array, skipping the header data
            x,y = np.loadtxt(fname, unpack=True, skiprows = 13, delimiter = ';')
            
            # Extract date and time string
            date_time = linecache.getline(fname, 2).strip() 
            
            # Unpack date_time string
            date_time = datetime.datetime.strptime(date_time, '%Y/%m/%d %H:%M:%S')
            
            # Get spectrum number
            spec_no = int(fname[-10:-4])
            
        if spec_type == 'Spectrasuite':
            
            # Load data into a numpy array, skipping the header data
            x,y = np.genfromtxt(fname, unpack=True, skip_header = 17, skip_footer = 2)
            
            # Extract date and time string
            date_time = linecache.getline(fname, 3).strip() 
            
            # Unpack date_time string
            date_time = datetime.datetime.strptime(date_time, '%a %b %Y %H:%M:%S')
            
            # Get spectrum number
            spec_no = int(fname[-9:-4])
            
        if spec_type == 'GSJ':
            
            # Different because the wavelength and timing is separate to the spectrum

            # Get folder of spectra
            n = [pos for pos, c in enumerate(fname) if c == '/']
            spec_folder = fname[:n[-1] + 1]

            # Read Wavelength File
            x = np.loadtxt(spec_folder + 'wavelength.dat')
            
            # Read spectrum file
            y = np.loadtxt(fname)
            
            # Get spectrum number
            spec_no = int(fname[-8:-4])
            
            # Get timing from so2.dat
            times = read_csv(spec_folder+'SO2.dat', delimiter='\t', skiprows=3)['Time']
            
            date_time = '1900-01-01 ' + str(times[spec_no])
            
            # Unpack date_time string
            date_time = datetime.datetime.strptime(date_time, '%Y-%m-%d %H:%M:%S')
            
        read_err = False, 'No Error'
        
        # Unpack date and time separately
        date = date_time.date()
        time = date_time.time()
        
    except Exception as e:
        # Something wrong with reading
        x        = None
        y        = None
        date     = None
        time     = None
        spec_no  = None
        read_err = True,e
    
    return x, y, date, time, spec_no, read_err
  
#========================================================================================
#=====================================average_spectra====================================
#========================================================================================
  
def average_spectra(files, spec_type):
    
    for fname in files:
        
        # Check if file exists, if true then read
        try:
          
            # Load spectrum
            grid,y,spec_date,spec_time,spec_no,read_err = read_spectrum(fname,spec_type)
            
            # Sum spectra
            if fname == files[0]:
                spec = y
            else:
                spec += y
         
        except FileNotFoundError:
            grid, spec = None, None
            
    # Divide to get average spectrum
    spec = np.divide(spec, len(files))
    
    return (grid, spec)

