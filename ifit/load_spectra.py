# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 15:58:19 2019

@author: mqbpwbe2
"""

import linecache
import numpy as np
import datetime

#==============================================================================
#================================ read_spectrum ===============================
#==============================================================================

def read_spectrum(fname, spec_type='iFit'):

    '''
    Function to read spectra and extract information, depending on the file
    format

    **Parameters**

    fname : str
        File path to spectrum file

    spec_type : str, optional, default = 'iFit'
        Format of spectrum file. Choices: iFit, Master.Scope, Spectrasuite or
        Basic

    **Returns**

    spectrum : 2D numpy array
        The spectrum as [wavelength, intensities]

    spec_info : dict
        Dictionary of information about the spectrum

    read_err, tuple
        Error message to flag if reading the spectrum fails. Has the format:
        (bool, str), where Boolian = False for no error and True if an error
        occurs. The string is 'No error' for no error, and the error message if
        an error occurs.
    '''

    try:

        # Check the file format exists
        if spec_type not in ['iFit', 'Master.Scope', 'Spectrasuite', 'Basic']:
            raise Exception('File format not recognised')

        if spec_type == 'iFit':

            # Load data into a numpy array
            grid, spec = np.loadtxt(fname, unpack = True)

            # Extract date and time string
            date_time = linecache.getline(fname, 5)[27:].strip()

            # Unpack date_time string
            try:
                date_time = datetime.datetime.strptime(date_time,
                                                       '%Y-%m-%d %H:%M:%S.%f')
            except ValueError:
                date_time = datetime.datetime.strptime(date_time,
                                                       '%Y-%m-%d %H:%M:%S')

            # Get spectrum number
            try:
                spec_no = int(fname[-9:-4])
            except:
                spec_no = 0

        if spec_type == 'Master.Scope':

            # Load data into a numpy array, skipping the header data
            grid, spec = np.genfromtxt(fname, unpack=True, skip_header = 19,
                                       skip_footer = 1)

            # Extract date and time string
            date_time = linecache.getline(fname, 3)[6:].strip()

            # Unpack date_time string
            date_time = datetime.datetime.strptime(date_time,
                                                   '%m-%d-%Y, %H:%M:%S')

            # Get spectrum number
            try:
                spec_no = int(fname[-18:-13])
            except:
                spec_no = 0

        if spec_type == 'Spectrasuite':

            # Load data into a numpy array, skipping the header data
            grid, spec = np.genfromtxt(fname, unpack=True, skip_header = 17,
                                       skip_footer = 2)

            # Extract date and time string
            date_time = linecache.getline(fname, 3).strip()

            # Unpack date_time string
            date_time = datetime.datetime.strptime(date_time,
                                                   '%a %b %Y %H:%M:%S')

            # Get spectrum number
            spec_no = int(fname[-9:-4])

        if spec_type == 'Basic':

            # Generic file format that only holds the essential information

            # Load data into a numpy array
            grid, spec = np.loadtxt(fname, unpack = True, skiprows = 2)

            # Extract date and time string
            read_date = linecache.getline(fname, 1).strip()
            spec_no = int(linecache.getline(fname, 2).strip())

            # Get the date
            date_time = datetime.datetime.strptime(read_date,
                                                   '%Y-%m-%d %H:%M:%S')

        # Report no error
        read_err = False, 'No Error'

        # Unpack date and time separately
        date = date_time.date()
        time = date_time.time()

        spec_info = {'date':    date,
                     'time':    time,
                     'spec_no': spec_no}

    except Exception as e:
        # Something wrong with reading
        grid, spec = np.row_stack([[],[]])
        spec_info = {}
        read_err = True, e

    return grid, spec, spec_info, read_err

#==============================================================================
#============================== average_spectra ===============================
#==============================================================================

def average_spectra(files, spec_type='iFit'):

    '''
    Fuction to average a selection of spectra

    **Parameters**

    files, list
        List of file paths to spectra to read in

    spec_type, str
        Format of spectrum file. Seeread_spectrum for choices

    **Returns**

    grid, array
        Wavelength grid of the spectra

    spec, array
        The averaged spectrum
    '''

    # Load the first spectrum to get the shape
    grid, y, spec_info, read_err = read_spectrum(files[0], spec_type)

    # Create an emptry array to hold the spectra
    spec = np.zeros([len(files), len(grid)])

    for n, fname in enumerate(files):

        # Load spectrum
        grid, y, spec_info, read_err = read_spectrum(fname, spec_type)

        # Add to the spectra array
        spec[n] = y

    # Divide to get average spectrum
    spec = np.average(spec, axis=0)

    return grid, spec

#==============================================================================
#============================== get_spec_details ==============================
#==============================================================================

def get_spec_details(fname, station_name):

    '''
    Function to read in the scanning station details

    **Parameters**
    fname : str
        File path to the file holding the scanner information

    station_name : str
        The name of the station

    **Returns**
    calib_coefs : np.array
        A list of the calibration coefficients
    '''

    # Open the spectrometer setup file
    with open(fname, 'r') as r:

        # Read in the file
        lines = r.readlines()

    # Create a dictionary to hold the spectrometer information
    spec_info = {}

    for line in lines[1:]:

        # Split the line by commas
        line_data = line.strip().split(',')

        # Get the spectrometer station name, serial number and pixel number
        name = line_data[0].strip()
        serial_no = line_data[1].strip()
        pixel_no = line_data[2].strip()

        # Unpack the polynomail coeficients
        poly_coefs = [float(n) for n in line_data[3:]]

        # Add to the dictionary
        spec_info[name] = [serial_no, pixel_no, poly_coefs]

    return spec_info

















