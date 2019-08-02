# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 14:30:41 2017

@author: mqbpwbe2
"""

import linecache
import numpy as np
import datetime

#========================================================================================
#===================================== read_spectrum ====================================
#========================================================================================

# Reads spectrum file from spectro_gui

def read_spectrum(fname, spec_type = 'iFit'):

    '''
    Function to read spectra and extract information, depending on the file format

    INPUTS
    ------
    fname, str
        File path to spectrum file

    spec_type, str
        Format of spectrum file. Choices: iFit, Master.Scope, Spectrasuite or Basic

    OUTPUTS
    -------
    x, numpy array
        Spectrum wavelength data

    y, numpy array
        Spectrum intensity data

    date, datetime.date object
        Date of measurement

    time, datetime.time object
        Time of measurement

    spec_no, float
        Spectrum number, used for display and recording results

    read_err, tuple
        Error message to flag if reading the spectrum fails. Has the format:
        (bool, str), where Boolian = False for no error and True if an error occurs. The
        string is 'No error' for no error, and the error message if an error occurs.
    '''

    try:

        # Check the file format exists
        if spec_type not in ['iFit', 'Master.Scope', 'Spectrasuite', 'Basic']:
            raise Exception('File format not recognised. Must be one of "iFit", ' +\
                            '"Master.Scope", "Spectrasuite" or "Basic"')

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
            try:
                spec_no = int(fname[-9:-4])
            except:
                spec_no = 0

        if spec_type == 'Master.Scope':

            # Load data into a numpy array, skipping the header data
            x,y = np.genfromtxt(fname, unpack=True, skip_header = 19, skip_footer = 1)

            # Extract date and time string
            date_time = linecache.getline(fname, 3)[6:].strip()

            # Unpack date_time string
            date_time = datetime.datetime.strptime(date_time, '%m-%d-%Y, %H:%M:%S')

            # Get spectrum number
            try:
                spec_no = int(fname[-18:-13])
            except:
                spec_no = 0

        if spec_type == 'Spectrasuite':

            # Load data into a numpy array, skipping the header data
            x,y = np.genfromtxt(fname, unpack=True, skip_header = 17, skip_footer = 2)

            # Extract date and time string
            date_time = linecache.getline(fname, 3).strip()

            # Unpack date_time string
            date_time = datetime.datetime.strptime(date_time, '%a %b %Y %H:%M:%S')

            # Get spectrum number
            spec_no = int(fname[-9:-4])

        if spec_type == 'Basic':

            # Generic file format that only holds the essential information

            # Load data into a numpy array
            x, y = np.loadtxt(fname, unpack = True, skiprows = 2)

            # Extract date and time string
            read_date = linecache.getline(fname, 1).strip()
            spec_no = int(linecache.getline(fname, 2).strip())

            # Get the date
            date_time = datetime.datetime.strptime(read_date, '%Y-%m-%d %H:%M:%S')

        # Report no error
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
        read_err = True, e

    return x, y, date, time, spec_no, read_err

#========================================================================================
#==================================== average_spectra ===================================
#========================================================================================

def average_spectra(files, spec_type):

    '''
    Fuction to average a selection of spectra

    INPUTS
    ------
    files, list
        List of file paths to spectra to read in

    spec_type, str
        Format of spectrum file. Choices: iFit, Master.Scope, Jai Spec, Spectrasuite, GSJ

    OUTPUTS
    -------
    grid, array
        Wavelength grid of the spectra

    spec, array
        The resulting averaged spectrum

    read_err, tuple
        Error message to flag if reading the spectrum fails. Has the format:
        (bool, str), where Boolian = False for no error and True if an error occurs. The
        string is 'No error' for no error, and the error message if an error occurs.

    '''

    for fname in files:

        # Load spectrum
        grid, y, spec_date, spec_time, spec_no, read_err = read_spectrum(fname,spec_type)

        # Check if spectra were read correctly
        if not read_err[0]:

            # Sum spectra
            if fname == files[0]:
                spec = y
            else:
                spec += y

        else:
            return None, None, read_err

    # Divide to get average spectrum
    spec = np.divide(spec, len(files))

    return grid, spec, read_err

