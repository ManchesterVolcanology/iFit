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
#================================= read_scan ==================================
#==============================================================================

def read_scan(fpath, scan_type):

    '''
    Function to read spectra and header info for a data block created by
    SpectraLan

    **parameters**

    fpath : str
        File path to data block

    **Returns**

    error : bool
        An error code, 0 if all is OK, 1 if an error was produced

    info_block : array
        Spectra info: spec no, hours, minutes, seconds, motor position

    spec_block : array
        Array of the measured spectra for the scan block
    '''

    if scan_type == 'FLAME':

        try:

            # Get spec details
            with open(fpath, 'rb') as rb:
                data = rb.read()

            # Strucure of file is as follows:
            # First have header as a string, separated by CR's:
            # Vbatt. Vpanel IBatt. Temp.
            # 12.9268 12.9116 -0.0323 29.7
            # nN_Acq Hour Min Sec MotorPos CoAdding Int_Time Ch_Num Scan_Numb
            # Scan_MemInt_Time_Count Pixel_Mode Pixel_Mode_Param -> Spectral_Data

            # Then spectra, each preceded by an info string
            # 00000 06 56 29 00000

            # Split header row using CR's (ASCII code 10)
            CR_ind = []
            for n, i in enumerate(data[:200]):
                if i == 10:
                    CR_ind.append(n)

            # Unpack indices
            #head_idx0 = CR_ind[0] + 1
            #head_idx1 = CR_ind[1] - 1
            start_data_idx = CR_ind[2] + 1

            # Extract the Vbatt. Vpanel IBatt. Temp.
            #header = data[head_idx0:head_idx1].decode('utf-8').split(' ')

            # Read in all scans
            all_scans = data[start_data_idx:]

            # Calc how many spectra based on file size
            n_spectra = int(len(all_scans)/4131.0)

            # Define arrays
            info_block = np.ndarray([5,n_spectra])
            spec_block = np.ndarray([n_spectra, 2046])

            # Loop through each scan and add data to array
            # Header is always 37 bytes and the spectrum is 4094 bytes
            # consisting of 2047 two byte small-endian unsigned integers, miss
            # out the last number

            for n in range(n_spectra):

                # Extract the entire scan from the data

                # Define indices for start and end
                idx0 = n * 4131
                idx1 = idx0 + 4130

                scan = all_scans[idx0:idx1]

                # Separate header from the data
                scan_info = scan[0:36]

                # Extract just string data from the header
                # n_aq, hour, minute, second, motor_pos
                info_block[0,n] = int(scan_info[0:5].decode('utf-8'))
                info_block[1,n] = int(scan_info[6:9].decode('utf-8'))
                info_block[2,n] = int(scan_info[9:12].decode('utf-8'))
                info_block[3,n] = int(scan_info[12:15].decode('utf-8'))
                info_block[4,n] = float(scan_info[15:20].decode('utf-8'))

                # Extract the spectral data and copy into the spec_block
                scan_data = scan[37:]
                spec = np.arange(2046)

                for i in range(2046):

                    spec[i] = int.from_bytes([scan_data[i*2],
                                              scan_data[i*2+1]],
                                             byteorder='big')

                spec_block[n] = spec


            return 0, info_block, spec_block

        except:
            return 1, 1, 1

    elif scan_type == 'OpenSO2':

        try:
            # Read in the numpy file
            data = np.load(fpath)

            # Create empty arrays to hold the spectra
            w, h = data.shape
            info = np.zeros((w, 7))
            spec = np.zeros((w, h - 7))

            # Unpack the data
            for n, line in enumerate(data):

                # Split the spectrum from the spectrum info
                info[n] = data[n][:7]
                spec[n] = data[n][7:]

            return 0, info, spec

        except Exception:
            return 1, 0, 0

#==============================================================================
#============================== get_spec_details ==============================
#==============================================================================

def get_station_info(fname):

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
        pixel_no = int(line_data[2].strip())

        # Unpack the polynomail coeficients
        poly_coefs = [float(n) for n in line_data[3:]]

        # Add to the dictionary
        spec_info[name] = [serial_no, pixel_no, poly_coefs]

    return spec_info

















