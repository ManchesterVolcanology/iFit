import os
import logging
import linecache
import numpy as np
from datetime import datetime

logger = logging.getLogger(__name__)


# =============================================================================
# read_spectrum
# =============================================================================

def read_spectrum(fname, spec_type='iFit', wl_calib_file=None,
                  stop_on_err=True):
    """Read spectra and extract metadata.

    Parameters
    ----------
    fname : str
        Spectrum file path.
    spec_type : str, optional
        Format of spectrum file. Choices: iFit, Master.Scope, Spectrasuite,
        mobileDOAS or Basic. Default is "iFit"
    wl_calib_file : str, optional
        The wavelength calibration file to use if not provided in the spectrum
        file. Must be provided for mobileDOAS files.
    stop_on_err : bool, optional
        If True then if an error is found while reading a spectrum an Exception
        is raised.

    Returns
    -------
    grid : numpy array
        The spectrum wavelengths
    spec : numpy array
        The spectrum intensities
    spec_info : dict
        Dictionary of information about the spectrum
    read_err, tuple
        Error message to flag if reading the spectrum fails. Has the format:
        (bool, str), where bool=False for no error and True if an error occurs.
        The string is 'No error' for no error, and the error message if an
        error occurs.
    """
    try:

        if spec_type == 'iFit':

            # Load data into a numpy array
            grid, spec = np.loadtxt(fname, unpack=True)

            # Extract date and time string
            date_time = linecache.getline(fname, 5).strip().split(': ')[-1]

            # Unpack date_time string
            try:
                date_time = datetime.strptime(date_time,
                                              '%Y-%m-%d %H:%M:%S.%f')
            except ValueError:
                date_time = datetime.strptime(date_time, '%Y-%m-%d %H:%M:%S')

            # Get spectrum number
            try:
                spec_no = int(fname[-9:-4])
            except ValueError:
                spec_no = 0

        elif spec_type == 'Master.Scope':

            # Load data into a numpy array, skipping the header data
            grid, spec = np.genfromtxt(fname, unpack=True, skip_header=19,
                                       skip_footer=1)

            # Extract date and time string
            date_time = linecache.getline(fname, 3)[6:].strip()

            # Unpack date_time string
            date_time = datetime.strptime(date_time, '%m-%d-%Y, %H:%M:%S')

            # Get spectrum number
            try:
                spec_no = int(fname[-18:-13])
            except ValueError:
                spec_no = 0

        elif spec_type == 'Spectrasuite':

            # Load data into a numpy array, skipping the header data
            grid, spec = np.genfromtxt(fname, unpack=True, skip_header=17,
                                       skip_footer=2)

            # Extract date and time string
            date_time = linecache.getline(fname, 3).strip()

            # Unpack date_time string
            date_time = datetime.strptime(date_time, '%a %b %Y %H:%M:%S')

            # Get spectrum number
            spec_no = int(fname[-9:-4])

        elif spec_type == 'mobileDOAS':

            # Read the wavelength calibration
            try:
                grid = np.loadtxt(wl_calib_file)
            except FileNotFoundError:
                logger.error(f'Wavelength calibration file {wl_calib_file} '
                             + 'not found')
                raise FileNotFoundError

            # Get the spectrum number
            head, tail = os.path.split(fname)
            try:
                spec_no = int(tail[:5])
            except ValueError:
                spec_no = 0

            # Format created by the mobileDOAS program
            with open(fname, 'r') as r:

                # Read in the file
                lines = [line.strip() for line in r.readlines()]

                # Get the number of pixels
                npixels = int(lines[2])

                # Pull out the spectral information
                spec = np.array([float(y) for y in lines[3:3+npixels]])

                # Pull out the metadata
                metadata = lines[npixels+3:]

                # Get the date
                year, month, day = [int(n) for n in metadata[3].split('.')]
                hour, minute, second = [int(n) for n in metadata[4].split(':')]
                date_time = datetime(year, month, day, hour, minute, second)

        elif spec_type == 'Basic':

            # Generic file format that only holds the essential information

            # Load data into a numpy array
            grid, spec = np.loadtxt(fname, unpack=True, skiprows=2)

            # Extract date and time string
            read_date = linecache.getline(fname, 1).strip()
            spec_no = int(linecache.getline(fname, 2).strip())

            # Get the date
            date_time = datetime.strptime(read_date, '%Y-%m-%d %H:%M:%S')

        # Check the file format exists
        else:
            raise Exception('File format not recognised')

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
        logger.warning(f'Error reading file {fname}:\n{e}')
        grid, spec = np.row_stack([[], []])
        spec_info = {}
        read_err = True, e

        if stop_on_err:
            raise Exception

    return grid, spec, spec_info, read_err


# =============================================================================
# average_spectra
# =============================================================================

def average_spectra(files, spec_type='iFit', wl_calib_file=None):
    """Average a selection of spectra.

    Parameters
    ----------
    files : list
        List of file paths to spectra to read in.
    spec_type : str, optional
        Format of spectrum file. See read_spectrum for choices. Default is
        "iFit"

    Returns
    -------
    spectrum : 2d np.array
        The averaged spectrum as [wavelength, intensity]
    """
    # Load the first spectrum to get the shape
    grid, y, spec_info, read_err = read_spectrum(files[0], spec_type,
                                                 wl_calib_file)

    # Create an emptry array to hold the spectra
    spec = np.zeros([len(files), len(grid)])

    for n, fname in enumerate(files):

        # Load spectrum
        grid, y, spec_info, read_err = read_spectrum(fname, spec_type,
                                                     wl_calib_file)

        # Check for a read error
        if read_err[0]:
            logger.warn(f'Error reading spectrum:\n{read_err[1]}')

        else:
            # Add to the spectra array
            spec[n] = y

    # Divide to get average spectrum
    spec = np.average(spec, axis=0)

    return grid, spec


# =============================================================================
# read_scan
# =============================================================================

def read_scan(fpath, scan_type):
    """Read spectra and header info for a scan block.

    Parameters
    ----------
    fpath : str
        File path to spectra block.
    scan_type : str
        Format of spectrum file. See read_spectrum for choices.

    Returns
    -------
    error : bool
        An error code, 0 if all is OK, 1 if an error was produced
    info_block : array
        Spectra info: spec no, hours, minutes, seconds, motor position
    spec_block : array
        Array of the measured spectra for the scan block
    """
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
            # Scan_MemInt_Time_Count Pixel_Mode Pixel_Mode_Param->Spectral_Data

            # Then spectra, each preceded by an info string
            # 00000 06 56 29 00000

            # Split header row using CR's (ASCII code 10)
            CR_ind = []
            for n, i in enumerate(data[:200]):
                if i == 10:
                    CR_ind.append(n)

            # Unpack indices
            start_data_idx = CR_ind[2] + 1

            # Extract the Vbatt. Vpanel IBatt. Temp.
            # header = data[head_idx0:head_idx1].decode('utf-8').split(' ')

            # Read in all scans
            all_scans = data[start_data_idx:]

            # Calc how many spectra based on file size
            n_spectra = int(len(all_scans)/4131.0)

            # Define arrays
            info_block = np.ndarray([5, n_spectra])
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
                info_block[0, n] = int(scan_info[0:5].decode('utf-8'))
                info_block[1, n] = int(scan_info[6:9].decode('utf-8'))
                info_block[2, n] = int(scan_info[9:12].decode('utf-8'))
                info_block[3, n] = int(scan_info[12:15].decode('utf-8'))
                info_block[4, n] = float(scan_info[15:20].decode('utf-8'))

                # Extract the spectral data and copy into the spec_block
                scan_data = scan[37:]
                spec = np.arange(2046)

                for i in range(2046):

                    spec[i] = int.from_bytes([scan_data[i*2],
                                              scan_data[i*2+1]],
                                             byteorder='big')

                spec_block[n] = spec

            return 0, info_block, spec_block

        except Exception:
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
