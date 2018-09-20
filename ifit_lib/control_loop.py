# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 11:30:20 2018

@author: mqbpwbe2
"""

import os
import glob
import numpy as np
import datetime as dt
from queue import Queue
from threading import Thread
from seabreeze.cseabreeze.wrapper import SeaBreezeError

from ifit_lib.fit import fit_spec
from ifit_lib.read_spectrum import read_spectrum, average_spectra
from ifit_lib.acquire_spectrum import acquire_spectrum
from ifit_lib.file_control import make_csv_file

#========================================================================================
#=====================================Real Time Setup====================================
#========================================================================================

def rt_setup(self, settings, common, mygui):
    
    '''
    Setup function for iFit in real time analysis mode.
    
    INPUTS
    ------
    settings: dictionary, contains program settings
    common:   dictionary, contains fitting info and parameters
    mygui:    tk.tk(),    Main GUI object
    
    OUTPUTS
    -------
    setup:    boolean,    Reports True is setup was succsessful
    '''
    
    # Check to see if a spectrometer is connected
    if settings['Spectrometer'] == 'Not Connected':
        self.print_output('No spectrometer connected')
        return False
    
    # Get spectrometer serial number to get flat and ILS
    common['spec_name'] = str(self.c_spec.get())
    
    # Check to see if dark spectra have been aquired. If acquired, set in common
    if self.rt_dark_flag == False:
        common['dark_flag'] = 0
        self.print_output('WARNING! No dark spectra aquired!')
    
    else:
        common['dark'] = self.dark_spec

#===================================Get wavelength info==================================

    # Read a single spectrum to get wavelength data
    try:
        x, y, header, t = acquire_spectrum(self, self.spec, 1, 1)
        read_date, read_time = t.split()
        
    except KeyError:
        self.print_output('No spectrometer connected')
        return False
    
    except SeaBreezeError:
        self.print_output('Spectrometer disconnected')
        return False
           
    # Find indices of desired wavelength window and add to common
    common['fit_idx'] = np.where(np.logical_and(settings['wave_start'] <= x, 
                                                x <= settings['wave_stop']))
    common['grid'] = x[common['fit_idx']]
    
    # Find stray light window
    common['stray_idx'] = np.where(np.logical_and(280 <= x, x <= 290))
    
    # If no stray light pixels available, turn off the flag
    if len(common['stray_idx'][0]) == 0:
        common['stray_flag'] = False     
    else:
        common['stray_flag'] = True

#===================================Create output folder=================================
            
    # Create new output folder
    if self.create_out_flag == True:
       
        # Reset loop counter
        self.loop = 0
        
        # Create filename for output file
        self.out_excel_fname = self.results_folder + 'iFit_out.csv'
        
        # Open excel file and write header line
        make_csv_file(self.out_excel_fname, common)
        
        # Create folder to hold spectra
        if not os.path.exists(self.results_folder + 'spectra/'):
                os.makedirs(self.results_folder + 'spectra/')
        
        # Turn off flag to limit folders to one per program run
        self.create_out_flag = False
        
    else:
        
        # Get final spectrum number in folder
        flist = glob.glob(self.results_folder + 'spectra/spectrum*')
        
        # Update loop number to append spectra to those in the folder
        self.loop = int(flist[-1][-9:-4]) + 1
        
        # Create filename for output file
        self.out_excel_fname = self.results_folder + 'iFit_out.csv'                

#==================================Begin analysis loop===================================

    # Set progress bar to the correct format   
    self.progress['mode'] = 'indeterminate'
    self.progress['value'] = 0
    
    return True


#========================================================================================
#===================================Post analysis setup==================================
#========================================================================================

def post_setup(self, settings, common, mygui):
    
    '''
    Setup function for iFit in post analysis mode.
    
    INPUTS
    ------
    settings: dictionary, contains program settings
    common:   dictionary, contains fitting info and parameters
    mygui:    tk.tk(),    Main GUI object
    
    OUTPUTS
    -------
    setup:    boolean,    Reports True is setup was succsessful
    '''
    
    # Get spectrometer serial number to get flat and ILS
    common['spec_name'] = str(self.spec_name.get())

#====================================Read dark spectra===================================
    
    # Read format of spectra
    spec_type = self.spec_type.get()

    # Read in dark spectra
    if common['dark_flag'] == True:
        x, common['dark'], read_err = average_spectra(self.dark_fpaths, spec_type)

        # If there is an error reading the spectrum exit loop
        if read_err[0]:
            self.print_output('Error reading dark spectrum:\n'+str(read_err[1]))
            return False

#===================================Get wavelength info==================================

    # Read first spectrum to get date of data and define stray light indices
    spectrum_data = read_spectrum(self.spec_fpaths[0], spec_type)
    x, y, read_date, read_time, spec_no, read_err = spectrum_data
    
    # If there is an error reading the spectrum exit loop
    if read_err[0]:
        self.print_output('Error reading spectrum:\n' + str(read_err[1]))
        return False
           
    # Find indices of desired wavelength window and add to common
    common['fit_idx'] = np.where(np.logical_and(settings['wave_start'] <= x, 
                                                x <= settings['wave_stop']))
    common['grid'] = x[common['fit_idx']]
    
    # Find stray light window
    common['stray_idx'] = np.where(np.logical_and(280 <= x, x <= 290))
    
    # If no stray light pixels available, turn off the flag
    if len(common['stray_idx'][0]) == 0:
        common['stray_flag'] = False     
    else:
        common['stray_flag'] = True

#===================================Create output folder=================================
            
     # Reset loop counter
    self.loop = 0
    
    # Create filepath to directory to hold program outputs
    self.results_folder = 'Results/iFit/' + str(read_date) + '/'
    
    # Create folder if it doesn't exist
    if not os.path.exists(self.results_folder):
            os.makedirs(self.results_folder)
    
    # Create filename for output file
    self.out_excel_fname = self.results_folder + 'iFit_out.csv'

    try:
        
        # Open excel file and write header line
        make_csv_file(self.out_excel_fname, common)
            
    except PermissionError:
        self.print_output('Please close iFit output file to continue')
        self.build_model_flag = True
        return False                              

#=====================================Set progress bar===================================

    # Set progress bar to the correct format   
    self.progress['mode'] = 'determinate'
    self.progress['value'] = 0
        
    return True



#========================================================================================
#====================================Real time Analysis==================================
#========================================================================================        

def rt_analyse(self, settings, common, mygui):
    
    '''
    Analysis function for iFit in real time mode.
    
    INPUTS
    ------
    common:   dict,     contains fitting info and parameters
    mygui:    tk.tk(),  main GUI object
    
    OUTPUTS
    -------
    spectrum: 2D array, recorded spectrum
    results:  list,     results of the fit
    '''

    # Simultaneously read a spectrum and fit the previous one
    if self.toggle_b.config('text')[-1]=='FITTING ON' and len(common['last_spec']) !=0:

        # Unpack spectrum info
        x, y, timestamp = common['last_spec']
        
        # Create results queue
        result_queue = Queue()
        
        # Create two threads, one to read a spectrum and one to fit
        t1 = Thread(target = acquire_spectrum, args=(self,
                                                     self.spec,
                                                     settings['int_time'],
                                                     int(self.coadds.get()),
                                                     True,
                                                     True,
                                                     result_queue))
        
        t2 = Thread(target = fit_spec, args = (common, [x, y], common['grid'], 
                                               result_queue))
        
        # Initiate threads
        t1.start()
        t2.start()
        
        # Join threads once finished
        t1.join()
        t2.join()
        
        # Get results from threads
        thread_out = {}
        while not result_queue.empty():
            result = result_queue.get()
            thread_out[result[0]] = result[1]
       
        # Unpack thread results
        x, y, header, t = thread_out['spectrum'] 
        fit_results = thread_out['fit']        
        
    # Otherwise just read a spectrum
    else:
        # Read spectrum
        x, y, header, t = acquire_spectrum(self, self.spec, settings['int_time'],
                                           int(self.coadds.get()))
        
        # No fitting
        fit_results = []
        
    # Convert the spectrum read time into a datetime object
    date_time = dt.datetime.strptime(t, '%Y-%m-%d %H:%M:%S.%f')
    date = date_time.date()
    time = date_time.time()
    
    # Build file name
    n = str('{num:05d}'.format(num=self.loop))
    fname = self.results_folder + 'spectra/spectrum_' + n + '.txt'
    
    # Save spectrum
    np.savetxt(fname, np.column_stack((x, y)), header = header)
    
    # Update last spec variable
    common['last_spec'] = [x, y, time]        
    
    return [x,y], fit_results, [self.loop, date, time, fname]




