# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 10:43:17 2018

@author: mqbpwbe2
"""

import numpy as np
from tkinter import ttk
import tkinter as tk
import tkinter.messagebox as tkMessageBox
import datetime
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog as fd
import seabreeze.spectrometers as sb
from seabreeze.cseabreeze.wrapper import SeaBreezeError

from ifit_lib.file_control import make_directory
from ifit_lib.acquire_spectrum import acquire_spectrum
from ifit_lib.update_graph import update_graph
from ifit_lib.build_gui import make_input
from ifit_lib.read_spectrum import read_spectrum
from ifit_lib.find_nearest import find_nearest

#========================================================================================
#====================================== read_setttings ==================================
#========================================================================================

def read_settings(fname, settings):
    
    '''
    Fuction to read in the settings file
    
    INPUTS
    ------
    fname, str
        File path to settings file
    
    settings, dict
        Dictionary of GUI settings
        
    OUTPUTS
    -------
    settings, dict
        Setting dictionary updated with setings from the file
    '''
    
    # Open the settings file 
    with open(fname, 'r') as r:
                
        # Read data line by line
        data = r.readlines()
        
        # Unpack and save to dictionary
        for i in data:
            name, val, dtype = i.strip().split(';')
            
            # Get the parameter value and change to correct variable type
            if dtype == "<class 'float'>":
                settings[name] = float(val)
                
            
            if dtype == "<class 'int'>":
                settings[name] = int(val)
                
                
            if dtype == "<class 'bool'>":
                if val == 'True':
                    settings[name] = True
                if val == 'False':
                    settings[name] = False
                             
            if dtype == "<class 'str'>":
                settings[name] = str(val)
                    
    return settings

#========================================================================================
#====================================== fit_toggle ======================================
#========================================================================================
            
# Function to toggle fitting on and off
def fit_toggle(self, settings):
    
    '''
    Function to toggle fitting on or off for real time analysis
    
    INPUTS
    ------
    self,
        Program object containing parameters
        
    settings, dict
        Contains the GUI settings
        
    OUTPUTS
    -------
    None
    '''
    
    # Toggle button text and colour
    if self.toggle_b.config('text')[-1] == 'FITTING ON':
        self.toggle_b.config(text = 'FITTING OFF')
        self.toggle_b.config(bg = 'red')
        self.print_output('Fitting turned off\n' +\
                          'Spectrum number ' + str(self.loop))
    else:
        self.toggle_b.config(text = 'FITTING ON')
        self.toggle_b.config(bg = 'green')
        self.print_output('Fitting turned on\n' +\
                          'Spectrum number ' + str(self.loop))

#========================================================================================
#================================ Connect to Spectrometer ===============================
#========================================================================================
        
# Function to connect to the attached spectrometer
def connect_spec(self, settings):
    
    '''
    Fuction to connect to a spectrometer
    
    INPUTS
    ------
    self,
        Program object containing parameters
        
    settings, dict
        Contains the GUI settings
        
    OUTPUTS
    -------
    None
    '''
    
    # Find connected spectrometers
    devices = sb.list_devices()
    
    # If no devices are connected then set string to show. Else assign first to spec
    if len(devices) == 0:
        self.spec = 0
        settings['Spectrometer'] = 'Not Connected'
        devices = ['Not Connected']
        self.print_output('No devices found')
    else:
        try:
            # Connect to spectrometer
            self.spec = sb.Spectrometer(devices[0])
            
            # Set intial integration time
            self.spec.integration_time_micros(float(self.int_time.get())*1000)
            
            # Record serial number in settings
            settings['Spectrometer'] = str(self.spec.serial_number)
            
            self.print_output('Spectrometer '+settings['Spectrometer']+' Connected')
            
            # Create filepath to directory to hold program outputs
            results_folder = 'Results/iFit/'+str(datetime.date.today())+'/ifit_out/'
           
            # Create folder
            self.results_folder = make_directory(results_folder, overwrite = True)
            
            # Create notes file
            self.notes_fname = self.results_folder + 'notes.txt'
            with open(self.notes_fname, 'w') as w:
                w.write('Notes file for iFit\n\n')
    
            
        except SeaBreezeError:
            self.print_output('Spectrometer already open')
        
    # Update text to show spectrometer name
    self.c_spec.set(settings['Spectrometer'])
        
#========================================================================================
#==================================== Filepath Buttons ==================================
#========================================================================================
   
# Function to select spectra to analyse
def spec_fp(self):
    
    '''
    Function to get plume spectra file paths
    
    INPUTS
    ------
    self,
        Program object containing parameters
        
    OUTPUTS
    -------
    None
    '''
        
    # Open dialouge to get files
    fpaths = fd.askopenfilenames()

    if fpaths != '':
        self.spec_fpaths = []
        for i in fpaths:
            self.spec_fpaths.append(str(i))
        
        # Save output to input 
        self.spec_ent.set(str(len(self.spec_fpaths)) + ' spectra selected')

# Function to select dark spectra
def dark_fp(self):
    
    '''
    Function to get dark spectra file paths
    
    INPUTS
    ------
    self,
        Program object containing parameters
        
    OUTPUTS
    -------
    None
    '''
        
    # Open dialouge to get files
    fpaths = fd.askopenfilenames()
    
    if fpaths != '':
        self.dark_fpaths = []
        for i in fpaths:
            self.dark_fpaths.append(str(i))
        
        # Save output to input
        self.dark_ent.set(str(len(self.dark_fpaths)) + ' spectra selected')
        
#========================================================================================
#================================ Update Integration Time ===============================
#========================================================================================

def update_int_time(self, settings):
    
    '''
    Function to update the spectrometer integration time
    
    INPUTS
    ------
    self,
        Program object containing parameters
        
    settings, dict
        Dictionary of GUI settings
        
    OUTPUTS
    -------
    None
    '''
    
    try:
        # Update integration time on spectrometer
        settings['int_time'] = float(self.int_time.get())
        self.spec.integration_time_micros(settings['int_time']*1000)
        
        self.print_output('Integration time updated to '+str(settings['int_time']) +\
                          '\nSpectrum number ' + str(self.loop))
        
    except KeyError:
        self.print_output('No spectrometer conected')
        
#========================================================================================
#============================== Read a single test spectrum =============================
#========================================================================================

def test_spec(self, settings, mygui, line, ax):
    
    '''
    Function to take a test spectrum
    
    INPUTS
    ------
    self,
        Program object containing parameters
        
    settings, dict
        Dictionary of GUI settings
        
    my_gui, tk.Tk() object
        Main GUI object
        
    line, pyplot.plot object
        Line on which to plot the test spectrum
        
    ax, pyplot.axis object
        Axis containing the line
        
    OUTPUTS
    -------
    None
    '''
    
    # Update status
    self.status.set('Acquiring')
    mygui.update(self)
    
    x, y, header, read_time = acquire_spectrum(self,
                                               self.spec, 
                                               settings['int_time'], 
                                               int(self.coadds.get()))
    
    # Display the test spectrum
    line.set_data(x, y)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    self.canvas.draw()
    
    # Update status
    self.status.set('Standby')        

#========================================================================================
#======================================= Read Darks =====================================
#========================================================================================

def read_darks(self, settings, mygui, line, ax):
    
    '''
    Function to take a dark spectrum
    
    INPUTS
    ------
    self,
        Program object containing parameters
        
    settings, dict
        Dictionary of GUI settings
        
    my_gui, tk.Tk() object
        Main GUI object
        
    line, pyplot.plot object
        Line on which to plot the dark spectrum
        
    ax, pyplot.axis object
        Axis containing the line
        
    OUTPUTS
    -------
    None
    '''
    
    if tkMessageBox.askyesno('Record Darks', 'Ready to begin measuring darks?'):
    
        # Update status
        self.status.set('Acquiring darks')
        
        # Reset progress bar
        self.progress['mode'] = 'determinate'
        self.progress['value'] = 0
        mygui.update(self)
        
        # Create zero array to hold dark spectra
        dark = np.zeros(2048)
        
        # Define dark filepath
        dark_fp = self.results_folder + 'dark/'
        
        # Create the directory  
        dark_fp = make_directory(dark_fp)
        
        # Loop over number of darks to collect
        dark_n = int(self.no_darks.get())
        
        for i in range(dark_n):
            
            # Create filename
            n = str('{num:05d}'.format(num=i))
            filepath = dark_fp + 'spectrum_' + n + '.txt'
            
            # Read spectrometer
            try:
                spc = acquire_spectrum(self, self.spec, settings['int_time'],
                                       int(self.coadds.get()))
            except KeyError:
                self.print_output('No spectrometer connected')
                return
            
            except SeaBreezeError:
                self.print_output('Spectrometer disconnected')
                return 
            
            # Unpack spectrum data
            x, y, header, read_time = spc
            
            # Save
            np.savetxt(filepath, np.column_stack((x,y)), header = header)
            
            # Sum up the darks
            dark = np.add(dark, y)
            
            # Update the progress bar
            prog = ((i+1)/dark_n) * 100
            self.progress['value'] = prog
            
            # Force gui to update
            mygui.update(self)
        
        # Divide by number of darks to get average
        self.dark_spec = np.divide(dark, dark_n)

        # Display the dark spectrum
        lines = [line]
        axes =  [ax]
        
        # Build data array to pass to graphing function
        #                 x data    y data          x limits     y limits
        data = np.array(([x,        self.dark_spec, 'auto',      'auto']))
        
        # Update graph
        update_graph(lines, axes, self.canvas, data)
        
        # Update notes file
        self.print_output('Dark updated\n' + \
                          'Spectrum no: ' + str(self.loop) + '\n' + \
                          'No. darks: ' + str(dark_n) + '\n' + \
                          'Integration time (ms): '+str(settings['int_time'])+'\n'+\
                          'Coadds: ' + str(self.coadds.get()))
         
        # Set dark_flag to True
        self.rt_dark_flag = True
        
        # Update status
        self.status.set('Standby')
        
#========================================================================================
#========================================== Stop ========================================
#========================================================================================
        
def stop(self):
    
    '''
    Function to stop acquisition and/or analysis
    
    INPUTS
    ------
    self,
        Program object containing parameters
        
    OUTPUTS
    -------
    None
    '''
    
    self.stop_flag = True
    self.print_output('Loop Stopped\nSpectrum number ' + str(self.loop))

#========================================================================================
#=================================== Advanced Setings ===================================
#========================================================================================

# Define some fonts to use in the program
NORM_FONT = ('Verdana', 8)
MED_FONT  = ('Veranda', 11)
LARG_FONT = ('Verdana', 12, 'bold')

def adv_settings(self, settings, version):
    
    '''
    Function to take a test spectrum
    
    INPUTS
    ------
    self,
        Program object containing parameters
        
    settings, dict
        Dictionary of GUI settings
        
    version, str
        Version of the program. Controls whether plot labels are updated or not
        
    OUTPUTS
    -------
    None
    '''
    
    # Make popup window
    popup = tk.Tk()
    tk.Tk.wm_title(popup, 'Advanced Settings')
    
    # Make updating/closing function
    def update_settings(settings, close_flag):
        
        # Update model settings
        settings['wave_start']       = popup.wave_start.get()
        settings['wave_stop']        = popup.wave_stop.get()
        settings['model_pad']        = popup.model_pad.get()
        settings['model_res']        = popup.model_res.get()
        settings['dark_flag']        = popup.dark_b.get()
        settings['flat_flag']        = popup.flat_b.get()
        settings['fit_weight']       = popup.fit_weight.get()
        settings['solar_resid_flag'] = popup.resid_b.get()
        settings['update_params']    = popup.update_b.get()
        settings['good_fit_bound']   = popup.fit_bound.get()
        
        # Update parameter settings
        settings['so2_amt']          = popup.so2_amt.get()
        settings['no2_amt']          = popup.no2_amt.get()
        settings['o3_amt']           = popup.o3_amt.get()
        settings['bro_amt']          = popup.bro_amt.get()
        settings['shift']            = popup.shift.get()
        settings['ils_width']        = popup.ils_width.get()
        settings['gauss_weight']     = popup.gauss_weight.get()
        settings['ldf']              = popup.ldf.get()
        settings['poly_n']           = popup.poly_n.get()
        settings['stretch']          = popup.stretch.get()
        settings['ring_amt']         = popup.ring_amt.get()
        
        settings['Fit so2']          = popup.so2_amt_c.get()
        settings['Fit no2']          = popup.no2_amt_c.get()
        settings['Fit o3']           = popup.o3_amt_c.get()
        settings['Fit bro']          = popup.bro_amt_c.get()
        settings['Fit shift']        = popup.shift_c.get()
        settings['Fit ILS']          = popup.ils_width_c.get()
        settings['Fit LDF']          = popup.ldf_c.get()
        settings['Fit stretch']      = popup.stretch_c.get()
        settings['Fit ring']         = popup.ring_amt_c.get()
        
        # Update data base settings
        settings['sol_path']         = popup.sol_path.get()
        settings['ring_path']        = popup.ring_path.get()
        settings['so2_path']         = popup.so2_path.get()
        settings['no2_path']         = popup.no2_path.get()
        settings['o3_path']          = popup.o3_path.get()
        settings['o3_temp']          = popup.o3_temp.get()
        settings['bro_path']         = popup.bro_path.get()
        settings['solar_resid_path'] = popup.solar_resid_path.get()
        
        # Update graph settings
        settings['Show Graphs']      = graph_b.get()
        settings['scroll_flag']      = scroll_b.get()
        settings['scroll_spec_no']   = popup.spec_no.get()
        settings['x_plot']           = x_plot.get()
        settings['resid_type']       = resid_type.get()
        settings['analysis_gas']     = gas.get()
        
        if version == 'iFit':
        
            # Translate gas choice to parameter names and graph print
            gas_choice = {'so2' : r'SO$_2$',
                          'no2' : r'NO$_2$',
                          'o3'  : r'O$_3$' ,
                          'bro' : 'BrO'    ,
                          'ring': 'Ring'   }
            
            # Update graph
            if settings['resid_type'] == 'Percentage':
                self.ax2.set_ylabel('Fit residual (%)', fontsize=10)
            if settings['resid_type'] == 'Absolute':
                self.ax2.set_ylabel('Fit residual (Abs)', fontsize=10)
            if settings['resid_type'] == 'Spec/Fit':
                self.ax2.set_ylabel('Fit residual (Spec/Fit)', fontsize=10)
                
            self.ax3.set_ylabel(gas_choice[gas.get()] + ' Absorbance', fontsize=10)
            self.ax4.set_ylabel(gas_choice[gas.get()] + ' amt (ppm.m)', fontsize = 10)
            
            if settings['x_plot'] == 'Number':
                self.ax4.set_xlabel('Spectrum number', fontsize=10)
            if settings['x_plot'] == 'Time':
                self.ax4.set_xlabel('Time (decimal hours)', fontsize=10)
                
            self.canvas.draw()
        
        # Print message to confirm settings have been updated
        self.print_output('Settings updated')
        
        # Close the window
        if close_flag == True:
            popup.destroy()
        
    # Make function to select filepaths
    def update_fp(entry):
        
        # Open dialouge to get files
        fpath = fd.askopenfilenames()
        
        if fpath != '':
            entry.set(str(fpath[0]))   
        
            # Turn on making fwd model
            self.build_model_flag = True
        
    # Make function to turn on re-building the forward model
    def build_fwd_model(*args):
        self.build_model_flag = True
        
#========================================================================================
#===================================== Make frames ======================================
#========================================================================================        

    # Create notebook to hold different frames
    nb = ttk.Notebook(popup)
    model_frame = ttk.Frame(nb)
    datab_frame = ttk.Frame(nb)
    graph_frame = ttk.Frame(nb)
    
    # Add the frames to the notebook
    nb.add(model_frame, text = 'Model Settings')
    nb.add(graph_frame, text = 'Graph Settings')
    nb.add(datab_frame, text = 'Data Base Settings')
    
    # Add the notebook to the window
    nb.grid(row = 0, column=0, padx=10, pady=10, sticky = 'NW', columnspan = 10)

    # Add aditional frame for the buttons
    button_frame = ttk.Frame(popup)
    button_frame.grid(row=1, column=0, padx = 10, pady = 10, sticky="EW")

    # Button to apply changes and close
    b1 = ttk.Button(button_frame, text='Ok', command=lambda: update_settings(settings, 
                                                                             True))
    b1.grid(row = 0, column = 0, padx = 70, pady = 5, sticky = "EW")
    
    # Button to just apply changes
    b2 = ttk.Button(button_frame, text='Apply', command=lambda: update_settings(settings,
                                                                                False))
    b2.grid(row = 0, column = 1, padx = 70, pady = 5, sticky = "EW")
    
    # Buttong to cancel
    b3 = ttk.Button(button_frame, text='Cancel', command=lambda: popup.destroy())
    b3.grid(row = 0, column = 2, padx = 70, pady = 5, sticky = "EW")

#========================================================================================
#==================================== Model Settings ====================================
#========================================================================================

#===================================== Model Setup ======================================

    # Create frame to hold wavelength settings
    wl_frame = ttk.Frame(model_frame)
    wl_frame.grid(row=0, column=0, padx=10, pady=10, columnspan=10)

    # Create row and column number counters
    row_n = 0
    col_n = 0
    
    # Create entry for start and stop wavelengths
    fit_range_l = tk.Label(wl_frame, text = 'Fit Range:', font = NORM_FONT)
    fit_range_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    
    popup.wave_start = tk.DoubleVar(wl_frame, value = settings['wave_start'])
    make_input(frame = wl_frame, 
               text = 'Fit Range:', 
               row = row_n, column = col_n, 
               var = popup.wave_start, 
               input_type = 'Entry',
               width = 12,
               command = build_fwd_model)
    
    popup.wave_stop = tk.DoubleVar(wl_frame, value = settings['wave_stop'])
    make_input(frame = wl_frame, 
               text = 'to:', 
               row = row_n, column = col_n+2, 
               var = popup.wave_stop, 
               input_type = 'Entry',
               width = 12,
               command = build_fwd_model)
    
    # Set model grid padding
    popup.model_pad = tk.DoubleVar(wl_frame, value = settings['model_pad'])
    make_input(frame = wl_frame, 
               text = 'Model Grid\nPadding (nm):', 
               row = row_n, column = col_n+4, 
               var = popup.model_pad, 
               input_type = 'Entry',
               width = 12,
               command = build_fwd_model)
    
    # Set resolution of model grid
    popup.model_res = tk.DoubleVar(wl_frame, value = settings['model_res'])
    make_input(frame = wl_frame, 
               text = 'Model Grid\nSpacing (nm):', 
               row = row_n, column = col_n+6, 
               var = popup.model_res, 
               input_type = 'Entry',
               width = 12,
               command = build_fwd_model)
    row_n += 1
    
    # Control whether or not to remove dark spectra
    popup.dark_b = tk.BooleanVar(model_frame, value = settings['dark_flag'])
    make_input(frame = model_frame, 
               text = 'Remove Dark\nSpectrum?', 
               row = row_n, column = col_n, 
               var = popup.dark_b, 
               input_type = 'Checkbutton')
    row_n += 1
    
    # Control whether or not to remove flat spectra
    popup.flat_b = tk.BooleanVar(model_frame, value = settings['flat_flag'])
    make_input(frame = model_frame, 
               text = 'Remove Flat\nSpectrum?', 
               row = row_n, column = col_n, 
               var = popup.flat_b, 
               input_type = 'Checkbutton',
               command = build_fwd_model)
    row_n += 1
    
    # Control how fit weighting is applied
    weight_options = [settings['fit_weight'],
                      'SO2 xsec',
                      'Noise',
                      'None']
    popup.fit_weight = tk.StringVar(model_frame, value = settings['fit_weight'])
    make_input(frame = model_frame, 
               text = 'Fit\nWeighting:', 
               row = row_n, column = col_n, 
               var = popup.fit_weight, 
               input_type = 'OptionMenu',
               options = weight_options,
               width = 8)
    row_n += 1
    
    # Control whether or not to remove, form or ignore the solar residual spectrum
    resid_options = [settings['solar_resid_flag'],
                     'Ignore',
                     'Generate',
                     'Remove']
    popup.resid_b = tk.StringVar(model_frame, value = settings['solar_resid_flag'])
    make_input(frame = model_frame, 
               text = 'Solar\nresidual:', 
               row = row_n, column = col_n, 
               var = popup.resid_b, 
               input_type = 'OptionMenu',
               options = resid_options,
               command = build_fwd_model,
               width = 8)
    row_n += 1
    
    # Control whether or not to update fit parameter guesses with the last fit values
    popup.update_b = tk.BooleanVar(model_frame, value = settings['update_params'])
    make_input(frame = model_frame, 
               text = 'Auto-update\nfit parameters?', 
               row = row_n, column = col_n, 
               var = popup.update_b, 
               input_type = 'Checkbutton')
    row_n += 1
    
    # Control bound of goodness of fit
    popup.fit_bound = tk.DoubleVar(model_frame, value = settings['good_fit_bound'])
    make_input(frame = model_frame, 
               text = 'Good Fit\nBound (%):', 
               row = row_n, column = col_n, 
               var = popup.fit_bound, 
               input_type = 'Entry',
               width = 12)
    
    # Create separator
    sep = ttk.Separator(model_frame, orient='vertical')
    sep.grid(column=2, row = 1, rowspan=7, sticky='ns', padx=10, pady=5)

#====================================== Parameters ======================================
    
    # Reset row counter
    row_n = 1
    col_n = 3
    
    # Create array of fitting options
    fit_options = ['', 'Fit', 'Fix', 'N/A']
    
    # Gas Amounts
    popup.so2_amt = tk.DoubleVar(model_frame, value = settings['so2_amt'])
    make_input(frame = model_frame, 
               text = 'SO2:', 
               row = row_n, column = col_n, 
               var = popup.so2_amt, 
               input_type = 'Entry',
               width = 12)
    popup.so2_amt_c = tk.StringVar(model_frame, value = settings['Fit so2'])
    fit_options[0] = settings['Fit so2']
    so2_c = ttk.OptionMenu(model_frame, popup.so2_amt_c, *fit_options)
    so2_c.config(width = 6)
    so2_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1
    
    popup.no2_amt = tk.DoubleVar(model_frame, value = settings['no2_amt'])
    make_input(frame = model_frame, 
               text = 'NO2:', 
               row = row_n, column = col_n, 
               var = popup.no2_amt, 
               input_type = 'Entry',
               width = 12)
    popup.no2_amt_c = tk.StringVar(model_frame, value = settings['Fit no2'])
    fit_options[0] = settings['Fit no2']
    no2_c = ttk.OptionMenu(model_frame, popup.no2_amt_c, *fit_options)
    no2_c.config(width = 6)
    no2_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1
    
    popup.o3_amt = tk.DoubleVar(model_frame, value = settings['o3_amt'])
    make_input(frame = model_frame, 
               text = 'O3:', 
               row = row_n, column = col_n, 
               var = popup.o3_amt, 
               input_type = 'Entry',
               width = 12)
    popup.o3_amt_c = tk.StringVar(model_frame, value = settings['Fit o3'])
    fit_options[0] = settings['Fit o3']
    o3_c = ttk.OptionMenu(model_frame, popup.o3_amt_c, *fit_options)
    o3_c.config(width = 6)
    o3_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1
   
    popup.bro_amt = tk.DoubleVar(model_frame, value = settings['bro_amt'])
    make_input(frame = model_frame, 
               text = 'BrO:', 
               row = row_n, column = col_n, 
               var = popup.bro_amt, 
               input_type = 'Entry',
               width = 12)
    popup.bro_amt_c = tk.StringVar(model_frame, value = settings['Fit bro'])
    fit_options[0] = settings['Fit bro']
    bro_c = ttk.OptionMenu(model_frame, popup.bro_amt_c, *fit_options)
    bro_c.config(width = 6)
    bro_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5) 
    row_n += 1   
    
    # Spectrometer ILS width
    popup.ils_width = tk.DoubleVar(model_frame, value = settings['ils_width'])
    ils_width_l, ils_width_e = make_input(frame = model_frame, 
                                          text = 'ILS\nWidth:', 
                                          row = row_n, column = col_n, 
                                          var = popup.ils_width, 
                                          input_type = 'Entry',
                                          width = 12)
    ils_width_e.config(state='disabled')
    
    def ils_activate_check(*args):
        
        # Enable/disable the ils width input
        if popup.ils_width_c.get() in ['File', 'Shape']:
            ils_width_e.config(state='disabled')
        
            # Turn on fwd model creation
            self.build_model_flag = True
            
        else:
            ils_width_e.config(state='normal')
    
    popup.ils_width_c = tk.StringVar(model_frame, value = settings['Fit ILS'])
    ils_fit_options = [settings['Fit ILS'], 'Fit', 'Fix', 'File', 'Shape']
    ils_width_c = ttk.OptionMenu(model_frame, popup.ils_width_c, *ils_fit_options,
                                 command = ils_activate_check)
    ils_width_c.config(width = 6)
    ils_width_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1
    ils_activate_check()
     
    # ILS Gaussian weighting
    popup.gauss_weight = tk.DoubleVar(model_frame, value = settings['gauss_weight'])
    make_input(frame = model_frame, 
               text = 'ILS Gauss\nWeight:', 
               row = row_n, column = col_n, 
               var = popup.gauss_weight, 
               input_type = 'Entry',
               width = 12)    
    row_n += 1 
    
    
    
    # Create separator
    sep = ttk.Separator(model_frame, orient='vertical')
    sep.grid(column=6, row = 1, rowspan=7, sticky='ns', padx=10, pady=5)
    
#====================================== New Column ======================================
    
    # Reset row counter
    row_n = 1
    col_n = 7
    
    # Polynomial coefficents
    popup.poly_n = tk.IntVar(model_frame, value = int(settings['poly_n']))
    make_input(frame = model_frame, 
               text = 'Poly\nOrder:', 
               row = row_n, column = col_n, 
               var = popup.poly_n, 
               input_type = 'Spinbox',
               vals = [0, 10],
               width = 12)
    popup.poly_n.set(settings['poly_n'])
    row_n += 1

    # Spectrometer wavelength shift
    popup.shift = tk.DoubleVar(model_frame, value = settings['shift'])
    shift_l, shift_e = make_input(frame = model_frame, 
                                  text = 'Shift:', 
                                  row = row_n, column = col_n, 
                                  var = popup.shift, 
                                  input_type = 'Entry',
                                  width = 12)
    
    def shift_activate_check(*args):
        if popup.shift_c.get() == 'Pre-calc':
            shift_e.config(state='disabled')
        else:
            shift_e.config(state='normal')
    
    popup.shift_c = tk.StringVar(model_frame, value = settings['Fit shift'])
    shift_options = [settings['Fit shift'], 'Fit', 'Fix', 'N/A', 'Pre-calc']
    shift_c = ttk.OptionMenu(model_frame, popup.shift_c, *shift_options,
                             command = shift_activate_check)
    shift_c.config(width = 6)
    shift_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1 
    shift_activate_check() 
    
    # Spectrometer wavelength stretch
    popup.stretch = tk.DoubleVar(model_frame, value = settings['stretch'])
    make_input(frame = model_frame, 
               text = 'Stretch:', 
               row = row_n, column = col_n, 
               var = popup.stretch, 
               input_type = 'Entry',
               width = 12)
    popup.stretch_c = tk.StringVar(model_frame, value = settings['Fit stretch'])
    fit_options[0] = settings['Fit stretch']
    stretch_c = ttk.OptionMenu(model_frame, popup.stretch_c, *fit_options)
    stretch_c.config(width = 6)
    stretch_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1
    
    # Ring effect
    popup.ring_amt = tk.DoubleVar(model_frame, value = settings['ring_amt'])
    make_input(frame = model_frame, 
               text = 'Ring:', 
               row = row_n, column = col_n, 
               var = popup.ring_amt, 
               input_type = 'Entry',
               width = 12)
    popup.ring_amt_c = tk.StringVar(model_frame, value = settings['Fit ring'])
    fit_options[0] = settings['Fit ring']
    ring_c = ttk.OptionMenu(model_frame, popup.ring_amt_c, *fit_options)
    ring_c.config(width = 6)
    ring_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1
    
    # Light dilution factor
    popup.ldf = tk.DoubleVar(model_frame, value = settings['ldf'])
    make_input(frame = model_frame, 
               text = 'LDF:', 
               row = row_n, column = col_n, 
               var = popup.ldf, 
               input_type = 'Entry',
               width = 12)
    popup.ldf_c = tk.StringVar(model_frame, value = settings['Fit LDF'])
    fit_options[0] = settings['Fit LDF']
    ldf_c = ttk.OptionMenu(model_frame, popup.ldf_c, *fit_options)
    ldf_c.config(width = 6)
    ldf_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1 
    
    
    
    
    
    
#========================================================================================
#================================= Data base file paths =================================
#========================================================================================
     
    # Create row number counter
    row_n = 0
    col_n = 0
    
    # Solar spectrum
    popup.sol_path = tk.StringVar(datab_frame, value = settings['sol_path'])
    make_input(frame = datab_frame, 
               text = 'Solar Spectrum:', 
               row = row_n, column = col_n, 
               var = popup.sol_path, 
               input_type = 'Entry',
               width = 40,
               command = build_fwd_model)
    sol_path_b = ttk.Button(datab_frame, text = "Browse", 
                            command = lambda: update_fp(popup.sol_path))
    sol_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # Ring spectrum
    popup.ring_path = tk.StringVar(datab_frame, value = settings['ring_path'])
    make_input(frame = datab_frame, 
               text = 'Ring Spectrum:', 
               row = row_n, column = col_n, 
               var = popup.ring_path, 
               input_type = 'Entry',
               width = 40,
               command = build_fwd_model)
    ring_path_b = ttk.Button(datab_frame, text = "Browse", 
                             command = lambda: update_fp(popup.ring_path))
    ring_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # SO2 xsec
    popup.so2_path = tk.StringVar(datab_frame, value = settings['so2_path'])
    make_input(frame = datab_frame, 
               text = 'SO2 xsec:', 
               row = row_n, column = col_n, 
               var = popup.so2_path, 
               input_type = 'Entry',
               width = 40,
               command = build_fwd_model)
    so2_path_b = ttk.Button(datab_frame, text = "Browse", 
                            command = lambda: update_fp(popup.so2_path))
    so2_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # NO2 xsec
    popup.no2_path = tk.StringVar(datab_frame, value = settings['no2_path'])
    make_input(frame = datab_frame, 
               text = 'NO2 xsec:', 
               row = row_n, column = col_n, 
               var = popup.no2_path, 
               input_type = 'Entry',
               width = 40,
               command = build_fwd_model)
    no2_path_b = ttk.Button(datab_frame, text = "Browse", 
                            command = lambda: update_fp(popup.no2_path))
    no2_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # O3 xsec
    popup.o3_path = tk.StringVar(datab_frame, value = settings['o3_path'])
    make_input(frame = datab_frame, 
               text = 'O3 Path:', 
               row = row_n, column = col_n, 
               var = popup.o3_path, 
               input_type = 'Entry',
               width = 40,
               command = build_fwd_model)
    o3_path_b = ttk.Button(datab_frame, text = "Browse", 
                           command = lambda: update_fp(popup.o3_path))
    o3_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    
    # O3 Temp
    temp_options = [settings['o3_temp'], '298K', '283K', '273K', '263K', '253K', '243K', 
                    '233K', '223K', '213K', '203K', '193K']
    popup.o3_temp = tk.StringVar(datab_frame, value = settings['o3_temp'])
    o3_temp_c = ttk.OptionMenu(datab_frame, popup.o3_temp, *temp_options, 
                               command = build_fwd_model)
    o3_temp_c.grid(row = row_n, column = 3, padx = 5, pady = 5)
    row_n += 1
       
    # Bro xsec
    popup.bro_path = tk.StringVar(datab_frame, value = settings['bro_path'])
    make_input(frame = datab_frame, 
               text = 'BrO xsec:', 
               row = row_n, column = col_n, 
               var = popup.bro_path, 
               input_type = 'Entry',
               width = 40,
               command = build_fwd_model) 
    bro_path_b = ttk.Button(datab_frame, text = "Browse", 
                            command = lambda: update_fp(popup.bro_path))
    bro_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # Solar residual
    popup.solar_resid_path = tk.StringVar(datab_frame,value=settings['solar_resid_path'])
    make_input(frame = datab_frame, 
               text = 'Solar Residual:', 
               row = row_n, column = col_n, 
               var = popup.solar_resid_path, 
               input_type = 'Entry',
               width = 40,
               command = build_fwd_model)
    solar_resid_path_b = ttk.Button(datab_frame, text = "Browse", 
                                    command = lambda: update_fp(popup.solar_resid_path))
    solar_resid_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1    
    
#========================================================================================
#=================================== Graph Settings =====================================
#========================================================================================


    # Create row number counter
    row_n = 0
    
    # Control graph display settings
    graph_b = tk.BooleanVar(graph_frame, value = settings['Show Graphs'])
    make_input(frame = graph_frame, 
               text = 'Show Graphs?', 
               row = row_n, column = col_n, 
               var = graph_b, 
               input_type = 'Checkbutton')
    row_n += 1
    
    # Select which gas to analyse
    gas_options = [settings['analysis_gas'],
                   'so2',
                   'no2',
                   'o3',
                   'bro',
                   'ring']
    gas = tk.StringVar(graph_frame, value = settings['scroll_flag'])
    make_input(frame = graph_frame, 
               text = 'Parameter\nto View:', 
               row = row_n, column = col_n, 
               var = gas, 
               input_type = 'OptionMenu',
               options = gas_options)
    row_n += 1
    
    # Set x-axis for time series
    x_plot_options = [settings['x_plot'],
                      'Time',
                      'Number']
    x_plot = tk.StringVar(graph_frame, value = settings['resid_type'])
    make_input(frame = graph_frame, 
               text = 'Time Series\nx-axis:', 
               row = row_n, column = col_n, 
               var = x_plot, 
               input_type = 'OptionMenu',
               options = x_plot_options)
    x_plot_l = tk.Label(graph_frame, text = 'Time Series\nx-axis', font = NORM_FONT)
    x_plot_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    x_plot_m = ttk.OptionMenu(graph_frame, x_plot, *x_plot_options)
    x_plot_m.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1 
    
    # Turn on/off graph scrolling
    scroll_b = tk.BooleanVar(graph_frame, value = settings['scroll_flag'])
    make_input(frame = graph_frame, 
               text = 'Scoll Graphs?', 
               row = row_n, column = col_n, 
               var = scroll_b, 
               input_type = 'Checkbutton')
    row_n += 1
    
    # Set number of spectra to display on graph
    popup.spec_no = tk.IntVar(graph_frame, value = settings['scroll_spec_no'])
    make_input(frame = graph_frame, 
               text = 'No. Spectra\nto display', 
               row = row_n, column = col_n, 
               var = popup.spec_no, 
               input_type = 'Entry')
    row_n += 1
    
    # Set format of residual
    resid_options = [settings['resid_type'],
                     'Absolute',
                     'Percentage',
                     'Spec/Fit']
    resid_type = tk.StringVar(graph_frame, value = settings['resid_type'])
    make_input(frame = graph_frame, 
               text = 'Residual Display:', 
               row = row_n, column = col_n, 
               var = resid_type, 
               input_type = 'OptionMenu',
               options = resid_options)
    row_n += 1

#========================================================================================   
#================================ Measure flat spectrum =================================
#======================================================================================== 

def meas_ils():
    
    # Define gaussian equation
    def gaussian_fit(x,a,b,c):
        
        # Calculate pure gaussian
        return np.multiply(np.exp(-(np.power(np.subtract(x, b), 2))/(2 * c**2)), a)
    
    # Define function to get the ILS
    def begin():
        
        # Read in spectrum
        x, y, date, time, spec_no, read_err = read_spectrum(fpath.get(), spec_type.get())
        
        # Plot the spectrum
        line0.set_data(x, y)
        ax0.relim()
        ax0.autoscale_view()
        
        # Find peaks
        p, props = find_peaks(y,
                              prominence = 500,
                              width = 5,
                              rel_height = 500)
        
        # Find peak nearest to the 302.728 nm line
        i, val = find_nearest(x[p], 302.728)  
        
        # Extract that line
        start, stop = x[p[i]] - 3, x[p[i]] + 3
        idx = np.where(np.logical_and(x > start, x < stop))
        
        grid = x[idx]
        spec = y[idx]
        
        # Centre the baseline and normalise
        spec = spec - min(spec)
        spec = np.divide(spec, np.sum(spec))
        
        # Fit a guassian
        guess = [max(spec), 302.728, 0.5]
        popt, pcov = curve_fit(gaussian_fit, grid, spec, p0=guess)
        
        # Unpack results
        a, b, c = popt
        a_err, b_err, c_err = np.sqrt(np.diag(pcov))
        
        # Plot
        line1.set_data(grid, spec)
        line2.set_data(grid, gaussian_fit(grid, *popt))
        ax1.relim()
        ax1.autoscale_view()
        
        canvas.draw()
        
        header = 'ILS for spectrometer ' + serial.get() + '\n' + \
                 'Wavelength (nm), Intensity'
        
        np.savetxt(save_path.get() + 'ils_' + serial.get() + '.txt', 
                   np.column_stack((x, y)), header = header)
        
        with open(save_path.get() + 'ils_width_' + serial.get() + '.txt', 'w') as w:
            w.write('# ILS Gaussian FWHM for spectrometer ' + serial.get() + '\n')
            w.write('# Measured ILS is ' + str(2.355*c) + ' +/- ' + str(2.355*c_err) + '\n')
            w.write(str(2.355*c))
    
    # Create window
    win = tk.Tk()
    tk.Tk.wm_title(win, 'Measure spectrometer ILS')
    
    # Create notebook to hold different frames
    nb = ttk.Notebook(win)
    page1 = ttk.Frame(nb)
    page2 = ttk.Frame(nb)
    
    # Create two frames, one for post analysis and one for real time acquisition
    nb.add(page1, text = 'Load')
    nb.add(page2, text = 'Measure')
    
    nb.grid(column=0, row = 0, padx=10, pady=10, sticky = 'NW')
    
    # Create frames
    control_frame = ttk.Frame(win)
    control_frame.grid(row = 1, column = 0, padx = 10, pady = 10)
    graph_frame = ttk.Frame(win)
    graph_frame.grid(row = 0, column = 1, padx = 10, pady = 10, rowspan = 10)
    
#==================================== Program setup =====================================
    
    ### LOADING SPECTRA ###
    
    # Filepath to the spectrum of the HG lamp
    fpath = tk.StringVar(page1, value = '')
    make_input(frame = page1,
               text = 'Spectrum\nPath:',
               row = 0, column = 0,
               var = fpath,
               input_type = 'Entry',
               width  = 40)
    
    def get_fpath():
        fpath.set(fd.askopenfilename())
        
    ttk.Button(page1, text = 'Select', command = get_fpath
               ).grid(row = 0, column = 2, padx = 5, pady = 5)
    
    # Set spectrometer serial no
    serial = tk.StringVar(win, value = '')
    make_input(frame = page1,
               text = 'Spectrometer\nSerial No.:',
               row = 1, column = 0,
               var = serial,
               input_type = 'Entry',
               width  = 40)
    
    # Create entry to select spectra type
    spec_options = ['iFit',
                    'iFit',
                    'Master.Scope',
                    'Jai Spec',
                    'Spectrasuite',
                    'GSJ',
                    'Ind']
    
    spec_type = tk.StringVar(page1, value = 'iFit')
    make_input(frame = page1, 
               text = 'Spectrum\nType:', 
               row = 2, column = 0, 
               var = spec_type, 
               input_type = 'OptionMenu',
               options = spec_options, 
               sticky = 'W')
    
    ### MEASURING SPECTRA ###
    
    # 
    
#=================================== Program control ====================================
    
    # Set the savepath
    save_path = tk.StringVar(control_frame, value = '')
    make_input(frame = control_frame,
               text = 'Save path:',
               row = 0, column = 0,
               var = save_path,
               input_type = 'Entry',
               width  = 40)
    
    def get_save_path():
        save_path.set(fd.askdirectory() + '/')
        
    ttk.Button(control_frame, text = 'Select', command = get_save_path
               ).grid(row = 0, column = 2, padx = 5, pady = 5)
    
    # Create button to begin
    ttk.Button(control_frame, text = 'Begin', command = begin
               ).grid(row = 2, column = 0, padx = 5, pady = 5, columnspan = 3)
    
#======================================== Graphs ========================================
    
    # Create graph for display
    fig = plt.figure(figsize = (6,4))
    ax0 = fig.add_subplot(211)
    ax1 = fig.add_subplot(212)
    
    ax0.set_ylabel('Intensity (counts)', fontsize = 10)
    ax1.set_xlabel('Wavelength (nm)', fontsize = 10)
    ax1.set_ylabel('Intensity (counts)', fontsize = 10)
    
    line0, = ax0.plot(0, 0)
    line1, = ax1.plot(0, 0)
    line2, = ax1.plot(0, 0)
    
    plt.tight_layout()
    
    # Create the canvas to hold the graph in the GUI
    canvas = FigureCanvasTkAgg(fig, graph_frame)
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady = 10)
    

#========================================================================================   
#===================================== Measure ILS ======================================
#======================================================================================== 

def meas_flat():
    
    print('Not yet operational')   
    
#========================================================================================   
#================================== Gas Xsec analyser ===================================
#======================================================================================== 

def conv_xsec():
    
    print('Not yet operational')
    '''
    Function to convert a gas cross section. Not functional yet
    '''
    
    # Make popup window
    win = tk.Tk()
    tk.Tk.wm_title(win, 'Cross Section Analyser')
    
    # Create frames
    fpath_frame = ttk.Frame(win)
    fpath_frame.grid(row=0, column=0, padx=10, sticky = 'N')
    
    graph_frame = ttk.Frame(win)
    graph_frame.grid(row=0, column=1, padx=10, rowspan=20, sticky = 'N')
    
    setup_frame = ttk.Frame(win)
    setup_frame.grid(row=1, column=0, padx=10, sticky = 'N')
    
    # Make function to select filepaths
    def update_fp(entry):
        
        # Open dialouge to get files
        fpath = fd.askopenfilenames()
        
        if fpath != '':
            entry.set(str(fpath[0])) 
            
    # Create array of intergers for spinboxes
    choices = []
    for i in range(101):
        choices.append(i)
   
#================================== Filepath controls ===================================
    
    # Load spectrum fpath
    win.load_path = tk.StringVar(fpath_frame)
    load_path_l = tk.Label(fpath_frame, text = 'Load File Path:', font = NORM_FONT)
    load_path_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
    load_path_e = ttk.Entry(fpath_frame, text=win.load_path, font=NORM_FONT, width=30)
    load_path_e.grid(row = 0, column = 1, padx = 5, pady = 5)
    load_path_b = ttk.Button(fpath_frame, text = "Browse", 
                             command = lambda: update_fp(win.load_path))
    load_path_b.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = 'W')
    
    # Save spectrum fpath
    win.save_path = tk.StringVar(fpath_frame)
    save_path_l = tk.Label(fpath_frame, text = 'Save file path:', font = NORM_FONT)
    save_path_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
    save_path_e = ttk.Entry(fpath_frame, text=win.load_path, font=NORM_FONT, width=30)
    save_path_e.grid(row = 1, column = 1, padx = 5, pady = 5)
    save_path_b = ttk.Button(fpath_frame, text = "Browse", 
                             command = lambda: update_fp(win.save_path))
    save_path_b.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
   
#=================================== Program controls ===================================
    
    # Create inputs for number of header and footer rows to ignore
    win.head_rows = tk.IntVar(setup_frame, value = 0)
    head_rows_l = tk.Label(setup_frame, text = 'Header\nrows:', font = NORM_FONT)
    head_rows_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
    head_rows_e = tk.Spinbox(setup_frame, values = choices, width = 12,
                             textvariable = win.head_rows)
    head_rows_e.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = 'W')
       
    win.foot_rows = tk.IntVar(setup_frame, value = 0)
    foot_rows_l = tk.Label(setup_frame, text = 'Footer\nrows:', font = NORM_FONT)
    foot_rows_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
    foot_rows_e = tk.Spinbox(setup_frame, values = choices, width = 12,
                             textvariable = win.foot_rows)
    foot_rows_e.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = 'W')
    
    # Create inputs for columns for wavelength and xsec
    win.wl_col = tk.IntVar(setup_frame, value = 0)
    wl_col_l = tk.Label(setup_frame, text = 'Wavelength\nColumn:', font = NORM_FONT)
    wl_col_l.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = 'W')
    wl_col_e = tk.Spinbox(setup_frame, values = choices, width = 12,
                          textvariable = win.wl_col)
    wl_col_e.grid(row = 0, column = 3, padx = 5, pady = 5, sticky = 'W')
    
    win.xsec_col = tk.IntVar(setup_frame, value = 0)
    xsec_col_l = tk.Label(setup_frame, text = 'Spectrum\nColumn:', font = NORM_FONT)
    xsec_col_l.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
    xsec_col_e = tk.Spinbox(setup_frame, values = choices, width = 12,
                            textvariable = win.xsec_col)
    xsec_col_e.grid(row = 1, column = 3, padx = 5, pady = 5, sticky = 'W')
    
    # Create control on wether the grid is nm or wavenumber
    grid_options = ['nanometer','nanometer', 'wavenumber']
    win.unit = tk.StringVar(setup_frame, value = grid_options[0])
    unit_l = tk.Label(setup_frame, text = 'Wavlength\nUnits:', font = NORM_FONT)
    unit_l.grid(row = 2, column = 2, padx = 5, pady = 5)
    unit_c = ttk.OptionMenu(setup_frame, win.unit, *grid_options)
    unit_c.grid(row = 2, column = 3, padx = 5, pady = 5)
    
    # Create input for column delimiter
    delim_options = ['tab','tab','space','comma','input']
    win.delimiter = tk.StringVar(setup_frame, value = 0)
    delimiter_l = tk.Label(setup_frame, text = 'Delimiter:', font = NORM_FONT)
    delimiter_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
    delimiter_e = ttk.OptionMenu(setup_frame, win.delimiter, *delim_options)
    delimiter_e.grid(row = 2, column = 1, padx = 5, pady = 5, sticky = 'W')

#=================================== Convert button =====================================
    
    # Create frame to hold the button
    button_frame = ttk.Frame(win)
    button_frame.grid(row=2, column=0, sticky='N')
    
    # Create button to convert the selected spectrum
    conv_b = ttk.Button(button_frame, text = 'Convert', command = convert)
    conv_b.grid()

#==================================== Graph control =====================================
    
    # Create graph for display
    fig = plt.figure(figsize = (6,4))
    ax = fig.add_subplot(111)
    
    ax.set_xlabel('Wavelength (nm)', fontsize = 10)
    ax.set_ylabel('Absorption (photons/cm2)', fontsize = 10)
    
    line, = ax.plot(0, 0)
    
    plt.tight_layout()
    
    # Create the canvas to hold the graph in the GUI
    canvas = FigureCanvasTkAgg(fig, graph_frame)
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=0, rowspan = 20, padx=10, pady = 10)
    
def convert():
    pass
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    