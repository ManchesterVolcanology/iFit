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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog as fd
import seabreeze.spectrometers as sb
from seabreeze.cseabreeze.wrapper import SeaBreezeError

from ifit_lib.file_control import make_directory
from ifit_lib.acquire_spectrum import acquire_spectrum
from ifit_lib.update_graph import update_graph

#========================================================================================
#=======================================Read Setttings===================================
#========================================================================================

def read_settings(fname, settings):
    
    # Open the settings file 
    with open(fname, 'r') as r:
                
        # Read data line by line
        data = r.readlines()
        
        # Unpack and save to dictionary
        for i in data:
            name, val = i.strip().split(';')
            
            # Check if float
            try:
                settings[name] = float(val)
            except ValueError:
                
                # Check if boolean
                if val in ['True', 'False']:
                    # Convert to boolian
                    settings[name] = bool(val)
               
                else:
                    # Set as string
                    settings[name] = val
                    
    return settings

#========================================================================================
#=======================================Toggle fitting===================================
#========================================================================================
            
# Function to toggle fitting on and off
def fit_toggle(self, settings):
    
    # Toggle button text and colour
    if self.toggle_button.config('text')[-1] == 'FITTING ON':
        self.toggle_button.config(text = 'FITTING OFF')
        self.toggle_button.config(bg = 'red')
        self.print_output('Fitting turned off\n' +\
                          'Spectrum number ' + str(self.loop))
    else:
        self.toggle_button.config(text = 'FITTING ON')
        self.toggle_button.config(bg = 'green')
        self.print_output('Fitting turned on\n' +\
                          'Spectrum number ' + str(self.loop))

#========================================================================================
#=================================Connect to Spectrometer================================
#========================================================================================
        
# Function to connect to the attached spectrometer
def connect_spec(self, settings):
    
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
            self.rt_folder = make_directory(results_folder, overwrite = True)
            
            # Create notes file
            self.notes_fname = self.rt_folder + 'notes.txt'
            with open(self.notes_fname, 'w') as w:
                w.write('Notes file for iFit\n\n')
    
            
        except SeaBreezeError:
            self.print_output('Spectrometer already open')
        
    # Update text to show spectrometer name
    self.c_spec.set(settings['Spectrometer'])
        
#========================================================================================
#=====================================Filepath Buttons===================================
#========================================================================================
   
# Function to select spectra to analyse
def spec_fp(self):
        
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
        
    # Open dialouge to get files
    fpaths = fd.askopenfilenames()
    
    if fpaths != '':
        self.dark_fpaths = []
        for i in fpaths:
            self.dark_fpaths.append(str(i))
        
        # Save output to input
        self.dark_ent.set(str(len(self.dark_fpaths)) + ' spectra selected')
        
#========================================================================================
#=================================Update Integration Time================================
#========================================================================================

def update_int_time(self, settings):
    
    try:
        # Update integration time on spectrometer
        settings['int_time'] = float(self.int_time.get())
        self.spec.integration_time_micros(settings['int_time']*1000)
        
        self.print_output('Integration time updated to '+str(settings['int_time']) +\
                          '\nSpectrum number ' + str(self.loop))
        
    except KeyError:
        self.print_output('No spectrometer conected')
        
#========================================================================================
#===============================Read a single test spectrum==============================
#========================================================================================

def test_spec(self, settings, mygui, line, ax):
    
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
#========================================Read Darks======================================
#========================================================================================

def read_darks(self, settings, mygui, line, ax):
    
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
        dark_fp = self.rt_folder + 'dark/'
        
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
#===========================================Stop=========================================
#========================================================================================
        
def stop(self, settings):
    self.stop_flag = True
    self.print_output('Loop Stopped\nSpectrum number ' + str(self.loop))

#========================================================================================
#====================================Advanced Setings====================================
#========================================================================================

# Define some fonts to use in the program
NORM_FONT = ('Verdana', 8)
MED_FONT  = ('Veranda', 11)
LARG_FONT = ('Verdana', 12, 'bold')

def adv_settings(self, settings, version):
    
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
            
            # Update graph
            if settings['resid_type'] == 'Percentage':
                self.ax2.set_ylabel('Fit residual (%)', fontsize=10)
            if settings['resid_type'] == 'Absolute':
                self.ax2.set_ylabel('Fit residual (Abs)', fontsize=10)
            if settings['resid_type'] == 'Spec/Fit':
                self.ax2.set_ylabel('Fit residual (Spec/Fit)', fontsize=10)
                
            self.ax3.set_ylabel(settings['analysis_gas'] + ' Absorbance', fontsize=10)
            
            self.ax4.set_ylabel(settings['analysis_gas'] + ' amt (ppm.m)', fontsize = 10)
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
#======================================Make frames=======================================
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
#=====================================Model Settings=====================================
#========================================================================================

#======================================Model Setup=======================================

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
    wave_start_e = ttk.Entry(wl_frame, textvariable = popup.wave_start, width = 12, 
                             validate="focusout", validatecommand = build_fwd_model)
    wave_start_e.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    
    popup.wave_stop = tk.DoubleVar(wl_frame, value = settings['wave_stop'])
    wave_stop_l = tk.Label(wl_frame, text = 'to:', font = NORM_FONT)
    wave_stop_l.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    wave_stop_e = ttk.Entry(wl_frame, textvariable = popup.wave_stop, width = 12, 
                            validate="focusout",validatecommand = build_fwd_model)
    wave_stop_e.grid(row = row_n, column = col_n+3, padx = 5, pady = 5)
    
    # Set model grid padding
    popup.model_pad = tk.DoubleVar(wl_frame, value = settings['model_pad'])
    model_pad_l = tk.Label(wl_frame, text='Model Grid\nPadding (nm):', font=NORM_FONT)
    model_pad_l.grid(row = row_n, column = col_n+4, padx = 5, pady = 5)
    model_pad_e = ttk.Entry(wl_frame, textvariable = popup.model_pad, width = 12,
                            validate="focusout", validatecommand = build_fwd_model)
    model_pad_e.grid(row = row_n, column = col_n+5, padx = 5, pady = 5)
    
    # Set resolution of model grid
    popup.model_res = tk.DoubleVar(wl_frame, value = settings['model_res'])
    model_res_l = tk.Label(wl_frame, text='Model Grid\nSpacing (nm):', font=NORM_FONT)
    model_res_l.grid(row = row_n, column = col_n+6, padx = 5, pady = 5)
    model_res_e = ttk.Entry(wl_frame, textvariable = popup.model_res, width = 12,
                            validate="focusout", validatecommand = build_fwd_model)
    model_res_e.grid(row = row_n, column = col_n+7, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove dark spectra
    popup.dark_b = tk.BooleanVar(model_frame, value = settings['dark_flag'])
    dark_l = tk.Label(model_frame, text = 'Remove Dark\nSpectrum?', font = NORM_FONT)
    dark_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    dark_c = ttk.Checkbutton(model_frame, variable = popup.dark_b)
    dark_c.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove flat spectra
    popup.flat_b = tk.BooleanVar(model_frame, value = settings['flat_flag'])
    flat_l = tk.Label(model_frame, text = 'Remove Flat\nSpectrum?', font = NORM_FONT)
    flat_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    flat_c = ttk.Checkbutton(model_frame, variable = popup.flat_b, 
                             command = build_fwd_model)
    flat_c.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    row_n += 1
    
    # Control how fit weighting is applied
    weight_options = [settings['fit_weight'],
                      'SO2 xsec',
                      'Noise',
                      'None']
    popup.fit_weight = tk.StringVar(model_frame, value = settings['fit_weight'])
    fit_weight_l = tk.Label(model_frame, text = 'Fit\nWeighting:', font = NORM_FONT)
    fit_weight_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    fit_weight_c = ttk.OptionMenu(model_frame, popup.fit_weight, *weight_options)
    fit_weight_c.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove, form or ignore the solar residual spectrum
    resid_options = [settings['solar_resid_flag'],
                     'Ignore',
                     'Generate',
                     'Remove']
    popup.resid_b = tk.StringVar(model_frame, value = settings['solar_resid_flag'])
    resid_l = tk.Label(model_frame, text = 'Solar\nresidual:', font = NORM_FONT)
    resid_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    resid_c = ttk.OptionMenu(model_frame, popup.resid_b, *resid_options, 
                             command = build_fwd_model)
    resid_c.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to update fit parameter guesses with the last fit values
    popup.update_b = tk.BooleanVar(model_frame, value = settings['update_params'])
    update_l = tk.Label(model_frame, text='Auto-update\nfit parameters?', font=NORM_FONT)
    update_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    update_c = ttk.Checkbutton(model_frame, variable = popup.update_b)
    update_c.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    row_n += 1
    
    # Control bound of goodness of fit
    popup.fit_bound = tk.DoubleVar(model_frame, value = settings['good_fit_bound'])
    fit_bound_l = tk.Label(model_frame, text = 'Good Fit\nBound (%):', font = NORM_FONT)
    fit_bound_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    fit_bound_e = ttk.Entry(model_frame, textvariable = popup.fit_bound, width = 12)
    fit_bound_e.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    row_n += 1
    
    # Create separator
    sep = ttk.Separator(model_frame, orient='vertical')
    sep.grid(column=2, row = 1, rowspan=7, sticky='ns', padx=10, pady=5)

#=======================================Parameters=======================================
    
    # Reset row counter
    row_n = 1
    col_n = 3
    
    # Create array of fitting options
    fit_options = ['', 'Fit', 'Fix', 'N/A']
    
    # Gas Amounts
    popup.so2_amt = tk.DoubleVar(model_frame, value = settings['so2_amt'])
    so2_amt_l = tk.Label(model_frame, text = 'SO2:', font = NORM_FONT)
    so2_amt_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    so2_amt_e = ttk.Entry(model_frame, textvariable = popup.so2_amt, width = 12)
    so2_amt_e.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    popup.so2_amt_c = tk.StringVar(model_frame, value = settings['Fit so2'])
    fit_options[0] = settings['Fit so2']
    so2_c = ttk.OptionMenu(model_frame, popup.so2_amt_c, *fit_options)
    so2_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1
    
    popup.no2_amt = tk.DoubleVar(model_frame, value = settings['no2_amt'])
    no2_amt_l = tk.Label(model_frame, text = 'NO2:', font = NORM_FONT)
    no2_amt_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    no2_amt_e = ttk.Entry(model_frame, textvariable = popup.no2_amt, width = 12)
    no2_amt_e.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    popup.no2_amt_c = tk.StringVar(model_frame, value = settings['Fit no2'])
    fit_options[0] = settings['Fit no2']
    no2_c = ttk.OptionMenu(model_frame, popup.no2_amt_c, *fit_options)
    no2_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1
    
    popup.o3_amt = tk.DoubleVar(model_frame, value = settings['o3_amt'])
    o3_amt_l = tk.Label(model_frame, text = 'O3:', font = NORM_FONT)
    o3_amt_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    o3_amt_e = ttk.Entry(model_frame, textvariable = popup.o3_amt, width = 12)
    o3_amt_e.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    popup.o3_amt_c = tk.StringVar(model_frame, value = settings['Fit o3'])
    fit_options[0] = settings['Fit o3']
    o3_c = ttk.OptionMenu(model_frame, popup.o3_amt_c, *fit_options)
    o3_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1
   
    popup.bro_amt = tk.DoubleVar(model_frame, value = settings['bro_amt'])
    bro_amt_l = tk.Label(model_frame, text = 'BrO:', font = NORM_FONT)
    bro_amt_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    bro_amt_e = ttk.Entry(model_frame, textvariable = popup.bro_amt, width = 12)
    bro_amt_e.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    popup.bro_amt_c = tk.StringVar(model_frame, value = settings['Fit bro'])
    fit_options[0] = settings['Fit bro']
    bro_c = ttk.OptionMenu(model_frame, popup.bro_amt_c, *fit_options)
    bro_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5) 
    row_n += 1   
    
    # Spectrometer ILS width
    popup.ils_width = tk.DoubleVar(model_frame, value = settings['ils_width'])
    ils_width_l = tk.Label(model_frame, text = 'ILS\nWidth:', font = NORM_FONT)
    ils_width_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    ils_width_e = ttk.Entry(model_frame, textvariable = popup.ils_width, width = 12)
    ils_width_e.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    ils_width_e.config(state='disabled')
    
    def ils_activate_check(*args):
        
        # Enable/disable the ils width input
        if popup.ils_width_c.get() == 'File':
            ils_width_e.config(state='disabled')
        
            # Turn on fwd model creation
            self.build_model_flag = True
            
        else:
            ils_width_e.config(state='normal')
    
    popup.ils_width_c = tk.StringVar(model_frame, value = settings['Fit ILS'])
    ils_fit_options = [settings['Fit ILS'], 'Fit', 'Fix', 'File']
    ils_width_c = ttk.OptionMenu(model_frame, popup.ils_width_c, *ils_fit_options,
                                 command = ils_activate_check)
    ils_width_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1
    ils_activate_check()
     
    # ILS Gaussian weighting
    popup.gauss_weight = tk.DoubleVar(model_frame, value = settings['gauss_weight'])
    gauss_weight_l = tk.Label(model_frame, text = 'ILS Gauss\nWeight:', 
                              font = NORM_FONT)
    gauss_weight_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    gauss_weight_e = ttk.Entry(model_frame, textvariable = popup.gauss_weight,
                               width = 12)
    gauss_weight_e.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    row_n += 1 
    
    # Create separator
    sep = ttk.Separator(model_frame, orient='vertical')
    sep.grid(column=6, row = 1, rowspan=7, sticky='ns', padx=10, pady=5)
    
#=======================================New Column=======================================
    
    # Reset row counter
    row_n = 1
    col_n = 7
    
    # Polynomial coefficents
    popup.poly_n = tk.IntVar(model_frame, value = int(settings['poly_n']))
    poly_n_vals = [0,1,2,3,4,5,6,7,8,9,10]
    poly_n_l = tk.Label(model_frame, text = 'Poly\nOrder:', font = NORM_FONT)
    poly_n_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    poly_n_e = tk.Spinbox(model_frame, values = poly_n_vals, width = 12,
                          textvariable = popup.poly_n)
    poly_n_e.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    popup.poly_n.set(settings['poly_n'])
    row_n += 1

    # Spectrometer wavelength shift
    popup.shift = tk.DoubleVar(model_frame, value = settings['shift'])
    shift_l = tk.Label(model_frame, text = 'Shift:', font = NORM_FONT)
    shift_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    shift_e = ttk.Entry(model_frame, textvariable = popup.shift, width = 12)
    shift_e.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    
    def shift_activate_check(*args):
        if popup.shift_c.get() == 'Pre-calc':
            shift_e.config(state='disabled')
        else:
            shift_e.config(state='normal')
    
    popup.shift_c = tk.StringVar(model_frame, value = settings['Fit shift'])
    shift_options = [settings['Fit shift'], 'Fit', 'Fix', 'N/A', 'Pre-calc']
    shift_c = ttk.OptionMenu(model_frame, popup.shift_c, *shift_options,
                             command = shift_activate_check)
    shift_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1 
    shift_activate_check() 
    
    # Spectrometer wavelength stretch
    popup.stretch = tk.DoubleVar(model_frame, value = settings['stretch'])
    stretch_l = tk.Label(model_frame, text = 'Stretch:', font = NORM_FONT)
    stretch_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    stretch_e = ttk.Entry(model_frame, textvariable = popup.stretch, width = 12)
    stretch_e.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    popup.stretch_c = tk.StringVar(model_frame, value = settings['Fit stretch'])
    fit_options[0] = settings['Fit stretch']
    stretch_c = ttk.OptionMenu(model_frame, popup.stretch_c, *fit_options)
    stretch_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1
    
    # Ring effect
    popup.ring_amt = tk.DoubleVar(model_frame, value = settings['ring_amt'])
    ring_amt_l = tk.Label(model_frame, text = 'Ring:', font = NORM_FONT)
    ring_amt_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    ring_amt_e = ttk.Entry(model_frame, textvariable = popup.ring_amt, width = 12)
    ring_amt_e.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    popup.ring_amt_c = tk.StringVar(model_frame, value = settings['Fit ring'])
    fit_options[0] = settings['Fit ring']
    ring_c = ttk.OptionMenu(model_frame, popup.ring_amt_c, *fit_options)
    ring_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1
    
    # Light dilution factor
    popup.ldf = tk.DoubleVar(model_frame, value = settings['ldf'])
    ldf_l = tk.Label(model_frame, text = 'LDF:', font = NORM_FONT)
    ldf_l.grid(row = row_n, column = col_n, padx = 5, pady = 5)
    ldf_e = ttk.Entry(model_frame, textvariable = popup.ldf, width = 12)
    ldf_e.grid(row = row_n, column = col_n+1, padx = 5, pady = 5)
    popup.ldf_c = tk.StringVar(model_frame, value = settings['Fit LDF'])
    fit_options[0] = settings['Fit LDF']
    ils_width_c = ttk.OptionMenu(model_frame, popup.ldf_c, *fit_options)
    ils_width_c.grid(row = row_n, column = col_n+2, padx = 5, pady = 5)
    row_n += 1 
    
    
    
    
    
    
#========================================================================================
#==================================Data base file paths==================================
#========================================================================================
     
    # Create row number counter
    row_n = 0
    
    # Solar spectrum
    popup.sol_path = tk.StringVar(datab_frame, value = settings['sol_path'])
    sol_path_l = tk.Label(datab_frame, text = 'Solar spectrum:', font = NORM_FONT)
    sol_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    sol_path_e = ttk.Entry(datab_frame, text=popup.sol_path, font=NORM_FONT, width=40)
    sol_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    sol_path_b = ttk.Button(datab_frame, text = "Browse", 
                            command = lambda: update_fp(popup.sol_path))
    sol_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # Ring spectrum
    popup.ring_path = tk.StringVar(datab_frame, value = settings['ring_path'])
    ring_path_l = tk.Label(datab_frame, text = 'Ring Spectrum:', font = NORM_FONT)
    ring_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    ring_path_e = ttk.Entry(datab_frame, textvariable=popup.ring_path, font=NORM_FONT, 
                            width=40)
    ring_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    ring_path_b = ttk.Button(datab_frame, text = "Browse", 
                             command = lambda: update_fp(popup.ring_path))
    ring_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # SO2 xsec
    popup.so2_path = tk.StringVar(datab_frame, value = settings['so2_path'])
    so2_path_l = tk.Label(datab_frame, text = 'SO2 xsec:', font = NORM_FONT)
    so2_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    so2_path_e = ttk.Entry(datab_frame, textvariable = popup.so2_path, font = NORM_FONT,
                           width = 40)
    so2_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    so2_path_b = ttk.Button(datab_frame, text = "Browse", 
                            command = lambda: update_fp(popup.so2_path))
    so2_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # NO2 xsec
    popup.no2_path = tk.StringVar(datab_frame, value = settings['no2_path'])
    no2_path_l = tk.Label(datab_frame, text = 'NO2 xsec:', font = NORM_FONT)
    no2_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    no2_path_e = ttk.Entry(datab_frame, textvariable = popup.no2_path, font = NORM_FONT,
                           width = 40)
    no2_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    no2_path_b = ttk.Button(datab_frame, text = "Browse", 
                            command = lambda: update_fp(popup.no2_path))
    no2_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # O3 xsec
    popup.o3_path = tk.StringVar(datab_frame, value = settings['o3_path'])
    o3_path_l = tk.Label(datab_frame, text = 'O3 xsec:', font = NORM_FONT)
    o3_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    o3_path_e = ttk.Entry(datab_frame, textvariable = popup.o3_path, font = NORM_FONT,
                          width = 40)
    o3_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
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
    bro_path_l = tk.Label(datab_frame, text = 'BrO xsec:', font = NORM_FONT)
    bro_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    bro_path_e = ttk.Entry(datab_frame, textvariable = popup.bro_path, font = NORM_FONT,
                           width = 40)
    bro_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5) 
    bro_path_b = ttk.Button(datab_frame, text = "Browse", 
                            command = lambda: update_fp(popup.bro_path))
    bro_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # Solar residual
    popup.solar_resid_path = tk.StringVar(datab_frame,value=settings['solar_resid_path'])
    solar_resid_path_l = tk.Label(datab_frame, text = 'Solar Residual:', font=NORM_FONT)
    solar_resid_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    solar_resid_path_e = ttk.Entry(datab_frame, textvariable = popup.solar_resid_path,
                                   font = NORM_FONT, width = 40)
    solar_resid_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5) 
    solar_resid_path_b = ttk.Button(datab_frame, text = "Browse", 
                                    command = lambda: update_fp(popup.solar_resid_path))
    solar_resid_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1    
    
    # Create a button to convert a gas cross-section to be converted
    conv_xsec_b = ttk.Button(datab_frame, text = 'Convert Xsec', command = conv_xsec)
    conv_xsec_b.grid(row = row_n, column = 1, padx=5, pady=5)
    
#========================================================================================
#====================================Graph Settings======================================
#========================================================================================


    # Create row number counter
    row_n = 0
    
    # Control graph display settings
    graph_b = tk.BooleanVar(graph_frame, value = settings['Show Graphs'])
    graph_l = tk.Label(graph_frame, text = 'Show Graphs?', font = NORM_FONT)
    graph_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    graph_c = ttk.Checkbutton(graph_frame, variable = graph_b)
    graph_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Select which gas to analyse
    gas_options = [settings['analysis_gas'],
                   'SO2',
                   'NO2',
                   'O3',
                   'BrO',
                   'Ring']
    gas = tk.StringVar(graph_frame, value = settings['scroll_flag'])
    gas_l = tk.Label(graph_frame, text = 'Parameter\nto analyse:', font = NORM_FONT)
    gas_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    gas_c = ttk.OptionMenu(graph_frame, gas, *gas_options)
    gas_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Set x-axis for time series
    x_plot_options = [settings['x_plot'],
                      'Time',
                      'Number']
    x_plot = tk.StringVar(graph_frame, value = settings['resid_type'])
    x_plot_l = tk.Label(graph_frame, text = 'Time Series\nx-axix', font = NORM_FONT)
    x_plot_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    x_plot_m = ttk.OptionMenu(graph_frame, x_plot, *x_plot_options)
    x_plot_m.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1 
    
    # Turn on/off graph scrolling
    scroll_b = tk.BooleanVar(graph_frame, value = settings['scroll_flag'])
    scroll_l = tk.Label(graph_frame, text = 'Scroll Graph?', font = NORM_FONT)
    scroll_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    scroll_c = ttk.Checkbutton(graph_frame, variable = scroll_b)
    scroll_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Set number of spectra to display on graph
    popup.spec_no = tk.IntVar(graph_frame, value = settings['scroll_spec_no'])
    spec_no_l = tk.Label(graph_frame, text = 'No. Spectra to display', font = NORM_FONT)
    spec_no_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    spec_no_e = ttk.Entry(graph_frame, textvariable = popup.spec_no, font = NORM_FONT)
    spec_no_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Set format of residual
    resid_options = [settings['resid_type'],
                     'Absolute',
                     'Percentage',
                     'Spec/Fit']
    resid_type = tk.StringVar(graph_frame, value = settings['resid_type'])
    resid_type_l = tk.Label(graph_frame, text = 'Residual Display', font = NORM_FONT)
    resid_type_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    resid_type_m = ttk.OptionMenu(graph_frame, resid_type, *resid_options)
    resid_type_m.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1

#========================================================================================   
#===================================Gas Xsec analyser====================================
#======================================================================================== 

def conv_xsec():
    
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
   
#===================================Filepath controls====================================
    
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
   
#====================================Program controls====================================
    
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

#====================================Convert button======================================
    
    # Create frame to hold the button
    button_frame = ttk.Frame(win)
    button_frame.grid(row=2, column=0, sticky='N')
    
    # Create button to convert the selected spectrum
    conv_b = ttk.Button(button_frame, text = 'Convert', command = convert)
    conv_b.grid()

#=====================================Graph control======================================
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    