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
from tkinter import filedialog as fd
import seabreeze.spectrometers as sb
from seabreeze.cseabreeze.wrapper import SeaBreezeError

from ifit_lib.file_control import make_directory
from ifit_lib.acquire_spectrum import acquire_spectrum
from ifit_lib.update_graph import update_graph


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
                          'Spectrum number ' + str(settings['loop']))
    else:
        self.toggle_button.config(text = 'FITTING ON')
        self.toggle_button.config(bg = 'green')
        self.print_output('Fitting turned on\n' +\
                          'Spectrum number ' + str(settings['loop']))

#========================================================================================
#=================================Connect to Spectrometer================================
#========================================================================================
        
# Function to connect to the attached spectrometer
def connect_spec(self, settings):
    
    # Find connected spectrometers
    devices = sb.list_devices()
    
    # If no devices are connected then set string to show. Else assign first to spec
    if len(devices) == 0:
        settings['spec'] = 0
        settings['Spectrometer'] = 'No devices found'
        devices = ['No devices found']
    else:
        try:
            # Connect to spectrometer
            settings['spec'] = sb.Spectrometer(devices[0])
            
            # Set intial integration time
            settings['spec'].integration_time_micros(float(self.int_time.get())*1000)
            
            # Record serial number in settings
            settings['Spectrometer'] = str(settings['spec'].serial_number)
            
            self.print_output('Spectrometer '+settings['Spectrometer']+' Connected')
            
            # Create filepath to directory to hold program outputs
            results_folder = 'Results/iFit/'+str(datetime.date.today())+'/ifit_out/'
           
            # Create folder
            settings['rt_folder'] = make_directory(results_folder, overwrite = True)
            
            # Create notes file
            settings['notes_fname'] = settings['rt_folder'] + 'notes.txt'
            with open(settings['notes_fname'], 'w') as w:
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
        settings['spec'].integration_time_micros(settings['int_time']*1000)
        
        self.print_output('Integration time updated to '+str(settings['int_time']) +\
                          '\nSpectrum number ' + str(settings['loop']))
        
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
                                               settings['spec'], 
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
        dark_fp = settings['rt_folder'] + 'dark/'
        
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
                spc = acquire_spectrum(self, settings['spec'], settings['int_time'],
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
                          'Spectrum no: ' + str(settings['loop']) + '\n' + \
                          'No. darks: ' + str(dark_n) + '\n' + \
                          'Integration time (ms): '+str(settings['int_time'])+'\n'+\
                          'Coadds: ' + str(self.coadds.get()))
         
        # Set dark_flag to True
        settings['rt_dark_flag'] = True
        
        # Update status
        self.status.set('Standby')
        
#========================================================================================
#===========================================Stop=========================================
#========================================================================================
        
def stop(self, settings):
    settings['stop_flag'] = True
    self.print_output('Loop Stopped\nSpectrum number ' + str(settings['loop']-1))

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
        settings['wave_start']       = float(popup.wave_start.get())
        settings['wave_stop']        = float(popup.wave_stop.get())
        settings['model_resolution'] = float(model_res_e.get())
        settings['dark_flag']        = popup.dark_b.get()
        settings['flat_flag']        = popup.flat_b.get()
        settings['solar_resid_flag'] = popup.resid_b.get()
        settings['update_params']    = popup.update_b.get()
        settings['good_fit_bound']   = float(fit_bound_e.get())
        
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
        settings['calc_shift_flag']  = popup.calc_shift_b.get()
        settings['Fit ILS']          = popup.ils_width_c.get()
        settings['get_ils_flag']     = popup.get_ils_b.get()
        settings['Fit LDF']          = popup.ldf_c.get()
        settings['Fit stretch']      = popup.stretch_c.get()
        settings['Fit ring']         = popup.ring_amt_c.get()
        
        # Update data base settings
        settings['sol_path']         = sol_path_e.get()
        settings['ring_path']        = ring_path_e.get()
        settings['so2_path']         = so2_path_e.get()
        settings['no2_path']         = no2_path_e.get()
        settings['o3_path']          = o3_path_e.get()
        settings['bro_path']         = bro_path_e.get()
        settings['solar_resid_path'] = solar_resid_path_e.get()
        
        # Update graph settings
        settings['Show Graphs']      = graph_b.get()
        settings['scroll_flag']      = scroll_b.get()
        settings['scroll_spec_no']   = int(spec_no_e.get())
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
                
            self.canvas.draw()
        
        # Print message to confirm settings have been updated
        self.print_output('Settings updated')
        
        # Turn on flag to re-build the forward model
        self.build_model_flag = True
        
        # Close the window
        if close_flag == True:
            popup.destroy()
        
    # Make function to select filepaths
    def update_fp(entry):
        
        # Open dialouge to get files
        fpath = fd.askopenfilenames()
        
        entry.set(str(fpath[0]))   

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
    b1.grid(row = 0, column = 0, padx = 90, pady = 5, sticky = "EW")
    
    # Button to just apply changes
    b2 = ttk.Button(button_frame, text='Apply', command=lambda: update_settings(settings,
                                                                                False))
    b2.grid(row = 0, column = 1, padx = 90, pady = 5, sticky = "EW")
    
    # Buttong to cancel
    b3 = ttk.Button(button_frame, text='Cancel', command=lambda: popup.destroy())
    b3.grid(row = 0, column = 2, padx = 90, pady = 5, sticky = "EW")

#========================================================================================
#=====================================Model Settings=====================================
#========================================================================================

#======================================Model Setup=======================================

    # Create row number counter
    row_n = 0
    
    # Create entry for start and stop wavelengths
    fit_range_l = tk.Label(model_frame, text = 'Fit Range:', font = NORM_FONT)
    fit_range_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'E')
    
    popup.wave_start = tk.DoubleVar(model_frame, value = settings['wave_start'])
    wave_start_e = ttk.Entry(model_frame, textvariable = popup.wave_start,
                                  width = 12)
    wave_start_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    popup.wave_stop = tk.DoubleVar(model_frame, value = settings['wave_stop'])
    wave_stop_l = tk.Label(model_frame, text = 'to:', font = NORM_FONT)
    wave_stop_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'E')
    wave_stop_e = ttk.Entry(model_frame, textvariable = popup.wave_stop, width = 12)
    wave_stop_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Set resolution of model grid
    popup.model_res = tk.DoubleVar(model_frame, value = settings['model_res'])
    model_res_l = tk.Label(model_frame, text='Model Grid\nSpacing (nm):', font=NORM_FONT)
    model_res_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'E')
    model_res_e = ttk.Entry(model_frame, textvariable = popup.model_res, width = 12)
    model_res_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove dark spectra
    popup.dark_b = tk.BooleanVar(model_frame, value = settings['dark_flag'])
    dark_l = tk.Label(model_frame, text = 'Remove Dark\nSpectrum?', font = NORM_FONT)
    dark_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'E')
    dark_c = ttk.Checkbutton(model_frame, variable = popup.dark_b)
    dark_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove flat spectra
    popup.flat_b = tk.BooleanVar(model_frame, value = settings['flat_flag'])
    flat_l = tk.Label(model_frame, text = 'Remove Flat\nSpectrum?', font = NORM_FONT)
    flat_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'E')
    flat_c = ttk.Checkbutton(model_frame, variable = popup.flat_b)
    flat_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove, form or ignore the solar residual spectrum
    resid_options = [settings['solar_resid_flag'],
                     'Ignore',
                     'Generate',
                     'Remove']
    popup.resid_b = tk.StringVar(model_frame, value = settings['solar_resid_flag'])
    resid_l = tk.Label(model_frame, text = 'Solar\nresidual:', font = NORM_FONT)
    resid_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'E')
    resid_c = ttk.OptionMenu(model_frame, popup.resid_b, *resid_options)
    resid_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to update fit parameter guesses with the last fit values
    popup.update_b = tk.BooleanVar(model_frame, value = settings['update_params'])
    update_l = tk.Label(model_frame, text='Auto-update\nfit parameters?', font=NORM_FONT)
    update_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'E')
    update_c = ttk.Checkbutton(model_frame, variable = popup.update_b)
    update_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control bound of goodness of fit
    popup.fit_bound = tk.DoubleVar(model_frame, value = settings['good_fit_bound'])
    fit_bound_l = tk.Label(model_frame, text = 'Good Fit\nBound (%):', font = NORM_FONT)
    fit_bound_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'E')
    fit_bound_e = ttk.Entry(model_frame, textvariable = popup.fit_bound, width = 12)
    fit_bound_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Create separator
    sep = ttk.Separator(model_frame, orient='vertical')
    sep.grid(column=2, row = 0, rowspan=8, sticky='ns', padx=10, pady=5)

#=======================================Parameters=======================================
    
    # Reset row counter
    row_n = 0
    
    # Create array of fitting options
    fit_options = ['', 'Fit', 'Fix', 'N/A']
    
    # Gas Amounts
    popup.so2_amt = tk.DoubleVar(model_frame, value = settings['so2_amt'])
    so2_amt_l = tk.Label(model_frame, text = 'SO2:', font = NORM_FONT)
    so2_amt_l.grid(row = row_n, column = 3, padx = 5, pady = 5, sticky = 'W')
    so2_amt_e = ttk.Entry(model_frame, textvariable = popup.so2_amt, width = 12)
    so2_amt_e.grid(row = row_n, column = 4, padx = 5, pady = 5)
    popup.so2_amt_c = tk.StringVar(model_frame, value = settings['Fit so2'])
    fit_options[0] = settings['Fit so2']
    so2_c = ttk.OptionMenu(model_frame, popup.so2_amt_c, *fit_options)
    so2_c.grid(row = row_n, column = 5, padx = 5, pady = 5)
    row_n += 1
    
    popup.no2_amt = tk.DoubleVar(model_frame, value = settings['no2_amt'])
    no2_amt_l = tk.Label(model_frame, text = 'NO2:', font = NORM_FONT)
    no2_amt_l.grid(row = row_n, column = 3, padx = 5, pady = 5, sticky = 'W')
    no2_amt_e = ttk.Entry(model_frame, textvariable = popup.no2_amt, width = 12)
    no2_amt_e.grid(row = row_n, column = 4, padx = 5, pady = 5)
    popup.no2_amt_c = tk.StringVar(model_frame, value = settings['Fit no2'])
    fit_options[0] = settings['Fit no2']
    no2_c = ttk.OptionMenu(model_frame, popup.no2_amt_c, *fit_options)
    no2_c.grid(row = row_n, column = 5, padx = 5, pady = 5)
    row_n += 1
    
    popup.o3_amt = tk.DoubleVar(model_frame, value = settings['o3_amt'])
    o3_amt_l = tk.Label(model_frame, text = 'O3:', font = NORM_FONT)
    o3_amt_l.grid(row = row_n, column = 3, padx = 5, pady = 5, sticky = 'W')
    o3_amt_e = ttk.Entry(model_frame, textvariable = popup.o3_amt, width = 12)
    o3_amt_e.grid(row = row_n, column = 4, padx = 5, pady = 5)
    popup.o3_amt_c = tk.StringVar(model_frame, value = settings['Fit o3'])
    fit_options[0] = settings['Fit o3']
    o3_c = ttk.OptionMenu(model_frame, popup.o3_amt_c, *fit_options)
    o3_c.grid(row = row_n, column = 5, padx = 5, pady = 5)
    row_n += 1
   
    popup.bro_amt = tk.DoubleVar(model_frame, value = settings['bro_amt'])
    bro_amt_l = tk.Label(model_frame, text = 'BrO:', font = NORM_FONT)
    bro_amt_l.grid(row = row_n, column = 3, padx = 5, pady = 5, sticky = 'W')
    bro_amt_e = ttk.Entry(model_frame, textvariable = popup.bro_amt, width = 12)
    bro_amt_e.grid(row = row_n, column = 4, padx = 5, pady = 5)
    popup.bro_amt_c = tk.StringVar(model_frame, value = settings['Fit bro'])
    fit_options[0] = settings['Fit bro']
    bro_c = ttk.OptionMenu(model_frame, popup.bro_amt_c, *fit_options)
    bro_c.grid(row = row_n, column = 5, padx = 5, pady = 5) 
    row_n += 1

    # Spectrometer wavelength shift
    popup.shift = tk.DoubleVar(model_frame, value = settings['shift'])
    shift_l = tk.Label(model_frame, text = 'Shift:', font = NORM_FONT)
    shift_l.grid(row = row_n, column = 3, padx = 5, pady = 5, sticky = 'W')
    shift_e = ttk.Entry(model_frame, textvariable = popup.shift, width = 12)
    shift_e.grid(row = row_n, column = 4, padx = 5, pady = 5)
    popup.shift_c = tk.StringVar(model_frame, value = settings['Fit shift'])
    fit_options[0] = settings['Fit shift']
    shift_c = ttk.OptionMenu(model_frame, popup.shift_c, *fit_options)
    shift_c.grid(row = row_n, column = 5, padx = 5, pady = 5)
    
    # When shift is pre-calculated, deactivate input
    popup.calc_shift_b = tk.BooleanVar(model_frame, value = settings['calc_shift_flag'])
    def shift_activate_check():
        if popup.calc_shift_b.get() == 0:
            shift_e.config(state='normal')
        elif popup.calc_shift_b.get() == 1:
            shift_e.config(state='disabled')
    shift_activate_check()
    
    # Control whether or not to pre-calculate the wavelength shift
    calc_shift_l = tk.Label(model_frame, text = 'Pre-calc\nshift?', font = NORM_FONT)
    calc_shift_l.grid(row = row_n, column = 6, padx = 5, pady = 5, sticky = 'E')
    calc_shift_c = ttk.Checkbutton(model_frame, variable = popup.calc_shift_b, 
                                command = shift_activate_check)
    calc_shift_c.grid(row = row_n, column = 7, padx = 5, pady = 5)
    row_n += 1 
    
    # Spectrometer ILS width
    popup.ils_width = tk.DoubleVar(model_frame, value = settings['ils_width'])
    ils_width_l = tk.Label(model_frame, text = 'ILS\nWidth:', font = NORM_FONT)
    ils_width_l.grid(row = row_n, column = 3, padx = 5, pady = 5, sticky = 'W')
    ils_width_e = ttk.Entry(model_frame, textvariable = popup.ils_width, width = 12)
    ils_width_e.grid(row = row_n, column = 4, padx = 5, pady = 5)
    popup.ils_width_c = tk.StringVar(model_frame, value = settings['Fit ILS'])
    ils_fit_options = [settings['Fit ILS'], 'Fit', 'Fix']
    ils_width_c = ttk.OptionMenu(model_frame, popup.ils_width_c, *ils_fit_options)
    ils_width_c.grid(row = row_n, column = 5, padx = 5, pady = 5)
    
    # When ILS width is read in, deactivate input    
    popup.get_ils_b = tk.BooleanVar(model_frame, value = settings['get_ils_flag'])
    def ils_activate_check():
        if popup.get_ils_b.get() == 0:
            ils_width_e.config(state='normal')
        elif popup.get_ils_b.get() == 1:
            ils_width_e.config(state='disabled')
    ils_activate_check()
    
    # Control whether or not to read the ILS width from the file
    get_ils_l = tk.Label(model_frame, text = 'Use ILS\nfrom file?', font = NORM_FONT)
    get_ils_l.grid(row = row_n, column = 6, padx = 5, pady = 5, sticky = 'W')
    get_ils_c = ttk.Checkbutton(model_frame, variable = popup.get_ils_b, 
                                command = ils_activate_check)
    get_ils_c.grid(row = row_n, column = 7, padx = 5, pady = 5)
    row_n += 1 
        
    # ILS Gaussian weighting
    popup.gauss_weight = tk.DoubleVar(model_frame, value = settings['gauss_weight'])
    gauss_weight_l = tk.Label(model_frame, text = 'ILS Gauss\nWeight:', 
                              font = NORM_FONT)
    gauss_weight_l.grid(row = row_n, column = 3, padx = 5, pady = 5, sticky = 'W')
    gauss_weight_e = ttk.Entry(model_frame, textvariable = popup.gauss_weight,
                               width = 12)
    gauss_weight_e.grid(row = row_n, column = 4, padx = 5, pady = 5)
    row_n += 1 
    
    # Light dilution factor
    popup.ldf = tk.DoubleVar(model_frame, value = settings['ldf'])
    ldf_l = tk.Label(model_frame, text = 'LDF:', font = NORM_FONT)
    ldf_l.grid(row = row_n, column = 3, padx = 5, pady = 5, sticky = 'W')
    ldf_e = ttk.Entry(model_frame, textvariable = popup.ldf, width = 12)
    ldf_e.grid(row = row_n, column = 4, padx = 5, pady = 5)
    popup.ldf_c = tk.StringVar(model_frame, value = settings['Fit LDF'])
    fit_options[0] = settings['Fit LDF']
    ils_width_c = ttk.OptionMenu(model_frame, popup.ldf_c, *fit_options)
    ils_width_c.grid(row = row_n, column = 5, padx = 5, pady = 5)
    row_n += 1 
    
    # Create separator
    sep = ttk.Separator(model_frame, orient='vertical')
    sep.grid(column=6, row = 0, rowspan=4, sticky='ns', padx=10, pady=5)
    
#=======================================New Column=======================================
    
    # Reset row counter
    row_n = 0
    
    # Polynomial coefficents
    popup.poly_n = tk.IntVar(model_frame, value = int(settings['poly_n']))
    poly_n_vals = [0,1,2,3,4,5,6,7,8,9,10]
    poly_n_l = tk.Label(model_frame, text = 'Poly\nOrder:', font = NORM_FONT)
    poly_n_l.grid(row = row_n, column = 7, padx = 5, pady = 5, sticky = 'W')
    poly_n_e = tk.Spinbox(model_frame, values = poly_n_vals, width = 12,
                          textvariable = popup.poly_n)
    poly_n_e.grid(row = row_n, column = 8, padx = 5, pady = 5)
    popup.poly_n.set(settings['poly_n'])
    row_n += 1
    
    # Spectrometer wavelength stretch
    popup.stretch = tk.DoubleVar(model_frame, value = settings['stretch'])
    stretch_l = tk.Label(model_frame, text = 'Stretch:', font = NORM_FONT)
    stretch_l.grid(row = row_n, column = 7, padx = 5, pady = 5, sticky = 'W')
    stretch_e = ttk.Entry(model_frame, textvariable = popup.stretch, width = 12)
    stretch_e.grid(row = row_n, column = 8, padx = 5, pady = 5)
    popup.stretch_c = tk.StringVar(model_frame, value = settings['Fit stretch'])
    fit_options[0] = settings['Fit stretch']
    stretch_c = ttk.OptionMenu(model_frame, popup.stretch_c, *fit_options)
    stretch_c.grid(row = row_n, column = 9, padx = 5, pady = 5)
    row_n += 1
    
    # Ring effect
    popup.ring_amt = tk.DoubleVar(model_frame, value = settings['ring_amt'])
    ring_amt_l = tk.Label(model_frame, text = 'Ring:', font = NORM_FONT)
    ring_amt_l.grid(row = row_n, column = 7, padx = 5, pady = 5, sticky = 'W')
    ring_amt_e = ttk.Entry(model_frame, textvariable = popup.ring_amt, width = 12)
    ring_amt_e.grid(row = row_n, column = 8, padx = 5, pady = 5)
    popup.ring_amt_c = tk.StringVar(model_frame, value = settings['Fit ring'])
    fit_options[0] = settings['Fit ring']
    ring_c = ttk.OptionMenu(model_frame, popup.ring_amt_c, *fit_options)
    ring_c.grid(row = row_n, column = 9, padx = 5, pady = 5)
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
    sol_path_b = ttk.Button(datab_frame, text = "Update", 
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
    ring_path_b = ttk.Button(datab_frame, text = "Update", 
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
    so2_path_b = ttk.Button(datab_frame, text = "Update", 
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
    no2_path_b = ttk.Button(datab_frame, text = "Update", 
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
    o3_path_b = ttk.Button(datab_frame, text = "Update", 
                           command = lambda: update_fp(popup.o3_path))
    o3_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # Bro xsec
    popup.bro_path = tk.StringVar(datab_frame, value = settings['bro_path'])
    bro_path_l = tk.Label(datab_frame, text = 'BrO xsec:', font = NORM_FONT)
    bro_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    bro_path_e = ttk.Entry(datab_frame, textvariable = popup.bro_path, font = NORM_FONT,
                           width = 40)
    bro_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5) 
    bro_path_b = ttk.Button(datab_frame, text = "Update", 
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
    solar_resid_path_b = ttk.Button(datab_frame, text = "Update", 
                                    command = lambda: update_fp(popup.solar_resid_path))
    solar_resid_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1    
    
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
                   'O3',
                   'BrO']
    gas = tk.StringVar(graph_frame, value = settings['scroll_flag'])
    gas_l = tk.Label(graph_frame, text = 'Gas to analyse:', font = NORM_FONT)
    gas_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    gas_c = ttk.OptionMenu(graph_frame, gas, *gas_options)
    gas_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
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