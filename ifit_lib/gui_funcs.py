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
        settings['int_time'] = float(self.int_time.get())*1000
        settings['spec'].integration_time_micros(float(settings['int_time']))
        
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
        
        # Update model settings in common
        settings['model_resolution'] = float(model_res_e.get())
        settings['dark_flag']        = dark_b.get()
        settings['flat_flag']        = flat_b.get()
        settings['update_params']    = update_b.get()
        settings['solar_resid_flag'] = resid_b.get()
        settings['good_fit_bound']   = float(fit_bound_e.get())
        
        # Update data base settings in common
        settings['sol_path']         = sol_path_e.get()
        settings['ring_path']        = ring_path_e.get()
        settings['so2_path']         = so2_path_e.get()
        settings['no2_path']         = no2_path_e.get()
        settings['o3_path']          = o3_path_e.get()
        settings['bro_path']         = bro_path_e.get()
        settings['solar_resid_path'] = solar_resid_path_e.get()
        
        # Update graph settings in common
        settings['Show Graphs']      = graph_b.get()
        #settings['Show Error Bars']  = err_b.get()
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
    
    nb.add(model_frame, text = 'Model Settings')
    nb.add(graph_frame, text = 'Graph Settings')
    nb.add(datab_frame, text = 'Data Base Settings')
    
    nb.grid(row = 0, column=0, padx=10, pady=10, sticky = 'NW', columnspan = 10)
    '''
    # Create frames for diferent sections for post analysis and buttons
    model_frame = tk.LabelFrame(popup, text = 'Model Settings', font = LARG_FONT)
    model_frame.grid(row=0, column=0, padx = 10, pady = 10, sticky="NW")
    
    datab_frame = tk.LabelFrame(popup, text='Data Base Settings', font = LARG_FONT)
    datab_frame.grid(row=0, column=1, padx = 10, pady = 10, sticky="NW", rowspan=2)
    
    graph_frame = tk.LabelFrame(popup, text='Graph Settings', font = LARG_FONT)
    graph_frame.grid(row=1, column=0, padx = 10, pady = 10, sticky="NW")
    '''
    button_frame = ttk.Frame(popup)
    button_frame.grid(row=1, column=0, padx = 10, pady = 10, sticky="EW")

    # Button to apply changes and close
    b1 = ttk.Button(button_frame, text='Ok', command=lambda: update_settings(settings, 
                                                                             True))
    b1.grid(row = 0, column = 0, padx = 40, pady = 5, sticky = "EW")
    
    # Button to just apply changes
    b2 = ttk.Button(button_frame, text='Apply', command=lambda: update_settings(settings,
                                                                                False))
    b2.grid(row = 0, column = 1, padx = 40, pady = 5, sticky = "EW")
    
    # Buttong to cancel
    b3 = ttk.Button(button_frame, text='Cancel', command=lambda: popup.destroy())
    b3.grid(row = 0, column = 2, padx = 40, pady = 5, sticky = "EW")

#========================================================================================
#=====================================Model Settings=====================================
#========================================================================================

    # Create row number counter
    row_n = 0
    
    # Set resolution of model grid
    popup.model_res = tk.DoubleVar(model_frame, value = settings['model_res'])
    model_res_l = tk.Label(model_frame, text = 'Model Grid spacing:', font = NORM_FONT)
    model_res_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    model_res_e = ttk.Entry(model_frame, textvariable = popup.model_res)
    model_res_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove dark spectra
    dark_b = tk.BooleanVar(model_frame, value = settings['dark_flag'])
    dark_l = tk.Label(model_frame, text = 'Remove Dark Spectra?', font = NORM_FONT)
    dark_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    dark_c = ttk.Checkbutton(model_frame, variable = dark_b)
    dark_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove flat spectra
    flat_b = tk.BooleanVar(model_frame, value = settings['flat_flag'])
    flat_l = tk.Label(model_frame, text = 'Remove Flat Spectra?', font = NORM_FONT)
    flat_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    flat_c = ttk.Checkbutton(model_frame, variable = flat_b)
    flat_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove, form or ignore the solar residual spectrum
    resid_options = [settings['solar_resid_flag'],
                     'Ignore',
                     'Generate',
                     'Remove']
    resid_b = tk.StringVar(model_frame, value = settings['solar_resid_flag'])
    resid_l = tk.Label(model_frame, text = 'Solar residual:', font = NORM_FONT)
    resid_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    resid_c = ttk.OptionMenu(model_frame, resid_b, *resid_options)
    resid_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to update fit parameter guesses with the last fit values
    update_b = tk.BooleanVar(model_frame, value = settings['update_params'])
    update_l = tk.Label(model_frame, text = 'Auto-update fit parameters?', font = NORM_FONT)
    update_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    update_c = ttk.Checkbutton(model_frame, variable = update_b)
    update_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control bound of goodness of fit
    popup.fit_bound = tk.DoubleVar(model_frame, value = settings['good_fit_bound'])
    fit_bound_l = tk.Label(model_frame, text = 'Good Fit Bound (%):', font = NORM_FONT)
    fit_bound_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    fit_bound_e = ttk.Entry(model_frame, textvariable = popup.fit_bound)
    fit_bound_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
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
    graph_b = tk.IntVar(graph_frame, value = settings['Show Graphs'])
    graph_l = tk.Label(graph_frame, text = 'Show Graphs?', font = NORM_FONT)
    graph_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    graph_c = ttk.Checkbutton(graph_frame, variable = graph_b)
    graph_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    '''
    # Turn on/off error bars
    err_b = tk.IntVar(graph_frame, value = settings['Show Error Bars'])
    err_l = tk.Label(graph_frame, text = 'Show Error Bars?', font = NORM_FONT)
    err_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    err_c = ttk.Checkbutton(graph_frame, variable = err_b)
    err_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    '''
    
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
    scroll_b = tk.IntVar(graph_frame, value = settings['scroll_flag'])
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