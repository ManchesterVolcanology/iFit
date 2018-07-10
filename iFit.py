# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 09:24:05 2018

@author: mqbpwbe2
"""

# Import required libraries
import matplotlib
matplotlib.use('TkAgg')
import os
import traceback
import tkinter.messagebox as tkMessageBox
import tkinter.scrolledtext as tkst
import glob
import numpy as np
from tkinter import ttk
import tkinter as tk
from queue import Queue
from threading import Thread
from seabreeze.cseabreeze.wrapper import SeaBreezeError
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import OrderedDict
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

from ifit_lib.build_fwd_data import build_fwd_data
from ifit_lib.read_spectrum import read_spectrum, average_spectra
from ifit_lib.fit import fit_spec, ifit_fwd
from ifit_lib.acquire_spectrum import acquire_spectrum
from ifit_lib.find_nearest import extract_window
from ifit_lib.update_graph import update_graph
from ifit_lib.file_control import make_csv_file
from ifit_lib.gui_funcs import adv_settings, fit_toggle, spec_fp, dark_fp, stop, \
                               connect_spec, update_int_time, test_spec, read_darks

# Define some fonts to use in the program
NORM_FONT = ('Verdana', 8)
MED_FONT  = ('Veranda', 11)
LARG_FONT = ('Verdana', 12, 'bold')

class mygui(tk.Tk):
    
    def __init__(self, *args, **kwargs):
        
#========================================================================================
#=================================Build GUI Containers===================================
#========================================================================================
        
        # Create GUI in the backend
        tk.Tk.__init__(self, *args, **kwargs)
        
        # Cause exceptions to report in a new window
        tk.Tk.report_callback_exception = self.report_callback_exception
               
        # Close program on closure of window
        self.protocol("WM_DELETE_WINDOW", self.handler)
        
        # Button Style
        ttk.Style().configure('TButton', width = 20, height = 20, relief="flat") 
        
        # Add a title and icon
        tk.Tk.wm_title(self, 'iFit-2-2')
        tk.Tk.iconbitmap(self, default = 'data_bases/icon.ico')
        
        # Create notebook to hold different frames
        nb = ttk.Notebook(self)
        page1 = ttk.Frame(nb)
        page2 = ttk.Frame(nb)
        page3 = ttk.Frame(nb)
        
        # Create two frames, one for post analysis and one for real time acquisition
        nb.add(page2, text = 'Real Time Acquisition')
        nb.add(page1, text = 'Post Analysis')
        nb.add(page3, text = 'Forward Model')
        
        nb.grid(column=0, padx=10, pady=10, sticky = 'NW')
        
        # Create frame to hold graphs
        graph_frame = ttk.Frame(self, relief = 'groove')
        graph_frame.grid(row=0, column=1, padx=10, pady=10, rowspan=10, sticky="NW")
        graph_frame.columnconfigure(index=0, weight=1)
        graph_frame.rowconfigure(index = 0, weight = 1)
        
        # Create frame to hold text output
        text_frame = ttk.Frame(self, relief = 'groove')
        text_frame.grid(row=2, column=0, padx=10, pady=10, rowspan=10, sticky="NW")
        
        # Frame for quick analysis
        quick_frame = tk.LabelFrame(self, text = 'Quick Analysis', font = LARG_FONT)
        quick_frame.grid(row=1, column=0, padx=10, pady=10, sticky="NW")
        
        mygui.columnconfigure(index=1, weight=1, self = self)
        mygui.rowconfigure(index = 5, weight = 1, self = self)
        
#========================================================================================
#====================================Create text output==================================
#========================================================================================        
                 
        # Build text box
        self.text_box = tkst.ScrolledText(text_frame, width = 55, height = 10)
        self.text_box.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W',
                           columnspan = 4)
        self.text_box.insert('1.0', 'Welcome to iFit! Written by Ben Esse\n\n')  
        
        # Create button for advanced settings
        adv_set_b = ttk.Button(text_frame, text = 'Adv. Settings', 
                            command = lambda: adv_settings(self, settings, 'iFit'))
        adv_set_b.grid(row = 0, column = 0, padx = 5, pady = 5)

        # Create button to save settings
        save_b = ttk.Button(text_frame, text = 'Save Settings', command = self.save)
        save_b.grid(row = 0, column = 1, padx = 5, pady = 5)
        
        # Create a button to exit
        exit_b = ttk.Button(text_frame, text = 'Exit', command = self.quit)
        exit_b.grid(row = 0, column = 2, padx = 5, pady = 5)
        
#========================================================================================
#==============================Create quick analysis outputs=============================
#========================================================================================
        
        # Create progress bar
        self.progress = ttk.Progressbar(quick_frame, orient = tk.HORIZONTAL, length=350,
                                        mode = 'determinate')
        self.progress.grid(row = 0, column = 1, padx = 5, pady = 5, columnspan = 3)
        
        # Create status indicator
        self.status = tk.StringVar(quick_frame, value = 'Standby')
        self.status_e = tk.Label(quick_frame, textvariable = self.status)
        self.status_e.grid(row=0, column=0, padx=5, pady=5, sticky="EW")
        
        # Create ouput for last so2 amount
        self.last_so2_amt = tk.StringVar(self, value = '-')
        last_so2_amt_l = tk.Label(quick_frame, text = 'Last amt:', font = NORM_FONT)
        last_so2_amt_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        last_so2_amt_e = tk.Label(quick_frame, textvariable = self.last_so2_amt)
        last_so2_amt_e.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        # Create ouput for last so2 error
        self.last_so2_err = tk.StringVar(self, value = '-')
        last_so2_err_l = tk.Label(quick_frame, text = 'Last error:', 
                                  font = NORM_FONT)
        last_so2_err_l.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
        last_so2_err_e = tk.Label(quick_frame, textvariable = self.last_so2_err)
        last_so2_err_e.grid(row = 1, column = 3, padx = 5, pady = 5, sticky = 'W')
        
#========================================================================================
#===================================Set program settings=================================
#========================================================================================

        # Create settings dictionary
        global settings
        settings = {}
        
        # Read in settings file
        try:
            with open('data_bases/ifit_settings.txt', 'r') as r:
                
                # Read data
                data = r.readlines()
                
                # Unpack and save to dictionary
                for i in data:
                    name, val = i.strip().split(';')
                    settings[name] = val
   
        except FileNotFoundError:
            self.print_output('No settings file found, reverting to origional')
            settings['Wave Start']        = 305
            settings['Wave Stop']         = 318
            settings['Spectrometer']      = '-select-'
            settings['Spectra Type']      = '-select-'
            settings['int_time']          = 100
            settings['coadds']            = 10
            settings['no_darks']          = 10
            settings['ILS Width']         = 0.52
            settings['Gauss Weight']      = 1.0
            settings['Fit ILS']           = 'Fix'
            settings['LDF']               = 0.0
            settings['Fit LDF']           = 'N/A'
            settings['dark_flag']         = 1
            settings['flat_flag']         = 1
            settings['update_params']     = 1
            settings['good_fit_bound']    = 10
            settings['Show Graphs']       = 1
            settings['Show Error Bars']   = 0
            settings['analysis_gas']      = 'SO2'
            settings['scroll_flag']       = 1
            settings['scroll_spec_no']    = 200
            settings['resid_type']        = 'Spec/Fit'
            settings['solar_resid_flag']  = 'Ignore'
            settings['poly_n']            = 3
            settings['shift']             = -0.2
            settings['stretch']           = 0.05
            settings['ring']              = 1.0
            settings['SO2']               = 1e+16
            settings['NO2']               = 1e+17
            settings['O3']                = 1e+19
            settings['BrO']               = 1e+15
            settings['Fit shift']         = 'Fit'
            settings['Fit stretch']       = 'Fit'
            settings['Fit ring']          = 'Fit'
            settings['Fit so2']           = 'Fit'
            settings['Fit no2']           = 'Fit'
            settings['Fit o3']            = 'Fit'
            settings['Fit bro']           = 'Fit'
            settings['model_resolution']  = 0.02
            settings['sol_path']          = 'data_bases/gas data/sao2010.txt'
            settings['ring_path']         = 'data_bases/gas data/ring.dat'
            settings['so2_path']          = 'data_bases/gas data/SO2_293K.dat'
            settings['no2_path']          = 'data_bases/gas data/No2_223l.dat'
            settings['o3_path']           = 'data_bases/gas data/O3_293K.dat'
            settings['bro_path']          = 'data_bases/gas data/BrO_Cross_298K.txt'
            settings['solar_resid_path']  = 'data_bases/gas data/solar_resid.txt'
            settings['Spectra Filepaths'] = ''
            settings['Dark Filepaths']    = ''
 
        # Create loop counter to keep track of the analysis
        settings['loop'] = 0
        
        # Create flag to ensure only one output file is created per program launch
        settings['create_out_flag'] = True
        
        # Create flag to see if darks have been measured yet
        settings['rt_dark_flag'] = False
        
#========================================================================================
#====================================Create plot canvas==================================
#========================================================================================        
                    
        # Create figure to hold the graphs
        plt.rcParams.update({'font.size': 8} )
        self.fig = plt.figure(figsize = (8,6))
        gs = gridspec.GridSpec(3,2)
        
        # Create plot axes
        self.ax0 = self.fig.add_subplot(gs[0,0])
        self.ax1 = self.fig.add_subplot(gs[0,1])
        self.ax2 = self.fig.add_subplot(gs[1,0])
        self.ax3 = self.fig.add_subplot(gs[1,1])
        self.ax4 = self.fig.add_subplot(gs[2,:])
        
        # Axes: 1) Spectrum and fit
        #       2) Residual
        #       3) Full spectrum
        #       4) SO2 amount time series
        
        
        # Set axis labels
        self.ax0.set_ylabel('Intensity (arb)', fontsize=10)
        
        self.ax1.set_ylabel('Intensity (arb)', fontsize=10)
        
        if settings['resid_type'] == 'Percentage':
            self.ax2.set_ylabel('Fit residual (%)', fontsize=10)
        if settings['resid_type'] == 'Absolute':
            self.ax2.set_ylabel('Fit residual (Abs)', fontsize=10)
        if settings['resid_type'] == 'Spec/Fit':
            self.ax2.set_ylabel('Fit residual (Spec/Fit)', fontsize=10)
        self.ax2.set_xlabel('Wavelength (nm)', fontsize=10)
        
        self.ax3.set_ylabel(settings['analysis_gas'] + ' Absorbance', fontsize = 10)
        self.ax3.set_xlabel('Wavelength (nm)', fontsize=10)
        
        self.ax4.set_ylabel(settings['analysis_gas'] + ' amt (ppm.m)', fontsize = 10)
        self.ax4.set_xlabel('Spectrum number', fontsize=10)
        
        # Create lines to plot data series
        
        # Spectral data
        self.line0, = self.ax0.plot(0, 0, 'b', label = 'Spectrum')
        self.line1, = self.ax0.plot(0, 0, 'darkorange', label = 'Fit')
        self.ax0.legend(loc = 0)
        
        # Full spectrum
        self.line2, = self.ax1.plot(0, 0, 'b')
        
        # Residual
        self.line3, = self.ax2.plot(0, 0, 'r')
        
        # SO2 transmittance data
        self.line4, = self.ax3.plot(0, 0, 'b', label = 'Spec / F_no_gas')
        self.line5, = self.ax3.plot(0, 0, 'r', label = 'Gas abs')
        self.ax3.legend(loc = 0)
        
        # SO2 Time series and error bars
        self.line6, = self.ax4.plot(0, 0, 'g')
        
        # Make it look nice
        plt.tight_layout()
        
        # Create the canvas to hold the graph in the GUI
        self.canvas = FigureCanvasTkAgg(self.fig, graph_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady = 10)
        
        # Add matplotlib toolbar above the plot canvas
        toolbar_frame = tk.Frame(graph_frame, bg = 'black')  
        toolbar_frame.grid(row=1,column=0, sticky = 'W', padx = 5, pady = 5)                             
        toolbar = NavigationToolbar2TkAgg(self.canvas, toolbar_frame)
        toolbar.update()
        



          



#========================================================================================         
#========================================================================================
#===================================Post Analysis Frame==================================
#======================================================================================== 
#======================================================================================== 





#========================================================================================
#====================================Create GUI frames===================================
#========================================================================================

        # Create frames for diferent sections for post analysis     
        setup_frame = tk.LabelFrame(page1, text = 'Program Setup', font = LARG_FONT)
        setup_frame.grid(row=0, column=0, padx = 10, pady = 10, sticky="ew")
        
        button_frame = tk.LabelFrame(page1, text='Control', font = LARG_FONT)
        button_frame.grid(row=1, column=0, padx = 10, pady = 10, sticky="ew")
      
#========================================================================================
#==================================Create control inputs=================================
#======================================================================================== 
               
        # Find available flat spectra and form into list
        options = [settings['Spectrometer']]       

        for i, name in enumerate(glob.glob('data_bases/Spectrometer/flat_*')):
            options.append(name[29:-4])
        
        # Create entry to select spectrometer
        self.spec_name = tk.StringVar(setup_frame, value = options[0])
        spectro_l = tk.Label(setup_frame, text = 'Spectrometer:', font = NORM_FONT)
        spectro_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        spectro_c = ttk.OptionMenu(setup_frame, self.spec_name, *options)
        spectro_c.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        # Create entry to select spectra type
        spec_options = [settings['Spectra Type'],
                        'iFit',
                        'Jai Spec',
                        'Spectrasuite',
                        'GSJ']
        
        self.spec_type = tk.StringVar(setup_frame, value = spec_options[0])
        spec_l = tk.Label(setup_frame, text = 'Spectra Type:', font = NORM_FONT)
        spec_l.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = 'W')
        spec_c = ttk.OptionMenu(setup_frame, self.spec_type, *spec_options)
        spec_c.grid(row = 0, column = 3, padx = 5, pady = 5, sticky = 'W')
        
#========================================================================================
#==================================Create file dialouges=================================
#========================================================================================        
        
        # Reformat filepath strings
        self.spec_fpaths = []
        for i in settings['Spectra Filepaths'][1:-1].split(', '):
            self.spec_fpaths.append(str(i[1:-1]))
        self.dark_fpaths = []
        for i in settings['Dark Filepaths'][1:-1].split(', '):
            self.dark_fpaths.append(str(i[1:-1]))
        
        # File dialouge for spectra
        if self.spec_fpaths == ['']:
            message = 'No spectra selected'
        else:
            message = str(len(self.spec_fpaths))+' spectra selected'
        self.spec_ent = tk.StringVar(value = message)
        self.specfp_l = tk.Entry(setup_frame, font = NORM_FONT, width = 40, 
                                 text = self.spec_ent)
        self.specfp_l.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = 'W', 
                      columnspan = 3)
        specfp_b = ttk.Button(setup_frame, text="Select Spectra", 
                              command = lambda: spec_fp(self))
        specfp_b.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        
        # File dialouge for darks
        if self.dark_fpaths == ['']:
            message = 'No spectra selected'
        else:
            message = str(len(self.dark_fpaths))+' spectra selected'
        self.dark_ent = tk.StringVar(value = message)
        self.darkfp_l = tk.Entry(setup_frame, font = NORM_FONT, width = 40, 
                                 text = self.dark_ent)
        self.darkfp_l.grid(row = 2, column = 1, padx = 5, pady = 5, sticky = 'W', 
                      columnspan = 3)
        darkfp_b = ttk.Button(setup_frame, text = "Select Darks",
                              command = lambda: dark_fp(self))
        darkfp_b.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
               
#========================================================================================
#===================================Create start button==================================
#========================================================================================         
        
        # Create button to start
        start_b = ttk.Button(button_frame, text = 'Begin!', 
                             command = lambda: self.begin('post_analysis'))
        start_b.grid(row = 0, column = 0, padx = 40, pady = 5, columnspan = 2)
        
        # Create button to stop
        stop_b = ttk.Button(button_frame, text = 'Stop', 
                            command = lambda: stop(self, settings))
        stop_b.grid(row = 0, column = 2, padx = 40, pady = 5, columnspan = 2)
        

        
        
     

#========================================================================================         
#========================================================================================
#================================Real Time Analysis Frame================================
#======================================================================================== 
#======================================================================================== 





#========================================================================================
#====================================Create GUI frames===================================
#========================================================================================

        # Create frames for diferent sections for post analysis     
        setup_frame2 = tk.LabelFrame(page2, text='Spectrometer Setup', font = LARG_FONT)
        setup_frame2.grid(row=0, column=0, padx = 10, pady = 10, sticky="ew")
        
        button_frame2 = tk.LabelFrame(page2, text='Control', font = LARG_FONT)
        button_frame2.grid(row=1, column=0, padx = 10, pady = 10, sticky="ew")

#========================================================================================
#==================================Create control inputs=================================
#========================================================================================
        
        # Create label to display the spectrometer name
        self.c_spec = tk.StringVar(setup_frame2, value = 'None connected')
        c_spec_l = ttk.Label(setup_frame2, text="Device: ", font = NORM_FONT)
        c_spec_l.grid(row=0, column=0, pady=5, padx=5, sticky='W')
        c_spec_e = ttk.Entry(setup_frame2, textvariable = self.c_spec)
        c_spec_e.grid(row = 0, column = 1, padx = 5, pady = 5)
        
        # Integration Time
        self.int_time = tk.DoubleVar(self, value = settings['int_time'])
        int_time_l = tk.Label(setup_frame2, text = 'Integration time (ms):', 
                              font = NORM_FONT)
        int_time_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        int_time_e = ttk.Entry(setup_frame2, textvariable = self.int_time)
        int_time_e.grid(row = 1, column = 1, padx = 5, pady = 5)
        
        # Coadds
        self.coadds = tk.DoubleVar(self, value = settings['coadds'])
        coadds_l = tk.Label(setup_frame2, text = 'Coadds:', 
                              font = NORM_FONT)
        coadds_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        coadds_e = ttk.Entry(setup_frame2, textvariable = self.coadds)
        coadds_e.grid(row = 2, column = 1, padx = 5, pady = 5)
        
        # Number of darks to get
        self.no_darks = tk.DoubleVar(self, value = settings['no_darks'])
        no_darks_l = tk.Label(setup_frame2, text = 'No. Darks:', 
                              font = NORM_FONT)
        no_darks_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
        no_darks_e = ttk.Entry(setup_frame2, textvariable = self.no_darks)
        no_darks_e.grid(row = 3, column = 1, padx = 5, pady = 5)
        
        # Create button to connect to spectrometer
        connect_spec_b = ttk.Button(setup_frame2, text = 'Connect',
                                    command = lambda: connect_spec(self, settings))
        connect_spec_b.grid(row = 0, column = 2, padx = 5, pady = 5)
        
        # Create button to update integration time
        update_int_time_b = ttk.Button(setup_frame2, text = 'Update',
                                       command = lambda: update_int_time(self, settings))
        update_int_time_b.grid(row = 1, column = 2, padx = 5, pady = 5)
        
        # Create button to read a single spectrum
        test_spec_b = ttk.Button(setup_frame2, text = 'Test Spectrum',
                                 command = lambda: test_spec(self, settings, mygui, 
                                                               self.line2, self.ax1))
        test_spec_b.grid(row = 2, column = 2, padx = 5, pady = 5)
        
        # Create button to read darks
        read_darks_b = ttk.Button(setup_frame2, text = 'Acquire Darks',
                                  command = lambda: read_darks(self, settings, mygui, 
                                                               self.line2, self.ax1))
        read_darks_b.grid(row = 3, column = 2, padx = 5, pady = 5)
        
#========================================================================================
#==============================Create start and exit buttons=============================
#========================================================================================         
        
        # Create button to start
        start_aq_b = ttk.Button(button_frame2, text = 'Begin!', 
                                command = lambda: self.begin('rt_analysis'))
        start_aq_b.grid(row = 0, column = 0, padx = 40, pady = 5)
        
        # Create button to stop
        stop_aq_b = ttk.Button(button_frame2, text = 'Stop', 
                               command = lambda: stop(self, settings))
        stop_aq_b.grid(row = 0, column = 1, padx = 40, pady = 5)
        
        # Create switch to toggle fitting on or off
        self.toggle_button = tk.Button(button_frame2, text = 'FITTING OFF', width = 12, 
                                       height = 1, bg = 'red', font = LARG_FONT,
                                       command = lambda: fit_toggle(self, settings))
        self.toggle_button.grid(row=1, column=0, padx=5, pady=5, columnspan=2)
        
        
        
        



#========================================================================================         
#========================================================================================
#=================================Model Parameters Frame=================================
#======================================================================================== 
#======================================================================================== 






#========================================================================================
#====================================Create GUI frames===================================
#========================================================================================

        # Frame for parameters
        model_frame = tk.LabelFrame(page3, text = 'Model Settings', font = LARG_FONT)
        model_frame.grid(row=0, column=0, padx=10, pady=5, sticky="NW")
        
        # Frame for parameters
        param_frame = tk.LabelFrame(page3, text = 'Fit Paramters', font = LARG_FONT)
        param_frame.grid(row=1, column=0, padx=10, pady=5, sticky="NW")

#========================================================================================
#=======================================Model Setup======================================
#========================================================================================

        # Create entry for start and stop wavelengths
        self.wave_start = tk.DoubleVar(model_frame, value = settings['Wave Start'])
        self.wave_start_l = tk.Label(model_frame, text = 'Wave Start:', font = NORM_FONT)
        self.wave_start_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.wave_start_e = ttk.Entry(model_frame, textvariable = self.wave_start,
                                      width = 12)
        self.wave_start_e.grid(row = 0, column = 1, padx = 5, pady = 5)
        
        self.wave_stop = tk.DoubleVar(model_frame, value = settings['Wave Stop'])
        self.wave_stop_l = tk.Label(model_frame, text = 'Wave Stop:', font = NORM_FONT)
        self.wave_stop_l.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.wave_stop_e = ttk.Entry(model_frame, textvariable = self.wave_stop, 
                                     width = 12)
        self.wave_stop_e.grid(row = 0, column = 3, padx = 5, pady = 5)
        
        # Create array of fitting options
        fit_options = ['', 'Fit', 'Fix', 'N/A']
        
        # Instrument lineshape width
        self.ils_width = tk.DoubleVar(model_frame, value = settings['ILS Width'])
        ils_width_l = tk.Label(model_frame, text = 'ILS Width:', font = NORM_FONT)
        ils_width_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        ils_width_e = ttk.Entry(model_frame, textvariable = self.ils_width, width = 12)
        ils_width_e.grid(row = 1, column = 1, padx = 5, pady = 5)
        self.ils_width_c = tk.StringVar(model_frame, value = settings['Fit ILS'])
        ils_fit_options = [settings['Fit ILS'], 'Fit', 'Fix']
        ils_width_c = ttk.OptionMenu(model_frame, self.ils_width_c, *ils_fit_options)
        ils_width_c.grid(row = 1, column = 2, padx = 5, pady = 5)
        
        # ILS Gaussian weighting
        self.gauss_weight = tk.DoubleVar(model_frame, value = settings['Gauss Weight'])
        gauss_weight_l = tk.Label(model_frame, text = 'ILS Gauss\nWeight:', 
                                  font = NORM_FONT)
        gauss_weight_l.grid(row = 1, column = 3, padx = 5, pady = 5, sticky = 'W')
        gauss_weight_e = ttk.Entry(model_frame, textvariable = self.gauss_weight,
                                   width = 12)
        gauss_weight_e.grid(row = 1, column = 4, padx = 5, pady = 5)
        
        # Light dilution factor
        self.ldf = tk.DoubleVar(model_frame, value = settings['LDF'])
        ldf_l = tk.Label(model_frame, text = 'LDF:', font = NORM_FONT)
        ldf_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        ldf_e = ttk.Entry(model_frame, textvariable = self.ldf, width = 12)
        ldf_e.grid(row = 2, column = 1, padx = 5, pady = 5)
        self.ldf_c = tk.StringVar(model_frame, value = settings['Fit LDF'])
        fit_options[0] = settings['Fit LDF']
        ils_width_c = ttk.OptionMenu(model_frame, self.ldf_c, *fit_options)
        ils_width_c.grid(row = 2, column = 2, padx = 5, pady = 5)

#========================================================================================
#================================Set initial fit parameters==============================
#========================================================================================
        
        # Polynomial coefficents
        self.poly_n = tk.DoubleVar(self, value = settings['poly_n'])
        poly_n_vals = [1,2,3,4,5,6,7,8,9,10]
        poly_n_l = tk.Label(param_frame, text = 'Poly Order:', font = NORM_FONT)
        poly_n_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        poly_n_e = tk.Spinbox(param_frame, values = poly_n_vals, width = 12,
                              textvariable = self.poly_n)
        poly_n_e.grid(row = 1, column = 1, padx = 5, pady = 5)
        self.poly_n.set(settings['poly_n'])

        # Spectrometer wavelength shift parameters
        self.shift = tk.DoubleVar(param_frame, value = settings['shift'])
        shift_l = tk.Label(param_frame, text = 'Shift:', font = NORM_FONT)
        shift_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        shift_e = ttk.Entry(param_frame, textvariable = self.shift, width = 12)
        shift_e.grid(row = 2, column = 1, padx = 5, pady = 5)
        self.shift_c = tk.StringVar(param_frame, value = settings['Fit shift'])
        fit_options[0] = settings['Fit shift']
        shift_c = ttk.OptionMenu(param_frame, self.shift_c, *fit_options)
        shift_c.grid(row = 2, column = 2, padx = 5, pady = 5)
        
        self.stretch = tk.DoubleVar(param_frame, value = settings['stretch'])
        stretch_l = tk.Label(param_frame, text = 'Stretch:', font = NORM_FONT)
        stretch_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
        stretch_e = ttk.Entry(param_frame, textvariable = self.stretch, width = 12)
        stretch_e.grid(row = 3, column = 1, padx = 5, pady = 5)
        self.stretch_c = tk.StringVar(param_frame, value = settings['Fit stretch'])
        fit_options[0] = settings['Fit stretch']
        stretch_c = ttk.OptionMenu(param_frame, self.stretch_c, *fit_options)
        stretch_c.grid(row = 3, column = 2, padx = 5, pady = 5)
        
        # Ring effect
        self.ring_amt = tk.DoubleVar(param_frame, value = settings['ring'])
        ring_amt_l = tk.Label(param_frame, text = 'Ring:', font = NORM_FONT)
        ring_amt_l.grid(row = 4, column = 0, padx = 5, pady = 5, sticky = 'W')
        ring_amt_e = ttk.Entry(param_frame, textvariable = self.ring_amt, width = 12)
        ring_amt_e.grid(row = 4, column = 1, padx = 5, pady = 5)
        self.ring_c = tk.StringVar(param_frame, value = settings['Fit ring'])
        fit_options[0] = settings['Fit ring']
        ring_c = ttk.OptionMenu(param_frame, self.ring_c, *fit_options)
        ring_c.grid(row = 4, column = 2, padx = 5, pady = 5)
        
        # Gas amounts
        self.so2_amt = tk.DoubleVar(param_frame, value = settings['SO2'])
        so2_amt_l = tk.Label(param_frame, text = 'SO2:', font = NORM_FONT)
        so2_amt_l.grid(row = 1, column = 3, padx = 5, pady = 5, sticky = 'W')
        so2_amt_e = ttk.Entry(param_frame, textvariable = self.so2_amt, width = 12)
        so2_amt_e.grid(row = 1, column = 4, padx = 5, pady = 5)
        self.so2_c = tk.StringVar(param_frame, value = settings['Fit so2'])
        fit_options[0] = settings['Fit so2']
        so2_c = ttk.OptionMenu(param_frame, self.so2_c, *fit_options)
        so2_c.grid(row = 1, column = 5, padx = 5, pady = 5)
        
        self.no2_amt = tk.DoubleVar(param_frame, value = settings['NO2'])
        no2_amt_l = tk.Label(param_frame, text = 'NO2:', font = NORM_FONT)
        no2_amt_l.grid(row = 2, column = 3, padx = 5, pady = 5, sticky = 'W')
        no2_amt_e = ttk.Entry(param_frame, textvariable = self.no2_amt, width = 12)
        no2_amt_e.grid(row = 2, column = 4, padx = 5, pady = 5)
        self.no2_c = tk.StringVar(param_frame, value = settings['Fit no2'])
        fit_options[0] = settings['Fit no2']
        no2_c = ttk.OptionMenu(param_frame, self.no2_c, *fit_options)
        no2_c.grid(row = 2, column = 5, padx = 5, pady = 5)
        
        self.o3_amt = tk.DoubleVar(param_frame, value = settings['O3'])
        o3_amt_l = tk.Label(param_frame, text = 'O3:', font = NORM_FONT)
        o3_amt_l.grid(row = 3, column = 3, padx = 5, pady = 5, sticky = 'W')
        o3_amt_e = ttk.Entry(param_frame, textvariable = self.o3_amt, width = 12)
        o3_amt_e.grid(row = 3, column = 4, padx = 5, pady = 5)
        self.o3_c = tk.StringVar(param_frame, value = settings['Fit o3'])
        fit_options[0] = settings['Fit o3']
        o3_c = ttk.OptionMenu(param_frame, self.o3_c, *fit_options)
        o3_c.grid(row = 3, column = 5, padx = 5, pady = 5)
       
        self.bro_amt = tk.DoubleVar(param_frame, value = settings['BrO'])
        bro_amt_l = tk.Label(param_frame, text = 'BrO:', font = NORM_FONT)
        bro_amt_l.grid(row = 4, column = 3, padx = 5, pady = 5, sticky = 'W')
        bro_amt_e = ttk.Entry(param_frame, textvariable = self.bro_amt, width = 12)
        bro_amt_e.grid(row = 4, column = 4, padx = 5, pady = 5)
        self.bro_c = tk.StringVar(param_frame, value = settings['Fit bro'])
        fit_options[0] = settings['Fit bro']
        bro_c = ttk.OptionMenu(param_frame, self.bro_c, *fit_options)
        bro_c.grid(row = 4, column = 5, padx = 5, pady = 5)
        





        
#========================================================================================         
#========================================================================================
#======================================GUI Operations====================================
#======================================================================================== 
#======================================================================================== 

    # Report exceptions in a new window
    def report_callback_exception(self, *args):
        err = traceback.format_exception(*args)
        tkMessageBox.showerror('Exception', err)

    # Close program on 'x' button
    def handler(self):
        self.quit()
        
    # Function to print text to the output box          
    def print_output(self, text, add_line = True):
        
        if add_line == True:
            # Add new line return to text
            text = text + '\n\n'
            
        else:
            text = text + '\n'
        
        # Write text with a new line
        self.text_box.insert(tk.END, text)
        
        # Scroll if needed
        self.text_box.see(tk.END)
        
        # Write to notes file
        try:
            with open(settings['notes_fname'], 'a') as a:
                a.write(text)
        
        except KeyError:
            pass
        
        # Force gui to update
        mygui.update(self)  

#========================================================================================
#========================================Begin iFit======================================
#========================================================================================
     
    # Function to begin analysis loop
    def begin(self, rt_flag):
        
        # Turn off stopping flag
        settings['stop_flag'] = False
 
        # Create flag to skip an analyis
        skip_flag = [False, 'No error']
       
        # Create common dictionary
        common = {}
        
        # Populate common with other data from the GUI
        common['wave_start']       = float(self.wave_start_e.get())
        common['wave_stop']        = float(self.wave_stop_e.get())
        common['poly_n']           = int(self.poly_n.get()) + 1
        common['ils_width']        = float(self.ils_width.get())
        common['ils_gauss_weight'] = float(self.gauss_weight.get())
        common['ldf']              = float(self.ldf.get())
        common['spectra_files']    = self.spec_fpaths
        common['dark_files']       = self.dark_fpaths
        common['spec_name']        = self.spec_name.get()
        common['dark_flag']        = int(settings['dark_flag'])
        common['flat_flag']        = int(settings['flat_flag'])
        common['solar_resid_flag'] = settings['solar_resid_flag']

        # Turn of dark flag if in real time and no darks have been taken
        if rt_flag == 'rt_analysis' and settings['rt_dark_flag'] == False:
            common['dark_flag'] = 0
            self.print_output('WARNING! No dark spectra aquired!')

#========================================================================================
#================================Build parameter dictionary==============================
#========================================================================================
        
        # Create parameter array
        params = OrderedDict()

        for i in range(common['poly_n']):
            params['p'+str(i)] = [1.0, 'Fit']
            
        # Add other parameters
        params['shift']     = [float(self.shift.get())   , self.shift_c.get()]
        params['stretch']   = [float(self.stretch.get()) , self.stretch_c.get()]
        params['ring_amt']  = [float(self.ring_amt.get()), self.ring_c.get()]
        params['so2_amt']   = [float(self.so2_amt.get()) , self.so2_c.get()]
        params['no2_amt']   = [float(self.no2_amt.get()) , self.no2_c.get()]
        params['o3_amt']    = [float(self.o3_amt.get())  , self.o3_c.get()]
        params['bro_amt']   = [float(self.bro_amt.get()) , self.bro_c.get()]
        params['ils_width'] = [self.ils_width.get()      , self.ils_width_c.get()]
        params['ldf']       = [self.ldf.get()            , self.ldf_c.get()]
        
        # Create empty last spec
        common['last_spec'] = np.array([0])

        # Add to common
        common['params'] = params
        
        # Save initial guess parameters
        initial_params = common['params'].copy()

#========================================================================================
#=================================Read in xsecs and flat=================================
#========================================================================================
        
        # Update status
        self.status.set('Building Model')
        mygui.update(self)

        # Build filepath to flat spectrum from spectrometer serial number
        if rt_flag == 'post_analysis':
            settings['flat_path'] = 'data_bases/Spectrometer/flat_' + \
                                     str(self.spec_name.get())+'.txt'
        else:
            settings['flat_path'] = 'data_bases/Spectrometer/flat_' + \
                                     str(self.c_spec.get())+'.txt'
        
        # Load fitting data files
        common = build_fwd_data(common, settings, self)

#========================================================================================
#===================================Build dark spectrum==================================
#========================================================================================

        if rt_flag == 'post_analysis':
            
            # Read format of spectra
            spec_type = self.spec_type.get()
    
            # Read in dark spectra
            if common['dark_flag'] == True:
                x, common['dark'] = average_spectra(common['dark_files'], spec_type)
                
        elif settings['rt_dark_flag'] == True and common['dark_flag'] == True:
            
            # Assign dark spectrum
            common['dark'] = self.dark_spec

#========================================================================================
#========================Read test spectrum to get wavelength grid=======================
#========================================================================================
 
        if rt_flag == 'post_analysis':
            
            # Read first spectrum to get date of data and define stray light indices
            spectrum_data = read_spectrum(common['spectra_files'][0], spec_type)
            x, y, read_date, read_time, spec_no, fit_flag = spectrum_data
            
        else:

            # Read a single spectrum to get wavelength data
            try:
                x, y, header, t = acquire_spectrum(self, settings['spec'], 1, 1)
                read_date, read_time = t.split(' ')
                
            except KeyError:
                self.print_output('No spectrometer connected')
                return
            
            except SeaBreezeError:
                self.print_output('Spectrometer disconnected')
                return 
            
        # Find indices of desired wavelength window and add to common
        grid,common['ind1'],common['ind2'] = extract_window(x, common['wave_start'],
                                                            common['wave_stop'])
        
        # Find stray light window
        stray_grid, common['stray_i1'], common['stray_i2'] = extract_window(x, 280, 290)
        
        # If no stray light pixels available, turn off the flag
        if common['stray_i1'] == common['stray_i2']:
            common['stray_flag'] = False     
        else:
            common['stray_flag'] = True
            
        # If forming solar residual create empty array and counter
        if common['solar_resid_flag'] == 'Generate':
            common['solar_resid'] = np.zeros(len(grid))
            resid_count = 0
 
#========================================================================================
#===================================Create ouput folder==================================
#========================================================================================
            
        # Create output folder to hold analysis results
        if rt_flag == 'rt_analysis':
            
            # Create new output folder
            if settings['create_out_flag'] == True:
               
                # Reset loop counter
                settings['loop'] = 0
                
                # Create filename for output file
                out_excel_fname = settings['rt_folder'] + 'iFit_out.csv'
                
                # Open excel file and write header line
                make_csv_file(out_excel_fname, common)
                
                # Create folder to hold spectra
                if not os.path.exists(settings['rt_folder'] + 'spectra/'):
                        os.makedirs(settings['rt_folder'] + 'spectra/')
                
                # Turn off flag to limit folders to one per program run
                settings['create_out_flag'] = False
                
            else:
                
                # Get final spectrum number in folder
                flist = glob.glob(settings['rt_folder'] + 'spectra/spectrum*')
                
                # Update loop number to append spectra to those in the folder
                settings['loop'] = int(flist[-1][-9:-4]) + 1
                
                # Create filename for output file
                out_excel_fname = settings['rt_folder'] + 'iFit_out.csv'
                
        else:
            
            # Reset loop counter
            settings['loop'] = 0
            
            # Create filepath to directory to hold program outputs
            post_results_folder = 'Results/iFit/' + str(read_date) + '/'
            
            # Create folder if it doesn't exist
            if not os.path.exists(post_results_folder):
                    os.makedirs(post_results_folder)
            
            # Create filename for output file
            out_excel_fname = post_results_folder + 'iFit_out.csv'
        
            try:
                
                # Open excel file and write header line
                make_csv_file(out_excel_fname, common)
                    
            except PermissionError:
                self.print_output('Please close iFit output file to continue')
                return                

#========================================================================================
#===========================Set Progress bar to correct format===========================
#========================================================================================
            
        if rt_flag == 'rt_analysis':
            self.progress['mode'] = 'indeterminate'
            self.progress['value'] = 0
        else:
            self.progress['mode'] = 'determinate'
            self.progress['value'] = 0
             
#========================================================================================
#===================================Start Analysis Loop==================================
#========================================================================================

        # Open excel file
        with open(out_excel_fname, 'a') as writer:

            # Print output message to begin
            self.print_output('Loop Started\n' +\
                              'Spectrum number ' + str(settings['loop']))  
            
            # Create empty arrays to hold the loop number and so2_amt values
            gas = {}
            spec_nos = []
            gas['SO2_amts'] = []
            gas['SO2_errs'] = []
            gas['O3_amts']  = []
            gas['O3_errs']  = []
            gas['BrO_amts'] = []
            gas['BrO_errs'] = []
        
            # Update status
            if rt_flag == 'rt_analysis':
                self.status.set('Acquiring')
            else:
                self.status.set('Analysing')
            mygui.update(self)

            # Begin analysis loop
            while True:
                            
                
#========================================================================================
#======================================Post analysis=====================================
#========================================================================================                
                
                # End loop if finished
                if settings['stop_flag'] == True:
                    break
                    
                # Read spectrum from file and fit
                if rt_flag == 'post_analysis':

                    try:
                        fname = common['spectra_files'][settings['loop']]
                        spec_data = read_spectrum(fname, spec_type)
                        x, y, read_date, read_time, spec_no, skip_flag = spec_data

                    except IndexError:
                        self.print_output('Fitting complete')
                        break 
                    
                    if skip_flag[0] == False:
                    
                        # Fit
                        results = fit_spec(common, y, grid)
                        fit_dict, err_dict, y_data, fit, gas_T, fit_flag = results
                        
                        now_fit_spec = True
                        
                        # Update progress bar
                        prog = (settings['loop']+1)/len(common['spectra_files']) * 100
                        self.progress['value'] = prog
                        
                
#========================================================================================
#=======================================RT analysis======================================
#========================================================================================  
                        
                
                # Read spectrum from spectrometer and fit
               # self.print_output(rt_flag, self.toggle_button.config('text')[-1], 
                #                  )
                elif rt_flag=='rt_analysis' and self.toggle_button.config('text')[-1]==\
                'FITTING ON' and common['last_spec'].all() != 0:

                    # Create results queue
                    result_queue = Queue()
                    
                    # Create two threads, one to read a spectrum and one to fit
                    t1 = Thread(target = acquire_spectrum, args=(self,
                                                                 settings['spec'],
                                                                 settings['int_time'],
                                                                 int(self.coadds.get()),
                                                                 True,
                                                                 True,
                                                                 result_queue))
                    
                    t2 = Thread(target = fit_spec, args = (common,
                                                           common['last_spec'],
                                                           grid,
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

                    # Get fit results
                    fit_dict, err_dict, y_data, fit, gas_T, fit_flag = thread_out['fit']
                   
                    # Get spectrum
                    x, y, header, t = thread_out['spectrum']
                    read_date, read_time = t.split(' ')
                    
                    # Build file name
                    n = str('{num:05d}'.format(num=settings['loop']))
                    fname = settings['rt_folder'] + 'spectra/spectrum_' + n + '.txt'
                    
                    # Save
                    np.savetxt(fname, np.column_stack((x, y)), header = header)
                    
                    # Update last spec variable and spec number
                    common['last_spec'] = y
                    spec_no = settings['loop']
                    
                    # Update progress bar
                    prog = settings['loop'] + 1
                    self.progress['value'] = prog
                    
                    
                    now_fit_spec = True
                
                
#========================================================================================
#====================================Just read spectrum==================================
#========================================================================================                
                
                
                # Read spectrum from spectrometer but do not fit
                else:

                    # Read spectrum
                    x, y, header, t = acquire_spectrum(self, settings['spec'], 
                                                      settings['int_time'],
                                                      int(self.coadds.get()))
                    read_date, read_time = t.split(' ')
                    
                    # Build file name
                    n = str('{num:05d}'.format(num=settings['loop']))
                    fname = settings['rt_folder'] + 'spectra/spectrum_' + n + '.txt'
                    
                    # Save
                    np.savetxt(fname, np.column_stack((x,y)), header = header)
                    
                    # Update last spec variable
                    common['last_spec'] = y
                    
                    # Update progress bar
                    prog = settings['loop']
                    self.progress['value'] = prog
                    
                    now_fit_spec = False
                
#========================================================================================
#========================Unpack fit results and add to output file=======================
#========================================================================================
                 
                if now_fit_spec == True and skip_flag[0] == False:                  
                    
                    # Calculate the residual of the fit
                    if settings['resid_type'] == 'Percentage':
                        resid=np.multiply(np.divide(np.subtract(y_data,fit),y_data), 100)
                        max_resid = resid.max() / 100
                        
                    if settings['resid_type'] == 'Absolute': 
                        resid = np.subtract(y_data, fit)
                        max_resid = np.divide(resid, fit.max()).max()
                        
                    if settings['resid_type'] == 'Spec/Fit':
                        resid = np.divide(y_data, fit)
                        max_resid = np.abs(np.subtract(resid, 1)).max()
                    
                    # Add to solar residual
                    if common['solar_resid_flag'] == 'Generate':
                        common['solar_resid'] = np.add(common['solar_resid'],
                                                       np.divide(y_data, fit))
                        resid_count += 1
                     
                    # If update params on, check fit quality. If fit fails or is bad, 
                    #   reset the first guess parameters
                    if bool(settings['update_params']) == True:
                        if fit_flag == False:
                            common['params'] = initial_params.copy()
                            self.print_output('Fitting for spectrum ' + str(spec_no) + \
                                              ' failed, resetting parameters')
                            fit_msg = 'Failed'
                           
                        elif max_resid > float(settings['good_fit_bound'])/100:
                            common['params'] = initial_params.copy()
                            self.print_output('Fitting for spectrum ' + str(spec_no) + \
                                              ' bad, resetting parameters')
                            fit_msg = 'Bad'
                       
                        else:
                            
                            # Update first guesses with last fitted params
                            for key, val in fit_dict.items():
                                common['params'][key][0] = val
                            fit_msg = 'Good'
    
                    # Write results to excel file, starting with spectrum info
                    writer.write(str(fname)          + ',' + \
                                 str(spec_no)        + ',' + \
                                 str(str(read_date)) + ',' + \
                                 str(str(read_time)))
                    
                    # Print so2 amount and error in ppm.m for ease
                    #writer.write(',' + str(fit_dict['so2_amt']/2.463e15) + ',' + \
                    #             str(err_dict['so2_amt']/2.463e15))            
        
                    # Print fit results and error for each parameter          
                    for key, val in common['params'].items():

                        if val[1] == 'Fit':
                            writer.write(','+str(fit_dict[key])+','+str(err_dict[key]))
                        if val[1] == 'Fix':
                            writer.write(','+str(common['params'][key][0])+',NaN')
                        if val[1] == 'N/A':
                            writer.write(',NaN,NaN')
                        
                    # Write fit quality and start new line
                    writer.write(',' + fit_msg + '\n')
                
                    # Add values to array for plotting
                    spec_nos.append(spec_no)
                    if common['params']['so2_amt'][1] == 'Fit':
                        gas['SO2_amts'].append(fit_dict['so2_amt']/2.463e15)
                        gas['SO2_errs'].append(err_dict['so2_amt']/2.463e15)
                    if common['params']['o3_amt'][1] == 'Fit':
                        gas['O3_amts'].append(fit_dict['o3_amt']/2.463e15)
                        gas['O3_errs'].append(err_dict['o3_amt']/2.463e15)
                    if common['params']['bro_amt'][1] == 'Fit':
                        gas['BrO_amts'].append(fit_dict['bro_amt']/2.463e15)
                        gas['BrO_errs'].append(err_dict['bro_amt']/2.463e15)
                    
                    # Update quick analysis with values
                    last_amt="{0:0.2f}".format(gas[settings['analysis_gas']+'_amts'][-1])
                    last_err="{0:0.2f}".format(gas[settings['analysis_gas']+'_errs'][-1])
                    self.last_so2_amt.set(last_amt + ' ppm.m')
                    self.last_so2_err.set(last_err + ' ppm.m')
                    
                    # Cut if too long to avoid slowing program
                    if bool(settings['scroll_flag']) == True:
                        lim = int(settings['scroll_spec_no'])
                        if len(spec_nos) > lim:
                            spec_nos = spec_nos[1:]
                            for m in gas:
                                gas[m] = gas[m][1:]
                                
#========================================================================================
#=======================================Update plot======================================
#========================================================================================
        
                # Replot data
                if int(settings['Show Graphs']) == 1 and skip_flag[0] == False:            
                    
                    if now_fit_spec == True:
                        
                        # Get selected transmittance data
                        gas_tran = -np.log(gas_T[settings['analysis_gas'] + '_tran'])
                        gas_spec = -np.log(gas_T[settings['analysis_gas'] + '_spec'])
                        gas_amts = gas[settings['analysis_gas'] + '_amts']
                    
                        # Build axes and lines arrays
                        lines = [self.line0, self.line1, self.line2, self.line3, 
                                 self.line4, self.line5, self.line6]
                        axes =  [self.ax0,   self.ax0,   self.ax1,   self.ax2,  
                                 self.ax3,   self.ax3,   self.ax4  ]
                        
                        # Calculate graph limits for spectrum fit
                        y_lo = min(y_data) - abs(max(y_data) - min(y_data)) * 0.1
                        y_hi = max(y_data) + abs(max(y_data) - min(y_data)) * 0.1
                        f_lo = min(fit) - abs(max(fit) - min(fit)) * 0.1
                        f_hi = max(fit) + abs(max(fit) - min(fit)) * 0.1
                        y_lo, y_hi = min([f_lo, y_lo]), max([f_hi, y_hi])
                        
                        # Limits for so2 trans spectrum
                        t_lo = min(gas_tran) - abs((0.1*min(gas_tran)))
                        t_hi = max(gas_tran) + (0.1*max(gas_tran)) 
                        s_lo = min(gas_spec) - abs((0.1*min(gas_spec)))
                        s_hi = max(gas_spec) + (0.1*max(gas_spec))
                        t_lo, t_hi = min([t_lo, s_lo]), max([t_hi, s_hi])
                        
                        
                        # Build data array to pass to graphing function
                        #                 x data    y data    x limits     y limits
                        data = np.array(([grid,     y_data,   'auto'     , [y_lo,y_hi]],
                                         [grid,     fit,      'auto'     , [y_lo,y_hi]],
                                         [x   ,     y,        'auto'     , 'auto'     ],
                                         [grid,     resid,    'auto'     , 'auto'     ],
                                         [grid,     gas_tran, 'auto'     , [t_lo,t_hi]],
                                         [grid,     gas_spec, 'auto'     , [t_lo,t_hi]],
                                         [spec_nos, gas_amts, 'auto'     , 'auto'     ]))
                    
                    else:
                        
                        # Build axes and lines arrays
                        lines = [self.line2]
                        axes =  [self.ax1  ]
                        
                        # Calculate limits
                        y_lo  = min(y) - abs((0.1*max(y)))
                        y_hi  = max(y) + abs((0.1*max(y)))
                        x_lo, x_hi = x.min() - 1, x.max() + 1
                        
                        # Build data array to pass to graphing function
                        #                 x data    y data    x limits     y limits
                        data = np.array(([x,        y,        [x_lo,x_hi], [y_lo,y_hi]]))
                    
                    # Update graph
                    update_graph(lines, axes, self.canvas, data)
                    
                    # Make it look nice
                    plt.tight_layout()
                
                if skip_flag[0] == True:
                    
                    self.print_output('Error reading spectrum:\n' + str(skip_flag[1]))
                               
                # Add to the count cycle
                settings['loop'] += 1
                
                # Force gui to update
                mygui.update(self)
                
            # Update solar residual
            if settings['solar_resid_flag'] == 'Generate':
                
                # Find average residual
                common['solar_resid'] = np.divide(common['solar_resid'], resid_count)
                
                # Save
                np.savetxt('data_bases/gas data/solar_resid.txt', 
                           np.column_stack((grid,common['solar_resid'])))
                
                self.print_output('Solar residual spectrum updated')
                
            # Update status
            self.status.set('Standby')
            
#========================================================================================
#=======================================Save Settings====================================
#========================================================================================
    
    # Function to save setting to the ifit_settings.txt file        
    def save(self):
        
        # Create or overright settings file
        with open('data_bases/ifit_settings.txt', 'w') as w:
            
            # Save each setting from the gui into settings
            settings['Wave Start']   = str(self.wave_start_e.get())   
            settings['Wave Stop']    = str(self.wave_stop_e.get())      
            settings['Spectrometer'] = str(self.spec_name.get())       
            settings['Spectra Type'] = str(self.spec_type.get())       
            settings['int_time']     = str(self.int_time.get())         
            settings['coadds']       = str(self.coadds.get())
            settings['no_darks']     = str(self.no_darks.get())
            settings['ILS Width']    = str(self.ils_width.get())
            settings['Gauss Weight'] = str(self.gauss_weight.get())
            settings['Fit ILS']      = str(self.ils_width_c.get())
            settings['LDF']          = str(self.ldf.get())
            settings['Fit LDF']      = str(self.ldf_c.get())
            settings['poly_n']       = str(self.poly_n.get())
            settings['shift']        = str(self.shift.get())
            settings['stretch']      = str(self.stretch.get())
            settings['ring']         = str(self.ring_amt.get())
            settings['SO2']          = str(self.so2_amt.get())
            settings['NO2']          = str(self.no2_amt.get())
            settings['O3']           = str(self.o3_amt.get())
            settings['BrO']          = str(self.bro_amt.get())
            settings['Fit shift']    = str(self.shift_c.get())
            settings['Fit stretch']  = str(self.stretch_c.get())
            settings['Fit ring']     = str(self.ring_c.get())
            settings['Fit so2']      = str(self.so2_c.get())
            settings['Fit no2']      = str(self.no2_c.get())
            settings['Fit o3']       = str(self.o3_c.get())
            settings['Fit bro']      = str(self.bro_c.get())
            
            # Add all of the settings dictionary
            for s in settings:
                w.write(s + ';' + str(settings[s]) + '\n')
            
            try:
                w.write('Spectra Filepaths;' + str(self.spec_fpaths) + '\n')
            except AttributeError:
                w.write('Spectra Filepaths; \n') 
            
            try:
                w.write('Dark Filepaths;' + str(self.dark_fpaths))
            except AttributeError:
                w.write('Dark Filepaths; ')
                
        self.print_output('Settings saved')

    
# Run the App!
if __name__ == '__main__':    
    mygui().mainloop()
