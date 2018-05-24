# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 09:24:05 2018

@author: mqbpwbe2
"""

# Import required libraries
import matplotlib
matplotlib.use('TkAgg')
import os
import glob
import numpy as np
from tkinter import ttk
import tkinter as tk
from tkinter import filedialog as fd
import seabreeze.spectrometers as sb
from queue import Queue
from threading import Thread
import datetime
from seabreeze.cseabreeze.wrapper import SeaBreezeError
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import OrderedDict
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

from ifit_lib.build_fwd_data import build_fwd_data
from ifit_lib.read_spectrum import read_spectrum, average_spectra
from ifit_lib.fit import fit_spec, ifit_fwd
from ifit_lib.aquire_spectrum import aquire_spectrum
from ifit_lib.find_nearest import extract_window
from ifit_lib.file_control import make_directory
from ifit_lib.update_graph import update_graph

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
        
        # Close program on closure of window
        #mygui.protocol("WM_DELETE_WINDOW", mygui.destroy(self))
        
        # Button Style
        ttk.Style().configure('TButton', width = 20, height = 20, relief="flat") 
        
        # Add a title and icon
        tk.Tk.wm_title(self, 'iFit-2-2')
        tk.Tk.iconbitmap(self, default = 'data_bases/icon.ico')
        
        # Build a menubar to hold options for the user
        menubar = tk.Menu(self)
        filemenu = tk.Menu(menubar, tearoff = 0)
        filemenu.add_command(label = 'Model Settings', command = model_settings)
        filemenu.add_command(label = 'Graph Settings', command = graph_settings)
        filemenu.add_command(label = 'Data Bases', command = data_bases)
        filemenu.add_separator()
        menubar.add_cascade(label = 'Options', menu = filemenu)
        tk.Tk.config(self, menu = menubar)
        
        # Create notebook to hold different frames
        nb = ttk.Notebook(self)
        page1 = ttk.Frame(nb)
        page2 = ttk.Frame(nb)
        page3 = ttk.Frame(nb)
        
        # Create two frames, one for post analysis and one for real time acquisition
        nb.add(page2, text = 'Real Time Acquisition')
        nb.add(page1, text = 'Post Analysis')
        nb.add(page3, text = 'Model Parameters')
        
        nb.grid(column=0, padx=5, pady=5)
        
        # Create frame to hold graphs
        graph_frame = ttk.Frame(self, relief = 'groove')
        graph_frame.grid(row=0, column=1, padx=10, pady=10, rowspan=10, sticky="NW")
        graph_frame.columnconfigure(index=0, weight=1)
        graph_frame.rowconfigure(index = 0, weight = 1)
        
        # Create frame to hold text output
        text_frame = ttk.Frame(self, relief = 'groove')
        text_frame.grid(row=1, column=0, padx=10, pady=10, rowspan=10, sticky="NW")
        
        mygui.columnconfigure(index=1, weight=1, self = self)
        mygui.rowconfigure(index = 5, weight = 1, self = self)
        
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
            settings['Wave Start']        = 306
            settings['Wave Stop']         = 318
            settings['Spectrometer']      = '-select-'
            settings['Spectra Type']      = '-select-'
            settings['int_time']          = 100
            settings['coadds']            = 10
            settings['no_darks']          = 10
            settings['ILS Width']         = 0.52
            settings['Gauss Weight']      = 1.0
            settings['Fit ILS']           = 0
            settings['LDF']               = 0.0
            settings['Fit LDF']           = 0
            settings['dark_flag']         = 1
            settings['flat_flag']         = 1
            settings['update_params']     = 1
            settings['Show Graphs']       = 1
            settings['Show Error Bars']   = 0
            settings['scroll_flag']       = 1
            settings['scroll_spec_no']    = 200
            settings['poly_n']            = 4
            settings['shift']             = -0.2
            settings['stretch']           = 0.05
            settings['ring']              = 1.0
            settings['SO2']               = 1e+16
            settings['NO2']               = 1e+17
            settings['O3']                = 1e+19
            settings['BrO']               = 1e+15
            settings['model_resolution']  = 0.52
            settings['sol_path']          = 'data_bases/sao2010_full.txt'
            settings['ring_path']         = 'data_bases/ring.dat'
            settings['so2_path']          = 'data_bases/SO2_293K.dat'
            settings['no2_path']          = 'data_bases/No2_223l.dat'
            settings['o3_path']           = 'data_bases/o3_223l.dat'
            settings['bro_path']          = 'data_bases/BrO_Cross_298K.txt'
            settings['Spectra Filepaths'] = ''
            settings['Dark Filepaths']    = ''
 
        # Create loop counter to keep track of the analysis
        settings['loop'] = 0
        
        # Create flag to ensure only one output file is created per program launch
        settings['create_out_flag'] = True
        
        # Create flag to see if darks have been measured yet
        settings['rt_dark_flag'] = False
        
        
#========================================================================================
#==================================Create output folder==================================
#======================================================================================== 

        # Create filepath to directory to hold program outputs
        rt_results_folder = 'Results/iFit/' + str(datetime.date.today()) + '/ifit_out/'
       
        # Create folder
        settings['rt_folder'] = make_directory(rt_results_folder, overwrite = True)
        
        # Create notes file
        settings['notes_fname'] = settings['rt_folder'] + 'notes.txt'
        with open(settings['notes_fname'], 'w') as w:
            w.write('Notes file for iFit\n\n')
        
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
        self.ax1.set_xlabel('Wavelength (nm)', fontsize=10)
        
        self.ax2.set_ylabel('Fit residual (%)', fontsize=10)
        self.ax2.set_xlabel('Wavelength (nm)', fontsize=10)
        
        self.ax3.set_ylabel('SO2 amt (ppm.m)', fontsize=10)
        self.ax3.set_xlabel('Spectrum number', fontsize=10)
        
        self.ax4.set_ylabel('SO2 amt (ppm.m)', fontsize=10)
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
        self.line4, = self.ax3.plot(0, 0, 'b', label = 'y / F_no_so2')
        self.line5, = self.ax3.plot(0, 0, 'r', label = 'so2 trans')
        self.ax3.legend(loc = 0)
        
        # SO2 Time series and error bars
        self.line6, = self.ax4.plot(0, 0, 'g')
        
        # Make it look nice
        plt.tight_layout()
        
        # Create the canvas to hold the graph in the GUI
        self.canvas = FigureCanvasTkAgg(self.fig, graph_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady = 10)
        
        # Add matplotlib toolbar above the plot canvas
        toolbar_frame = tk.Frame(graph_frame, bg = 'black')  
        toolbar_frame.grid(row=1,column=0, sticky = 'W', padx = 5, pady = 5)                             
        toolbar = NavigationToolbar2TkAgg(self.canvas, toolbar_frame)
        toolbar.update()
        
#========================================================================================
#====================================Create text output==================================
#========================================================================================        
                 
        # Create a scroll bar
        scrollbar = ttk.Scrollbar(text_frame)
        scrollbar.grid(row=1, column=4, padx = 5, pady = 5, sticky='NSE')
        
        # Build text box
        self.text_box = tk.Text(text_frame, width = 60, height = 10, 
                                yscrollcommand = scrollbar.set)
        self.text_box.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W',
                           columnspan = 4)
        self.text_box.insert('1.0', 'Welcome to iFit! Written by Ben Esse\n\n')
        
        # Connect the scrollbar to the textbox
        scrollbar.config(command = self.text_box.yview)
        
        # Create button to save settings
        save_b = ttk.Button(text_frame, text = 'Save Settings', command = self.save)
        save_b.grid(row = 0, column = 0, columnspan = 2, padx = 5, pady = 5)
        
        # Create a button to exit
        exit_b = ttk.Button(text_frame, text = 'Exit', command = self.quit)
        exit_b.grid(row = 0, column = 2, columnspan = 2, padx = 5, pady = 5)


          



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
        
        # Frame for quick analysis
        quick_frame = tk.LabelFrame(page1, text = 'Quick Analysis', font = LARG_FONT)
        quick_frame.grid(row=2, column=0, padx=10, pady=10, sticky="NW")
      
#========================================================================================
#==================================Create control inputs=================================
#======================================================================================== 
               
        # Find available flat spectra and form into list
        options = [settings['Spectrometer']]       
        for i, name in enumerate(glob.glob('data_bases/flat_*')):
            options.append(name[16:-4])
        
        # Create entry to select spectrometer
        self.spec_name = tk.StringVar(setup_frame, value = options[0])
        spectro_l = tk.Label(setup_frame, text = 'Spectrometer:', font = NORM_FONT)
        spectro_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        spectro_c = ttk.OptionMenu(setup_frame, self.spec_name, *options)
        spectro_c.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        # Create entry to select spectra type
        spec_options = [settings['Spectra Type'],
                        'IFRiT',
                        'Master.Scope',
                        'Jai Spec',
                        'Spectrasuite']
        
        self.spec_type = tk.StringVar(setup_frame, value = spec_options[0])
        self.spec_l = tk.Label(setup_frame, text = 'Spectra Type:', font = NORM_FONT)
        self.spec_l.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.spec_c = ttk.OptionMenu(setup_frame, self.spec_type, *spec_options)
        self.spec_c.grid(row = 0, column = 3, padx = 5, pady = 5, sticky = 'W')
        
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
        specfp_b = ttk.Button(setup_frame, text="Select Spectra", command = self.spec_fp)
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
        darkfp_b = ttk.Button(setup_frame, text = "Select Darks", command = self.dark_fp)
        darkfp_b.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
               
#========================================================================================
#===================================Create start button==================================
#========================================================================================         
        
        # Create button to start
        start_b = ttk.Button(button_frame, text = 'Begin!', 
                             command = lambda: self.begin('post_analysis'))
        start_b.grid(row = 0, column = 0, padx = 40, pady = 5, columnspan = 2)
        
        # Create button to stop
        stop_b = ttk.Button(button_frame, text = 'Stop', command = self.stop)
        stop_b.grid(row = 0, column = 2, padx = 40, pady = 5, columnspan = 2)
        
#========================================================================================
#==============================Create quick analysis outputs=============================
#========================================================================================
        
        # Create ouput for last so2 amount
        self.last_so2_amt = tk.DoubleVar(self, value = 0)
        last_so2_amt_l = tk.Label(quick_frame, text = 'SO2 amount:', 
                                  font = NORM_FONT)
        last_so2_amt_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        last_so2_amt_e = ttk.Entry(quick_frame, textvariable = self.last_so2_amt)
        last_so2_amt_e.grid(row = 0, column = 1, padx = 5, pady = 5)
        
        # Create ouput for last so2 error
        self.last_so2_err = tk.DoubleVar(self, value = 0)
        last_so2_err_l = tk.Label(quick_frame, text = 'SO2 error:', 
                                  font = NORM_FONT)
        last_so2_err_l.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = 'W')
        last_so2_err_e = ttk.Entry(quick_frame, textvariable = self.last_so2_err)
        last_so2_err_e.grid(row = 0, column = 3, padx = 5, pady = 5)
        
     

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
        
        # Frame for quick analysis
        quick_frame2 = tk.LabelFrame(page2, text = 'Quick Analysis', font = LARG_FONT)
        quick_frame2.grid(row=2, column=0, padx=10, pady=10, sticky="NW")

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
                                    command = self.connect_spec)
        connect_spec_b.grid(row = 0, column = 2, padx = 5, pady = 5)
        
        # Create button to update integration time
        update_int_time_b = ttk.Button(setup_frame2, text = 'Update',
                                       command = self.update_int_time)
        update_int_time_b.grid(row = 1, column = 2, padx = 5, pady = 5)
        
        # Create button to read a single spectrum
        test_spec_b = ttk.Button(setup_frame2, text = 'Test Spectrum',
                                 command = self.test_spec)
        test_spec_b.grid(row = 2, column = 2, padx = 5, pady = 5)
        
        # Create button to read darks
        read_darks_b = ttk.Button(setup_frame2, text = 'Read Darks',
                                  command = self.read_darks)
        read_darks_b.grid(row = 3, column = 2, padx = 5, pady = 5)
        
#========================================================================================
#==============================Create start and exit buttons=============================
#========================================================================================         
        
        # Create button to start
        start_aq_b = ttk.Button(button_frame2, text = 'Begin!', 
                                command = lambda: self.begin('rt_analysis'))
        start_aq_b.grid(row = 0, column = 0, padx = 40, pady = 5)
        
        # Create button to stop
        stop_aq_b = ttk.Button(button_frame2, text = 'Stop', command = self.stop)
        stop_aq_b.grid(row = 0, column = 1, padx = 40, pady = 5)
        
        # Create switch to toggle fitting on or off
        self.toggle_button = tk.Button(button_frame2, text = 'FITTING OFF', 
                                       command = self.fit_toggle, width = 14, height = 3,
                                       bg = 'red', font = LARG_FONT)
        self.toggle_button.grid(row=1, column=0, padx=5, pady=5, columnspan=2)


#========================================================================================
#==============================Create quick analysis outputs=============================
#========================================================================================
        
        # Create ouput for last so2 amount
        last_so2_amt_l = tk.Label(quick_frame2, text = 'SO2 amount:', 
                                  font = NORM_FONT)
        last_so2_amt_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        last_so2_amt_e = ttk.Entry(quick_frame2, textvariable = self.last_so2_amt)
        last_so2_amt_e.grid(row = 0, column = 1, padx = 5, pady = 5)
        
        # Create ouput for last so2 error
        last_so2_err_l = tk.Label(quick_frame2, text = 'SO2 error:', 
                                  font = NORM_FONT)
        last_so2_err_l.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = 'W')
        last_so2_err_e = ttk.Entry(quick_frame2, textvariable = self.last_so2_err)
        last_so2_err_e.grid(row = 0, column = 3, padx = 5, pady = 5)
        
        
        
        



#========================================================================================         
#========================================================================================
#=================================Model Parameters Frame=================================
#======================================================================================== 
#======================================================================================== 






#========================================================================================
#====================================Create GUI frames===================================
#========================================================================================

        # Frame for model settings
        model_frame = tk.LabelFrame(page3, text = 'Model Setup', font = LARG_FONT)
        model_frame.grid(row=0, column=0, padx=10, pady=10, sticky="NW")

        # Frame for parameters
        param_frame = tk.LabelFrame(page3, text = 'Fit Paramters', font = LARG_FONT)
        param_frame.grid(row=1, column=0, padx=10, pady=10, sticky="NW")

#========================================================================================
#=======================================Model Setup======================================
#========================================================================================

        # Create entry for start and stop wavelengths
        self.wave_start = tk.IntVar(model_frame, value = settings['Wave Start'])
        self.wave_start_l = tk.Label(model_frame, text = 'Wave Start:', font = NORM_FONT)
        self.wave_start_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.wave_start_e = ttk.Entry(model_frame, textvariable = self.wave_start)
        self.wave_start_e.grid(row = 0, column = 1, padx = 5, pady = 5)
        
        self.wave_stop = tk.IntVar(model_frame, value = settings['Wave Stop'])
        self.wave_stop_l = tk.Label(model_frame, text = 'Wave Stop:', font = NORM_FONT)
        self.wave_stop_l.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.wave_stop_e = ttk.Entry(model_frame, textvariable = self.wave_stop)
        self.wave_stop_e.grid(row = 0, column = 3, padx = 5, pady = 5)
        
        # Instrument lineshape width
        self.ils_width = tk.IntVar(model_frame, value = settings['ILS Width'])
        self.ils_width_b = tk.IntVar(model_frame, value = settings['Fit ILS'])
        self.ils_width_l = tk.Label(model_frame, text = 'ILS Width:', font = NORM_FONT)
        self.ils_width_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.ils_width_e = ttk.Entry(model_frame, textvariable = self.ils_width)
        self.ils_width_e.grid(row = 1, column = 1, padx = 5, pady = 5)
        self.ils_width_l2 = tk.Label(model_frame, text = 'Fit ILS?', font = NORM_FONT)
        self.ils_width_l2.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.ils_width_c = ttk.Checkbutton(model_frame, variable = self.ils_width_b)
        self.ils_width_c.grid(row = 1, column = 3, padx = 5, pady = 5)
        
        # ILS Gaussian weighting
        self.gauss_weight = tk.DoubleVar(model_frame, value = settings['Gauss Weight'])
        gauss_weight_l = tk.Label(model_frame, text = 'ILS Gauss Weight:', 
                                  font = NORM_FONT)
        gauss_weight_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        gauss_weight_e = ttk.Entry(model_frame, textvariable = self.gauss_weight)
        gauss_weight_e.grid(row = 2, column = 1, padx = 5, pady = 5)
        
        # Light dilution factor
        self.ldf = tk.IntVar(model_frame, value = settings['LDF'])
        self.ldf_b = tk.IntVar(model_frame, value = settings['Fit LDF'])
        self.ldf_l = tk.Label(model_frame, text = 'LDF:', font = NORM_FONT)
        self.ldf_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.ldf_e = ttk.Entry(model_frame, textvariable = self.ldf)
        self.ldf_e.grid(row = 3, column = 1, padx = 5, pady = 5)
        self.ldf_l2 = tk.Label(model_frame, text = 'Fit LDF?', font = NORM_FONT)
        self.ldf_l2.grid(row = 3, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.ldf_c = ttk.Checkbutton(model_frame, variable = self.ldf_b)
        self.ldf_c.grid(row = 3, column = 3, padx = 5, pady = 5)

#========================================================================================
#================================Set initial fit parameters==============================
#========================================================================================
        
        # Polynomial coefficents
        self.poly_n = tk.DoubleVar(self, value = settings['poly_n'])
        poly_n_l = tk.Label(param_frame, text = 'Polynomial Order:', font = NORM_FONT)
        poly_n_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        poly_n_e = ttk.Entry(param_frame, textvariable = self.poly_n)
        poly_n_e.grid(row = 1, column = 1, padx = 5, pady = 5)
        
        # Spectrometer wavelength shift parameters
        self.shift = tk.DoubleVar(param_frame, value = settings['shift'])
        shift_l = tk.Label(param_frame, text = 'shift:', font = NORM_FONT)
        shift_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        shift_e = ttk.Entry(param_frame, textvariable = self.shift)
        shift_e.grid(row = 2, column = 1, padx = 5, pady = 5)
        
        self.stretch = tk.DoubleVar(param_frame, value = settings['stretch'])
        stretch_l = tk.Label(param_frame, text = 'stretch:', font = NORM_FONT)
        stretch_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
        stretch_e = ttk.Entry(param_frame, textvariable = self.stretch)
        stretch_e.grid(row = 3, column = 1, padx = 5, pady = 5)
        
        # Ring effect
        self.ring_amt = tk.DoubleVar(param_frame, value = settings['ring'])
        ring_amt_l = tk.Label(param_frame, text = 'Ring:', font = NORM_FONT)
        ring_amt_l.grid(row = 4, column = 0, padx = 5, pady = 5, sticky = 'W')
        ring_amt_e = ttk.Entry(param_frame, textvariable = self.ring_amt)
        ring_amt_e.grid(row = 4, column = 1, padx = 5, pady = 5)
        
        # Gas amounts
        self.so2_amt = tk.DoubleVar(param_frame, value = settings['SO2'])
        so2_amt_l = tk.Label(param_frame, text = 'SO2:', font = NORM_FONT)
        so2_amt_l.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
        so2_amt_e = ttk.Entry(param_frame, textvariable = self.so2_amt)
        so2_amt_e.grid(row = 1, column = 3, padx = 5, pady = 5)
        
        self.no2_amt = tk.DoubleVar(param_frame, value = settings['NO2'])
        no2_amt_l = tk.Label(param_frame, text = 'NO2:', font = NORM_FONT)
        no2_amt_l.grid(row = 2, column = 2, padx = 5, pady = 5, sticky = 'W')
        no2_amt_e = ttk.Entry(param_frame, textvariable = self.no2_amt)
        no2_amt_e.grid(row = 2, column = 3, padx = 5, pady = 5)
        
        self.o3_amt = tk.DoubleVar(param_frame, value = settings['O3'])
        o3_amt_l = tk.Label(param_frame, text = 'O3:', font = NORM_FONT)
        o3_amt_l.grid(row = 3, column = 2, padx = 5, pady = 5, sticky = 'W')
        o3_amt_e = ttk.Entry(param_frame, textvariable = self.o3_amt)
        o3_amt_e.grid(row = 3, column = 3, padx = 5, pady = 5)
       
        self.bro_amt = tk.DoubleVar(param_frame, value = settings['BrO'])
        bro_amt_l = tk.Label(param_frame, text = 'BrO:', font = NORM_FONT)
        bro_amt_l.grid(row = 4, column = 2, padx = 5, pady = 5, sticky = 'W')
        bro_amt_e = ttk.Entry(param_frame, textvariable = self.bro_amt)
        bro_amt_e.grid(row = 4, column = 3, padx = 5, pady = 5)
        





        
#========================================================================================         
#========================================================================================
#=====================================Button Functions===================================
#======================================================================================== 
#======================================================================================== 

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
#=======================================Save Settings====================================
#========================================================================================
    
    # Function to save setting to the ifit_settings.txt file        
    def save(self):
        
        # Create or overright settings file
        with open('data_bases/ifit_settings.txt', 'w') as w:
            
            # Save each setting into the file
            w.write('Wave Start;'       + str(self.wave_start_e.get())      + '\n')
            w.write('Wave Stop;'        + str(self.wave_stop_e.get())       + '\n')
            w.write('Spectrometer;'     + str(self.spec_name.get())         + '\n')
            w.write('Spectra Type;'     + str(self.spec_type.get())         + '\n')
            w.write('int_time;'         + str(self.int_time.get())          + '\n')
            w.write('coadds;'           + str(self.coadds.get())            + '\n')
            w.write('no_darks;'         + str(self.no_darks.get())          + '\n')
            w.write('ILS Width;'        + str(self.ils_width_e.get())       + '\n')
            w.write('Gauss Weight;'     + str(self.gauss_weight.get())      + '\n')
            w.write('Fit ILS;'          + str(self.ils_width_b.get())       + '\n')
            w.write('LDF;'              + str(self.ldf_e.get())             + '\n')
            w.write('Fit LDF;'          + str(self.ldf_b.get())             + '\n')
            w.write('dark_flag;'        + str(settings['dark_flag'])        + '\n')
            w.write('flat_flag;'        + str(settings['flat_flag'])        + '\n')
            w.write('update_params;'    + str(settings['update_params'])    + '\n')
            w.write('Show Graphs;'      + str(settings['Show Graphs'])      + '\n')
            w.write('Show Error Bars;'  + str(settings['Show Error Bars'])  + '\n')
            w.write('scroll_flag;'      + str(settings['scroll_flag'])      + '\n')
            w.write('scroll_spec_no;'   + str(settings['scroll_spec_no'])   + '\n')
            w.write('poly_n;'           + str(self.poly_n.get())            + '\n')
            w.write('shift;'            + str(self.shift.get())             + '\n')
            w.write('stretch;'          + str(self.stretch.get())           + '\n')
            w.write('ring;'             + str(self.ring_amt.get())          + '\n')
            w.write('SO2;'              + str(self.so2_amt.get())           + '\n')
            w.write('NO2;'              + str(self.no2_amt.get())           + '\n')
            w.write('O3;'               + str(self.o3_amt.get())            + '\n')
            w.write('BrO;'              + str(self.bro_amt.get())           + '\n')
            w.write('model_resolution;' + str(settings['model_resolution']) + '\n')
            w.write('sol_path;'         + str(settings['sol_path'])         + '\n')
            w.write('ring_path;'        + str(settings['ring_path'])        + '\n')
            w.write('so2_path;'         + str(settings['so2_path'])         + '\n')
            w.write('no2_path;'         + str(settings['no2_path'])         + '\n')
            w.write('o3_path;'          + str(settings['o3_path'])          + '\n')
            w.write('bro_path;'         + str(settings['bro_path'])         + '\n')
            
            try:
                w.write('Spectra Filepaths;' + str(self.spec_fpaths) + '\n')
            except AttributeError:
                w.write('Spectra Filepaths; \n') 
            
            try:
                w.write('Dark Filepaths;' + str(self.dark_fpaths))
            except AttributeError:
                w.write('Dark Filepaths; ')
                
        self.print_output('Settings saved')

#========================================================================================
#======================================Print Outputs=====================================
#========================================================================================
     
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
#=================================Connect to Spectrometer================================
#========================================================================================
        
    # Function to connect to the attached spectrometer
    def connect_spec(self):
        
        # Find connected spectrometers
        devices = sb.list_devices()
        
        # If no devices are connected then set string to show. Else assign first to spec
        if len(devices) == 0:
            settings['spec'] = 0
            settings['Spectrometer'] = 'No devices connected'
            devices = ['No devices connected']
        else:
            try:
                # Connect to spectrometer
                settings['spec'] = sb.Spectrometer(devices[0])
                
                # Set intial integration time
                settings['spec'].integration_time_micros(float(self.int_time.get())*1000)
                
                # Record serial number in settings
                settings['Spectrometer'] = str(settings['spec'].serial_number)
                
                self.print_output('Spectrometer '+settings['Spectrometer']+' Connected')
                
            except SeaBreezeError:
                self.print_output('Spectrometer already open')
            
        # Update text to show spectrometer name
        self.c_spec.set(settings['Spectrometer'])
        
#========================================================================================
#=================================Connect to Spectrometer================================
#========================================================================================

    def update_int_time(self):
        
        try:
            # Update integration time on spectrometer
            settings['spec'].integration_time_micros(float(self.int_time.get())*1000)
            
            self.print_output('Integration time updated to '+str(settings['int_time']) +\
                              '\nSpectrum number ' + str(settings['loop']))
            
        except KeyError:
            self.print_output('No spectrometer conected')

#========================================================================================
#===============================Read a single test spectrum==============================
#========================================================================================

    def test_spec(self):
        x, y, header, read_time = aquire_spectrum(self,
                                                  settings['spec'], 
                                                  settings['int_time'], 
                                                  int(self.coadds.get()))
        
        # Display the test spectrum
        self.line2.set_data(x, y)
        self.ax1.set_xlim(x.min(), x.max())
        self.ax1.set_ylim(y.min(), y.max())
        self.canvas.draw()

#========================================================================================
#========================================Read Darks======================================
#========================================================================================

    def read_darks(self):
        
        # Update notes
        self.print_output('Begin reading darks')
        
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
                spc = aquire_spectrum(self, settings['spec'], settings['int_time'],
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
        
        # Divide by number of darks to get average
        settings['dark_spec'] = np.divide(dark, dark_n)
        
        # Display the dark spectrum
        lines = [self.line2]
        axes =  [self.ax1  ]
        
        # Calculate limits
        y_lo  = min(dark) - abs((0.1*max(dark)))
        y_hi  = max(dark) + abs((0.1*max(dark)))
        x_lo, x_hi = x.min() - 1, x.max() + 1
        
        # Build data array to pass to graphing function
        #                 x data    y data    x limits     y limits
        data = np.array(([x,        dark,     [x_lo,x_hi], [y_lo,y_hi]]))
        
        # Update graph
        update_graph(lines, axes, self, data)
        
        # Add dark to common
        settings['dark_spec']  = dark
        
        # Update notes file
        self.print_output('Dark updated\n' + \
                          'Spectrum no: ' + str(settings['loop']) + '\n' + \
                          'No. darks: ' + str(dark_n) + '\n' + \
                          'Integration time (ms): ' + str(settings['int_time']) + '\n' +\
                          'Coadds: ' + str(self.coadds.get()))
         
        # Set dark_flag to True
        settings['rt_dark_flag'] = True
        
#========================================================================================
#========================================Read Darks======================================
#========================================================================================
        
    def stop(self):
        settings['stop_flag'] = True
        self.print_output('Loop Stopped\nSpectrum number ' + str(settings['loop']-1))
        
#========================================================================================
#=======================================Toggle fitting===================================
#========================================================================================
            
    # Function to toggle fitting on and off
    def fit_toggle(self):
        
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
#========================================Begin iFit======================================
#========================================================================================
     
    # Function to begin analysis loop
    def begin(self, rt_flag):
        
        # Turn off stopping falg
        settings['stop_flag'] = False
        
        # Create common dictionary
        common = {}
        
        # Populate common with other data from the GUI
        common['wave_start']       = float(self.wave_start_e.get())
        common['wave_stop']        = float(self.wave_stop_e.get())
        common['poly_n']           = int(self.poly_n.get())
        common['ils_width']        = float(self.ils_width_e.get())
        common['ils_gauss_weight'] = float(self.gauss_weight.get())
        common['ldf']              = float(self.ldf_e.get())
        common['ils_flag']         = self.ils_width_b.get()
        common['ldf_flag']         = self.ldf_b.get()
        common['spectra_files']    = self.spec_fpaths
        common['dark_files']       = self.dark_fpaths
        common['spec_name']        = self.spec_name.get()
        common['dark_flag']        = int(settings['dark_flag'])
        common['flat_flag']        = int(settings['flat_flag'])
        
        # Turn of dark flag if in real time and no darks have been taken
        if rt_flag == 'rt_analysis' and settings['rt_dark_flag'] == False:
            common['dark_flag'] = 0
            self.print_output('WARNING! No dark spectra aquired')

#========================================================================================
#================================Build parameter dictionary==============================
#========================================================================================
        
        # Create parameter array
        common['params'] = []

        for i in range(common['poly_n']):
            common['params'].append(('p'+str(i), 1.0))
            
        # Add other parameters
        common['params'].append(('shift',    float(self.shift.get())   ))
        common['params'].append(('stretch',  float(self.stretch.get()) ))
        common['params'].append(('ring_amt', float(self.ring_amt.get())))
        common['params'].append(('so2_amt',  float(self.so2_amt.get()) ))
        common['params'].append(('no2_amt',  float(self.no2_amt.get()) ))       
        common['params'].append(('o3_amt',   float(self.o3_amt.get())  ))
        common['params'].append(('bro_amt',  float(self.bro_amt.get()) ))
        common['last_spec'] = np.array([0])
        
        # Include optional paramters
        if common['ils_flag'] == True:
            common['params'].append(('ils_width', float(self.ils_width_e.get())))
            
        if common['ldf_flag'] == True:
            common['params'].append(('ldf', float(self.ldf_e.get())))

        # Change to ordered dictionary
        common['params'] = OrderedDict(common['params'])
        
        # Save initial guess parameters
        common['initial_params'] = common['params']

#========================================================================================
#=================================Read in xsecs and flat=================================
#========================================================================================

        # Build filepath to flat spectrum from spectrometer serial number
        if rt_flag == 'post_analysis':
            settings['flat_path'] = 'data_bases/flat_'+str(self.spec_name.get())+'.txt'
        else:
            settings['flat_path'] = 'data_bases/flat_'+str(self.c_spec.get())+'.txt'
        
        # Load fitting data files
        common = build_fwd_data(self, common, settings)

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
            common['dark'] = settings['dark_spec']

#========================================================================================
#========================Read test spectrum to get wavelength grid=======================
#========================================================================================
 
        if rt_flag == 'post_analysis':
            
            # Read first spectrum to get date of data and define stray light indices
            x,y,read_date,read_time,spec_no = read_spectrum(common['spectra_files'][0],
                                                            spec_type)
            
        else:

            # Read a single spectrum to get wavelength data
            try:
                x, y, header, t = aquire_spectrum(self, settings['spec'], 1, 1)
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
        
        # Create flag to control whether stray light is removed
        common['stray_flag'] = True
        
        # If no stray light pixels available, turn off the flag
        if common['stray_i1'] == common['stray_i2']:
            common['stray_flag'] = False 
 
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
                with open(out_excel_fname, 'w') as writer:
                    
                    # Write header line
                    writer.write('File,Number,Date,Time,so2 (ppm.m),so2 error,')
                    
                    for i in common['params'].keys():
                        writer.write(i + ',' + i + '_e,')
                        
                    # Write other fit info
                    writer.write('Fit window: ' + str(common['wave_start']) + ' - ' + \
                                 str(common['wave_stop']) + ' nm,' + 'ILS width: '  + \
                                 str(common['ils_width']) + 'ILS Gauss Weight: '    + \
                                 str(common['ils_gauss_weight']) + '\n')
                
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
            post_results_folder = 'Results/iFit/' + read_date + '/'
            
            # Create folder if it doesn't exist
            if not os.path.exists(post_results_folder):
                    os.makedirs(post_results_folder)
            
            # Create filename for output file
            out_excel_fname = post_results_folder + 'iFit_out.csv'
        
            try:
                
                # Open excel file and write header line
                with open(out_excel_fname, 'w') as writer:
                    
                    # Write header line
                    writer.write('File,Number,Date,Time,so2 (ppm.m),so2 error,')
                    
                    for i in common['params'].keys():
                        writer.write(i + ',' + i + '_e,')
                        
                    # Write other fit info
                    writer.write('Fit window: ' + str(common['wave_start']) + ' - ' + \
                                 str(common['wave_stop']) + ' nm,' + 'ILS width: '  + \
                                 str(common['ils_width']) + 'ILS Gauss Weight: '    + \
                                 str(common['ils_gauss_weight']) + '\n')
                    
            except PermissionError:
                self.print_output('Please close iFit output file to continue')
                return                

#========================================================================================
#===================================Start Analysis Loop==================================
#========================================================================================

        # Open excel file and write header line
        with open(out_excel_fname, 'a') as writer:

            # Print output
            self.print_output('Loop Started\n' +\
                              'Spectrum number ' + str(settings['loop']))  
            
            # Create empty arrays to hold the loop number and so2_amt values
            spec_nos = []
            so2_amts = []
            amt_errs = []

            # Begin analysis loop
            while True:
                
                # End loop if finished
                if settings['stop_flag'] == True:
                    break
                    
                # Read spectrum from file and fit
                if rt_flag == 'post_analysis':

                    try:
                        fname = common['spectra_files'][settings['loop']]
                        x,y,read_date,read_time,spec_no = read_spectrum(fname,spec_type)
                        
                        # Fit
                        fit_p,err_dict,y_data,so2_T,so2_spec,fit_flag=fit_spec(common,
                                                                               y,grid)
                        
                        now_fit_spec = True
                
                    except IndexError:
                        self.print_output('Fitting complete')
                        break 
                
                # Read spectrum from spectrometer and fit
                elif rt_flag=='rt_analysis' and self.toggle_button.config('text')[-1]==\
                'FITTING ON' and common['last_spec'].all() != 0:

                    # Create results queue
                    result_queue = Queue()
                    
                    # Create two threads, one to read a spectrum and one to fit
                    t1 = Thread(target = aquire_spectrum, args = (self,
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
                    fit_p,err_dict,y_data,so2_T,so2_spec,fit_flag=thread_out['fit']
                   
                    # Get spectrum
                    x, y, header, t = thread_out['spectrum']
                    
                    # Build file name
                    n = str('{num:05d}'.format(num=settings['loop']))
                    fname = settings['rt_folder'] + 'spectra/spectrum_' + n + '.txt'
                    
                    # Save
                    np.savetxt(fname, np.column_stack((x,y)), header = header)
                    
                    # Update last spec variable and spec number
                    common['last_spec'] = y
                    spec_no = settings['loop']
                    
                    now_fit_spec = True
                
                # Read spectrum from spectrometer but do not fit
                else:

                    # Read spectrum
                    x, y, header, t = aquire_spectrum(self, settings['spec'], 
                                                      settings['int_time'],
                                                      int(self.coadds.get()))
                    
                    # Build file name
                    n = str('{num:05d}'.format(num=settings['loop']))
                    fname = settings['rt_folder'] + 'spectra/spectrum_' + n + '.txt'
                    
                    # Save
                    np.savetxt(fname, np.column_stack((x,y)), header = header)
                    
                    # Update last spec variable
                    common['last_spec'] = y
                    
                    now_fit_spec = False
                
#========================================================================================
#=========================Analyse spectrum and add to output file========================
#========================================================================================
                 
                if now_fit_spec == True:
                    
                    # Unpack fit results
                    fit_dict  = {}
                    for m, l in enumerate(common['params'].keys()):
                        fit_dict[l] = fit_p[m]
    
                    # Write results to excel file, starting with spectrum info
                    writer.write(str(fname)     + ',' + \
                                 str(spec_no)   + ',' + \
                                 str(read_date) + ',' + \
                                 str(read_time))
                    
                    # Print so2 amount and error in ppm.m for ease
                    writer.write(',' + str(fit_dict['so2_amt']/2.463e15) + ',' + \
                                 str(err_dict['so2_amt']/2.463e15))            
        
                    # Print fit results and error for each parameter          
                    for l in common['params'].keys():
                        writer.write(','+str(fit_dict[l])+','+str(err_dict[l]))
                        
                    # Start new line
                    writer.write('\n')
                
                    # Add values to array for plotting
                    spec_nos.append(spec_no)
                    so2_amts.append(fit_dict['so2_amt']/2.463e15)
                    amt_errs.append(err_dict['so2_amt']/2.463e15)
                    
                    # Update quick analysis with values
                    last_so2 = "{0:0.1f}".format(fit_dict['so2_amt']/2.463e15)
                    last_err = "{0:0.1f}".format(err_dict['so2_amt']/2.463e15)
                    self.last_so2_amt.set(last_so2 + ' ppm.m')
                    self.last_so2_err.set(last_err + ' ppm.m')
                    
                    # Cut if too long to avoid slowing program
                    if settings['scroll_flag'] == True:
                        if len(spec_nos) > settings['scroll_spec_no']:
                            spec_nos = spec_nos[1:]
                            so2_amts = so2_amts[1:]
                            amt_errs = amt_errs[1:]
                
                    # Feed fit params into forward
                    fit = ifit_fwd(grid, *fit_p)
    
                    # Calculate the residual of the fit
                    resid = np.multiply(np.divide(np.subtract(y_data, fit), y_data), 100)
                    
                    # If fit fails or max resid > 10% revert to initial fit parameters
                    
                    if settings['update_params'] == True:
                        if fit_flag == False:
                            common['params'] = common['initial_params']
                            self.print_output('Fitting for spectrum ' + str(spec_no) + \
                                              ' failed, resetting parameters')
                            
                        elif max((resid)**2)**0.5 > 10:
                            common['params'] = common['initial_params']
                            self.print_output('Fitting for spectrum ' + str(spec_no) + \
                                              ' bad, resetting parameters')
                       
                        else:
                            
                            # Update first guesses with last fitted params
                            for i in fit_dict:
                                common['params'][i] = fit_dict[i]
                
#========================================================================================
#=======================================Update plot======================================
#========================================================================================
        
                # Replot data
                if int(settings['Show Graphs']) == 1:            
                    
                    if now_fit_spec == True:
                    
                        # Build axes and lines arrays
                        lines = [self.line0, self.line1, self.line2, self.line3, 
                                 self.line4, self.line5, self.line6]
                        axes =  [self.ax0,   self.ax0,   self.ax1,   self.ax2,  
                                 self.ax3,   self.ax3,   self.ax4  ]
                        
                        # Calculate graph limits for spectrum fit
                        y_lo = min(y_data) - abs((0.1*max(y_data)))
                        y_hi = max(y_data) + abs((0.1*max(y_data)))
                        f_lo = min(fit) - abs((0.1*max(fit)))
                        f_hi = max(fit) + abs((0.1*max(fit)))
                        y_lo, y_hi = min([f_lo, y_lo]), max([f_hi, y_hi])
                        x_lo, x_hi = grid.min() - 1, grid.max() + 1
                        
                        # Limits for residual
                        r_lo = min(resid) - abs((0.1*max(resid)))
                        r_hi = max(resid) + abs((0.1*max(resid)))
                        
                        # Limits for full spectrum
                        g_lo = min(x) - abs((0.1*max(x)))
                        g_hi = max(x) + abs((0.1*max(x)))
                        i_lo = 0 #min(y) - abs((0.1*max(y)))
                        i_hi = max(y) + abs((0.1*max(y))) 
                        
                        # Limits for so2 trans spectrum
                        t_lo = min(so2_T) - abs((0.1*min(so2_T)))
                        t_hi = max(so2_T) + (0.1*max(so2_T)) 
                        s_lo = min(so2_spec) - abs((0.1*min(so2_spec)))
                        s_hi = max(so2_spec) + (0.1*max(so2_spec))
                        t_lo, t_hi = min([t_lo, s_lo]), max([t_hi, s_hi])
                        
                        # Limits for so2 time series
                        a_lo = min(so2_amts) - abs((0.1*max(so2_amts)))
                        a_hi = max(so2_amts) + abs((0.1*max(so2_amts)))
                        
                        # Build data array to pass to graphing function
                        #                 x data    y data    x limits     y limits
                        data = np.array(([grid,     y_data,   [x_lo,x_hi], [y_lo,y_hi]],
                                         [grid,     fit,      [x_lo,x_hi], [y_lo,y_hi]],
                                         [x   ,     y,        [g_lo,g_hi], [i_lo,i_hi]],
                                         [grid,     resid,    [x_lo,x_hi], [r_lo,r_hi]],
                                         [grid,     so2_T,    [x_lo,x_hi], [t_lo,t_hi]],
                                         [grid,     so2_spec, [x_lo,x_hi], [t_lo,t_hi]],
                                         [spec_nos, so2_amts, False,       [a_lo,a_hi]]))
                    
                    else:
                        
                        # Build axes and lines arrays
                        lines = [self.line2]
                        axes =  [self.ax1  ]
                        
                        # Calculate limits
                        y_lo  = min(y) - abs((0.1*max(y)))
                        y_hi  = max(y) + abs((0.1*max(y)))
                        x_lo, x_hi = grid.min() - 1, grid.max() + 1
                        
                        # Build data array to pass to graphing function
                        #                 x data    y data    x limits     y limits
                        data = np.array(([x,        y,        [x_lo,x_hi], [y_lo,y_hi]]))
                    
                    # Update graph
                    update_graph(lines, axes, self, data)
                                               
                # Add to the count cycle
                settings['loop'] += 1
                
                # Force gui to update
                mygui.update(self)
                    



















#========================================================================================
#========================================================================================
#==================================Option menu commands==================================
#========================================================================================
#========================================================================================
  

#========================================================================================
#=====================================Model Settings=====================================
#========================================================================================
            
            
# Define function to open another window to alter the graph settings
def model_settings():
    popup = tk.Tk()
    tk.Tk.wm_title(popup, 'Model Settings')
    
    def update_model_settings(settings):
        
        # Update graph settings in common
        settings['model_resolution'] = float(model_res_e.get())
        settings['dark_flag'] = dark_b.get()
        settings['flat_flag'] = flat_b.get()
        settings['update_params'] = update_b.get()
        
        # Close the window
        popup.destroy()
    
    # Create row number counter
    row_n = 0
    
    # Set resolution of model grid
    popup.model_res = tk.IntVar(popup, value = settings['model_resolution'])
    model_res_l = tk.Label(popup, text = 'Model Grid spacing:', font = NORM_FONT)
    model_res_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    model_res_e = ttk.Entry(popup, textvariable = popup.model_res)
    model_res_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove dark spectra
    dark_b = tk.IntVar(popup, value = settings['dark_flag'])
    dark_l = tk.Label(popup, text = 'Remove Dark Spectra?', font = NORM_FONT)
    dark_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    dark_c = ttk.Checkbutton(popup, variable = dark_b)
    dark_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove flat spectra
    flat_b = tk.IntVar(popup, value = settings['flat_flag'])
    flat_l = tk.Label(popup, text = 'Remove Flat Spectra?', font = NORM_FONT)
    flat_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    flat_c = ttk.Checkbutton(popup, variable = flat_b)
    flat_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to update fit parameter guesses with the last fit values
    update_b = tk.IntVar(popup, value = settings['update_params'])
    update_l = tk.Label(popup, text = 'Auto-update fit paramseters?', font = NORM_FONT)
    update_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    update_c = ttk.Checkbutton(popup, variable = update_b)
    update_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Button to apply changes and close
    b1 = ttk.Button(popup, text='Apply', command=lambda: update_model_settings(settings))
    b1.grid(row = row_n, column = 0, padx = 5, pady = 5)
    b2 = ttk.Button(popup, text='Cancel', command=lambda: popup.destroy())
    b2.grid(row = row_n, column = 1, padx = 5, pady = 5) 
    row_n += 1
    
    
#========================================================================================
#====================================Graph Settings======================================
#========================================================================================


# Define function to open another window to alter the graph settings
def graph_settings():
    popup = tk.Tk()
    tk.Tk.wm_title(popup, 'Graph Settings')
    
    def update_graph_settings(settings):
        
        # Update graph settings in common
        settings['Show Graphs']     = graph_b.get()
        settings['Show Error Bars'] = err_b.get()
        settings['scroll_flag']     = scroll_b.get()
        settings['scroll_spec_no']  = int(spec_no_e.get())
        
        # Close the window
        popup.destroy()
    
    # Create row number counter
    row_n = 0
    
    # Control graph display settings
    graph_b = tk.IntVar(popup, value = settings['Show Graphs'])
    graph_l = tk.Label(popup, text = 'Show Graphs?', font = NORM_FONT)
    graph_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    graph_c = ttk.Checkbutton(popup, variable = graph_b)
    graph_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Turn on/off error bars
    err_b = tk.IntVar(popup, value = settings['Show Error Bars'])
    err_l = tk.Label(popup, text = 'Show Error Bars?', font = NORM_FONT)
    err_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    err_c = ttk.Checkbutton(popup, variable = err_b)
    err_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Turn on/off graph scrolling
    scroll_b = tk.IntVar(popup, value = settings['scroll_flag'])
    scroll_l = tk.Label(popup, text = 'Scroll Graph?', font = NORM_FONT)
    scroll_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    scroll_c = ttk.Checkbutton(popup, variable = scroll_b)
    scroll_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Set number of spectra to display on graph
    popup.spec_no = tk.IntVar(popup, value = settings['scroll_spec_no'])
    spec_no_l = tk.Label(popup, text = 'No. Spectra to display', font = NORM_FONT)
    spec_no_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    spec_no_e = ttk.Entry(popup, textvariable = popup.spec_no, font = NORM_FONT)
    spec_no_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Button to apply changes and close
    b1 = ttk.Button(popup, text='Apply', command=lambda: update_graph_settings(settings))
    b1.grid(row = row_n, column = 0, padx = 5, pady = 5)
    b2 = ttk.Button(popup, text='Cancel', command=lambda: popup.destroy())
    b2.grid(row = row_n, column = 1, padx = 5, pady = 5) 
    row_n += 1
    

#========================================================================================
#==================================Data base file paths==================================
#========================================================================================
 
    
# Define function to open another window to alter the data base files
def data_bases():
    popup = tk.Tk()
    tk.Tk.wm_title(popup, 'Data Base Settings')
    
    def update_fp(entry):
        
        # Open dialouge to get files
        fpath = fd.askopenfilenames()
        
        entry.set(str(fpath[0]))
    
    def update_data_bases(settings):
        
        # Update graph settings in common
        settings['sol_path']   = sol_path_e.get()
        settings['ring_path']  = ring_path_e.get()
        settings['so2_path']   = so2_path_e.get()
        settings['no2_path']   = no2_path_e.get()
        settings['o3_path']    = o3_path_e.get()
        settings['bro_path']   = bro_path_e.get()
        
        # Close the window
        popup.destroy()
    
    # Create row number counter
    row_n = 0
    
    # Create inputs for data base filepaths

    # Solar spectrum
    popup.sol_path = tk.StringVar(popup, value = settings['sol_path'])
    sol_path_l = tk.Label(popup, text = 'Solar spectrum:', font = NORM_FONT)
    sol_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    sol_path_e = ttk.Entry(popup, text = popup.sol_path, font = NORM_FONT, width = 40)
    sol_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    sol_path_b = ttk.Button(popup, text = "Select", 
                            command = lambda: update_fp(popup.sol_path))
    sol_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # Ring spectrum
    popup.ring_path = tk.StringVar(popup, value = settings['ring_path'])
    ring_path_l = tk.Label(popup, text = 'Ring Spectrum:', font = NORM_FONT)
    ring_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    ring_path_e = ttk.Entry(popup, textvariable = popup.ring_path, font = NORM_FONT, 
                            width = 40)
    ring_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    ring_path_b = ttk.Button(popup, text = "Select", 
                             command = lambda: update_fp(popup.ring_path))
    ring_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # SO2 xsec
    popup.so2_path = tk.StringVar(popup, value = settings['so2_path'])
    so2_path_l = tk.Label(popup, text = 'SO2 xsec:', font = NORM_FONT)
    so2_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    so2_path_e = ttk.Entry(popup, textvariable = popup.so2_path, font = NORM_FONT,
                           width = 40)
    so2_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    so2_path_b = ttk.Button(popup, text = "Select", 
                            command = lambda: update_fp(popup.so2_path))
    so2_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # NO2 xsec
    popup.no2_path = tk.StringVar(popup, value = settings['no2_path'])
    no2_path_l = tk.Label(popup, text = 'NO2 xsec:', font = NORM_FONT)
    no2_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    no2_path_e = ttk.Entry(popup, textvariable = popup.no2_path, font = NORM_FONT,
                           width = 40)
    no2_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    no2_path_b = ttk.Button(popup, text = "Select", 
                            command = lambda: update_fp(popup.no2_path))
    no2_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # O3 xsec
    popup.o3_path = tk.StringVar(popup, value = settings['o3_path'])
    o3_path_l = tk.Label(popup, text = 'O3 xsec:', font = NORM_FONT)
    o3_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    o3_path_e = ttk.Entry(popup, textvariable = popup.o3_path, font = NORM_FONT,
                          width = 40)
    o3_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    o3_path_b = ttk.Button(popup, text = "Select", 
                           command = lambda: update_fp(popup.o3_path))
    o3_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # Bro xsec
    popup.bro_path = tk.StringVar(popup, value = settings['bro_path'])
    bro_path_l = tk.Label(popup, text = 'BrO xsec:', font = NORM_FONT)
    bro_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    bro_path_e = ttk.Entry(popup, textvariable = popup.bro_path, font = NORM_FONT,
                           width = 40)
    bro_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5) 
    bro_path_b = ttk.Button(popup, text = "Select", 
                            command = lambda: update_fp(popup.bro_path))
    bro_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # Button to apply changes and close
    b1 = ttk.Button(popup, text='Apply', command=lambda: update_data_bases(settings))
    b1.grid(row = row_n, column = 0, padx = 5, pady = 5, columnspan=2)
    b2 = ttk.Button(popup, text='Cancel', command=lambda: popup.destroy())
    b2.grid(row = row_n, column = 1, padx = 5, pady = 5, columnspan=2) 
    row_n += 1
          
# Tkinter stuff      
mygui().mainloop()
