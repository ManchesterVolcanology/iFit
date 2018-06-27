#!/home/pi/berryconda3/bin/python3.6

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
from collections import OrderedDict
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from ifit_lib.build_fwd_data import build_fwd_data
from ifit_lib.fit import fit_spec, ifit_fwd
from ifit_lib.aquire_spectrum import aquire_spectrum
from ifit_lib.find_nearest import extract_window
from ifit_lib.file_control import make_directory

# Define some fonts to use in the program
NORM_FONT = ('Verdana', 8)
MED_FONT  = ('Veranda', 8)
LARG_FONT = ('Verdana', 9, 'bold')

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
        ttk.Style().configure('TEntry', width = 10) 
        
        # Add a title and icon
        tk.Tk.wm_title(self, 'piFit-2-2')
        #tk.Tk.iconbitmap('data_bases/icon.ico')
        
        # Create notebook to hold different frames
        nb = ttk.Notebook(self)
        page1 = ttk.Frame(nb)
        page2 = ttk.Frame(nb)
        
        # Create two frames, one for post analysis and one for real time acquisition
        nb.add(page1, text = 'Control')
        nb.add(page2, text =  'Model Parameters')
        nb.grid(column=0, row = 0, padx=5, pady=5, rowspan=2)
        
        # Create frame to hold graphs
        graph_frame = ttk.Frame(self, relief = 'groove')
        graph_frame.grid( column=1,row=0, padx=10, pady=10, sticky="NW")
        graph_frame.columnconfigure(index=0, weight=1)
        graph_frame.rowconfigure(index = 0, weight = 1)
        
        # Create frame to hold text output
        text_frame = ttk.Frame(self, relief = 'groove')
        text_frame.grid(column=1, row=1, padx=10, pady=10, sticky="NW")
        
        # Frame for quick analysis
        quick_frame = tk.Frame(page1, relief = 'groove')
        quick_frame.grid(column=0, row=2, padx=10, pady=10, sticky="NW")
        
        mygui.columnconfigure(index=1, weight=1, self = self)
        mygui.rowconfigure(index = 5, weight = 1, self = self)
        
#========================================================================================
#====================================Create text output==================================
#========================================================================================        
                 
        # Build text box
        self.text_box = tkst.ScrolledText(text_frame, width = 45, height = 5)
        self.text_box.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W',
                           columnspan = 4)
        self.text_box.insert('1.0', 'Welcome to piFit! Written by Ben Esse\n\n')
        
        # Create progress bar
        self.progress = ttk.Progressbar(text_frame, orient = tk.HORIZONTAL, length=300,
                                        mode = 'determinate')
        self.progress.grid(row = 0, column = 1, padx = 5, pady = 5, columnspan = 4)
        
        # Create status indicator
        self.status = tk.StringVar(text_frame, value = 'Standby')
        self.status_e = tk.Label(text_frame, textvariable = self.status)
        self.status_e.grid(row=0, column=0, padx=5, pady=5, sticky="EW")
        
#========================================================================================
#==============================Create quick analysis outputs=============================
#========================================================================================
        
        # Create ouput for last so2 amount
        self.last_so2_amt = tk.DoubleVar(self, value = 0)
        last_so2_amt_l = tk.Label(quick_frame, text = 'Last amt:', 
                                  font = NORM_FONT)
        last_so2_amt_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        last_so2_amt_e = ttk.Entry(quick_frame, textvariable = self.last_so2_amt,
                                   width = 10)
        last_so2_amt_e.grid(row = 0, column = 1, padx = 5, pady = 5)
        
        # Create ouput for last so2 error
        self.last_so2_err = tk.DoubleVar(self, value = 0)
        last_so2_err_l = tk.Label(quick_frame, text = 'Last error:', 
                                  font = NORM_FONT)
        last_so2_err_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        last_so2_err_e = ttk.Entry(quick_frame, textvariable = self.last_so2_err,
                                   width = 10)
        last_so2_err_e.grid(row = 1, column = 1, padx = 5, pady = 5)

        # Create button for advanced settings
        adv_set_b = ttk.Button(quick_frame, text = 'Adv. Settings', width = 15, 
                            command=lambda: adv_settings(self))
        adv_set_b.grid(row = 0, column = 2, padx = 5, pady = 5)

        # Create button to save settings
        save_b = ttk.Button(quick_frame, text = 'Save Settings', width = 15,
                            command = self.save)
        save_b.grid(row = 1, column = 2, padx = 5, pady = 5)
        
#========================================================================================
#===================================Set program settings=================================
#========================================================================================

        # Create settings dictionary
        global settings
        settings = {}
        
        # Read in settings file
        try:
            with open('data_bases/pifit_settings.txt', 'r') as r:
                
                # Read data
                data = r.readlines()
                
                # Unpack and save to dictionary
                for i in data:
                    name, val = i.strip().split(';')
                    settings[name] = val
   
        except FileNotFoundError:
            self.print_output('No settings file found, reverting to default')
            settings['Wave Start']        = 305
            settings['Wave Stop']         = 318
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
            settings['good_fit_bound']    = 10
            settings['Show Graphs']       = 1
            settings['Show Error Bars']   = 0
            settings['analysis_gas']      = 'SO2'
            settings['scroll_flag']       = 1
            settings['scroll_spec_no']    = 200
            settings['resid_type']        = 'Spec/Fit'
            settings['solar_resid_flag']  = 'Ignore'
            settings['poly_n']            = 4
            settings['shift']             = -0.2
            settings['stretch']           = 0.05
            settings['ring']              = 1.0
            settings['SO2']               = 1e+16
            settings['NO2']               = 1e+17
            settings['O3']                = 1e+19
            settings['BrO']               = 1e+15
            settings['model_resolution']  = 0.02
            settings['sol_path']          = 'data_bases/gas data/sao2010.txt'
            settings['ring_path']         = 'data_bases/gas data/ring.dat'
            settings['so2_path']          = 'data_bases/gas data/SO2_293K.dat'
            settings['no2_path']          = 'data_bases/gas data/No2_223l.dat'
            settings['o3_path']           = 'data_bases/gas data/O3_293K.dat'
            settings['bro_path']          = 'data_bases/gas data/BrO_Cross_298K.txt'
            settings['solar_resid_path']  = 'data_bases/gas data/solar_resid.txt'
            settings['graph_view']        = 'Spectrum'
 
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
        results_folder = 'Results/iFit/' + str(datetime.date.today()) + '/ifit_out/'
       
        # Create folder
        settings['folder'] = make_directory(results_folder, overwrite = True)
        
        # Create notes file
        settings['notes_fname'] = settings['folder'] + 'notes.txt'
        with open(settings['notes_fname'], 'w') as w:
            w.write('Notes file for iFit\n\n')
        
#========================================================================================
#====================================Create plot canvas==================================
#========================================================================================        
                    
        # Create figure to hold the graphs
        plt.rcParams.update({'font.size': 8} )
        self.fig = plt.figure(figsize = (4,2.5))
        
        # Create plot axes
        self.ax = self.fig.add_subplot(111)
        
        # Set axis labels
        self.ax.set_ylabel('SO2 (ppm.m)')
        self.ax.set_xlabel('Spectrum Number')
        
        # SO2 time series
        self.line0, = self.ax.plot(0, 0, lw = 1.5)
        self.line1, = self.ax.plot(0, 0, lw = 1.5)
        
        # Make it look nice
        plt.tight_layout()
        
        # Create the canvas to hold the graph in the GUI
        self.canvas = FigureCanvasTkAgg(self.fig, graph_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady = 10)       




#========================================================================================         
#========================================================================================
#================================Real Time Analysis Frame================================
#======================================================================================== 
#======================================================================================== 





#========================================================================================
#====================================Create GUI frames===================================
#========================================================================================

        # Create frames for diferent sections for post analysis     
        setup_frame2 = tk.LabelFrame(page1, text='Spectrometer Setup', 
                                     font = LARG_FONT)
        setup_frame2.grid(row=0, column=0, padx = 10, pady = 10, sticky="ew")
        
        button_frame2 = tk.LabelFrame(page1, text='Control', font = LARG_FONT)
        button_frame2.grid(row=1, column=0, padx = 10, pady = 10, sticky="ew")

#========================================================================================
#==================================Create control inputs=================================
#========================================================================================
        
        # Create label to display the spectrometer name
        self.c_spec = tk.StringVar(setup_frame2, value = 'None connected')
        c_spec_l = ttk.Label(setup_frame2, text="Device: ", font = NORM_FONT)
        c_spec_l.grid(row=0, column=0, pady=5, padx=5, sticky='W')
        c_spec_e = ttk.Entry(setup_frame2, textvariable = self.c_spec, width = 10)
        c_spec_e.grid(row = 0, column = 1, padx = 5, pady = 5)
        
        # Integration Time
        self.int_time = tk.IntVar(setup_frame2, value = settings['int_time'])
        int_time_vals = [1,2,3,4,5,10,20,30,40,50,100,200,300,400,
                         500,600,700,800,900,1000]
        int_time_l = tk.Label(setup_frame2, text = 'Int. Time:', font = NORM_FONT)
        int_time_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        int_time_sb = tk.Spinbox(setup_frame2, textvariable=self.int_time,
                                   values = int_time_vals, width = 10)
        int_time_sb.grid(row = 1, column = 1, padx = 5, pady = 5)
        self.int_time.set(settings['int_time'])
        
        # Coadds
        self.coadds = tk.IntVar(setup_frame2, value = settings['coadds'])
        coadd_vals=[1,2,3,4,5,10,20,30,40,50,100,200,300,400,500,1000]
        coadds_l = tk.Label(setup_frame2, text = 'Coadds:', font = NORM_FONT)
        coadds_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        coadds_sb = tk.Spinbox(setup_frame2, textvariable = self.coadds, 
                               values = coadd_vals, width = 10)
        coadds_sb.grid(row = 2, column = 1, padx = 5, pady = 5)
        self.coadds.set(settings['coadds'])
        
        # Number of darks to get
        self.no_darks = tk.IntVar(setup_frame2, value = settings['no_darks'])
        no_dark_vals = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        no_darks_l = tk.Label(setup_frame2, text = 'No. Darks:', font = NORM_FONT)
        no_darks_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
        no_darks_sb = tk.Spinbox(setup_frame2, textvariable = self.no_darks,
                                 values = no_dark_vals, width = 10)
        no_darks_sb.grid(row = 3, column = 1, padx = 5, pady = 5)
        self.no_darks.set(settings['no_darks'])
        
        # Create button to connect to spectrometer
        connect_spec_b = ttk.Button(setup_frame2, text = 'Connect', width = 10,
                                    command = self.connect_spec)
        connect_spec_b.grid(row = 0, column = 2, pady = 5)
        
        # Create button to update integration time
        update_int_time_b = ttk.Button(setup_frame2, text = 'Update', width = 10,
                                       command = self.update_int_time)
        update_int_time_b.grid(row = 1, column = 2, pady = 5)
        
        # Create button to read a single spectrum
        test_spec_b = ttk.Button(setup_frame2, text = 'Test Spec', width = 10,
                                 command = self.test_spec)
        test_spec_b.grid(row = 2, column = 2, pady = 5)
        
        # Create button to read darks
        read_darks_b = ttk.Button(setup_frame2, text = 'Read Darks', width = 10,
                                  command = self.read_darks)
        read_darks_b.grid(row = 3, column = 2, pady = 5)
        
#========================================================================================
#==============================Create start and stop buttons=============================
#========================================================================================         
        
        # Create button to start
        start_aq_b = ttk.Button(button_frame2, text = 'Begin!', command = self.begin,
                                width = 10)
        start_aq_b.grid(row = 0, column = 0, padx = 5, pady = 5, columnspan=2)
        
        # Create button to stop
        stop_aq_b = ttk.Button(button_frame2, text = 'Stop', command = self.stop,
                               width = 10)
        stop_aq_b.grid(row = 0, column = 1, padx = 5, pady = 5, columnspan=2)
        
        # Create switch to toggle fitting on or off
        self.toggle_button = tk.Button(button_frame2, text = 'FITTING OFF', 
                                       command = self.fit_toggle, width = 12, height = 1,
                                       bg = 'red', font = LARG_FONT)
        self.toggle_button.grid(row=1, column=2, padx=5, pady=5)
        
        # Define function to handle changes to the graph view option
        def update(val):
            
            # Get labels and set
            labels = self.graph_labels[self.graph_view.get()]
            self.ax.set_xlabel(labels[0])
            self.ax.set_ylabel(labels[1])
            
            # Force gui to update
            self.canvas.draw()
            
        # Create optionmenu to decide what to show on the graph
        graph_options = [settings['graph_view'],
                         'Spectrum', 
                         'Fit',
                         'Residual',
                         'Absorbance',
                         'Gas amount']
        
        # Create graph label dictionary
        self.graph_labels = {'Spectrum':   ['Wavelength (nm)', 'Intensity (arb)'], 
                             'Fit':        ['Wavelength (nm)', 'Intensity (arb)'],
                             'Residual':   ['Wavelength (nm)', 'Residual'],
                             'Absorbance': ['Wavelength (nm)', 'Absorbance'],
                             'Gas amount': ['Spectrum Number', 'Amount (ppm.m)']}
                
        self.graph_view = tk.StringVar(button_frame2, value = settings['graph_view'])
        graph_view_l = tk.Label(button_frame2, text = 'Graph View:', font = NORM_FONT)
        graph_view_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        graph_view_m = ttk.OptionMenu(button_frame2, self.graph_view, *graph_options,
                                      command = update)
        graph_view_m.grid(row = 1, column = 1, padx = 5, pady = 5)
        



#========================================================================================         
#========================================================================================
#=================================Model Parameters Frame=================================
#======================================================================================== 
#======================================================================================== 






#========================================================================================
#====================================Create GUI frames===================================
#========================================================================================

        # Frame for model settings
        model_frame = tk.LabelFrame(page2, text = 'Model Setup', font = LARG_FONT)
        model_frame.grid(row=0, column=0, padx=10, pady=6, sticky="NW")

        # Frame for parameters
        param_frame = tk.LabelFrame(page2, text = 'Fit Paramters', font = LARG_FONT)
        param_frame.grid(row=1, column=0, padx=10, pady=5, sticky="NW")

#========================================================================================
#=======================================Model Setup======================================
#========================================================================================

        # Create entry for start and stop wavelengths
        self.wave_start = tk.DoubleVar(model_frame, value = settings['Wave Start'])
        self.wave_start_l = tk.Label(model_frame, text = 'Wave Start:', font = NORM_FONT)
        self.wave_start_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.wave_start_e = ttk.Entry(model_frame, textvariable = self.wave_start,
                                      width = 10)
        self.wave_start_e.grid(row = 0, column = 1, padx = 5, pady = 5)
        
        self.wave_stop = tk.DoubleVar(model_frame, value = settings['Wave Stop'])
        self.wave_stop_l = tk.Label(model_frame, text = 'Wave Stop:', font = NORM_FONT)
        self.wave_stop_l.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.wave_stop_e = ttk.Entry(model_frame, textvariable = self.wave_stop, 
                                     width = 10)
        self.wave_stop_e.grid(row = 0, column = 3, padx = 5, pady = 5)
        
        # Instrument lineshape width
        self.ils_width = tk.DoubleVar(model_frame, value = settings['ILS Width'])
        self.ils_width_b = tk.BooleanVar(model_frame, value = settings['Fit ILS'])
        self.ils_width_l = tk.Label(model_frame, text = 'ILS Width:', font = NORM_FONT)
        self.ils_width_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.ils_width_e = ttk.Entry(model_frame, textvariable = self.ils_width,
                                      width = 10)
        self.ils_width_e.grid(row = 1, column = 1, padx = 5, pady = 5)
        self.ils_width_l2 = tk.Label(model_frame, text = 'Fit ILS?', font = NORM_FONT)
        self.ils_width_l2.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.ils_width_c = ttk.Checkbutton(model_frame, variable = self.ils_width_b)
        self.ils_width_c.grid(row = 1, column = 3, padx = 5, pady = 5)
        
        # ILS Gaussian weighting
        self.gauss_weight = tk.DoubleVar(model_frame, value = settings['Gauss Weight'])
        gauss_weight_l = tk.Label(model_frame, text = 'Gauss Wt:', 
                                  font = NORM_FONT)
        gauss_weight_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        gauss_weight_e = ttk.Entry(model_frame, textvariable = self.gauss_weight,
                                      width = 10)
        gauss_weight_e.grid(row = 2, column = 1, padx = 5, pady = 5)
        
        # Light dilution factor
        self.ldf = tk.DoubleVar(model_frame, value = settings['LDF'])
        self.ldf_b = tk.BooleanVar(model_frame, value = settings['Fit LDF'])
        self.ldf_l = tk.Label(model_frame, text = 'LDF:', font = NORM_FONT)
        self.ldf_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.ldf_e = ttk.Entry(model_frame, textvariable = self.ldf,
                                      width = 10)
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
        poly_n_vals = [1,2,3,4,5,6,7,8,9,10]
        poly_n_l = tk.Label(param_frame, text = 'Poly Order:', font = NORM_FONT)
        poly_n_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        poly_n_e = tk.Spinbox(param_frame, values = poly_n_vals, width = 10)
        poly_n_e.grid(row = 1, column = 1, padx = 5, pady = 5)
        self.poly_n.set(settings['poly_n'])
        
        # Spectrometer wavelength shift parameters
        self.shift = tk.DoubleVar(param_frame, value = settings['shift'])
        shift_l = tk.Label(param_frame, text = 'shift:', font = NORM_FONT)
        shift_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        shift_e = ttk.Entry(param_frame, textvariable = self.shift, width = 10)
        
        shift_e.grid(row = 2, column = 1, padx = 5, pady = 5)
        
        self.stretch = tk.DoubleVar(param_frame, value = settings['stretch'])
        stretch_l = tk.Label(param_frame, text = 'stretch:', font = NORM_FONT)
        stretch_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
        stretch_e = ttk.Entry(param_frame, textvariable = self.stretch, width = 10)
        stretch_e.grid(row = 3, column = 1, padx = 5, pady = 5)
        
        # Ring effect
        self.ring_amt = tk.DoubleVar(param_frame, value = settings['ring'])
        ring_amt_l = tk.Label(param_frame, text = 'Ring:', font = NORM_FONT)
        ring_amt_l.grid(row = 4, column = 0, padx = 5, pady = 5, sticky = 'W')
        ring_amt_e = ttk.Entry(param_frame, textvariable = self.ring_amt, width = 10)
        ring_amt_e.grid(row = 4, column = 1, padx = 5, pady = 5)
        
        # Gas amounts
        self.so2_amt = tk.DoubleVar(param_frame, value = settings['SO2'])
        so2_amt_l = tk.Label(param_frame, text = 'SO2:', font = NORM_FONT)
        so2_amt_l.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
        so2_amt_e = ttk.Entry(param_frame, textvariable = self.so2_amt, width = 10)
        so2_amt_e.grid(row = 1, column = 3, padx = 5, pady = 5)
        
        self.no2_amt = tk.DoubleVar(param_frame, value = settings['NO2'])
        no2_amt_l = tk.Label(param_frame, text = 'NO2:', font = NORM_FONT)
        no2_amt_l.grid(row = 2, column = 2, padx = 5, pady = 5, sticky = 'W')
        no2_amt_e = ttk.Entry(param_frame, textvariable = self.no2_amt, width = 10)
        no2_amt_e.grid(row = 2, column = 3, padx = 5, pady = 5)
        
        self.o3_amt = tk.DoubleVar(param_frame, value = settings['O3'])
        o3_amt_l = tk.Label(param_frame, text = 'O3:', font = NORM_FONT)
        o3_amt_l.grid(row = 3, column = 2, padx = 5, pady = 5, sticky = 'W')
        o3_amt_e = ttk.Entry(param_frame, textvariable = self.o3_amt, width = 10)
        o3_amt_e.grid(row = 3, column = 3, padx = 5, pady = 5)
       
        self.bro_amt = tk.DoubleVar(param_frame, value = settings['BrO'])
        bro_amt_l = tk.Label(param_frame, text = 'BrO:', font = NORM_FONT)
        bro_amt_l.grid(row = 4, column = 2, padx = 5, pady = 5, sticky = 'W')
        bro_amt_e = ttk.Entry(param_frame, textvariable = self.bro_amt, width = 10)
        bro_amt_e.grid(row = 4, column = 3, padx = 5, pady = 5)
        





        
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






#========================================================================================         
#========================================================================================
#=====================================Button Functions===================================
#======================================================================================== 
#========================================================================================
        
        
        
        


#========================================================================================
#=======================================Save Settings====================================
#========================================================================================
    
    # Function to save setting to the ifit_settings.txt file        
    def save(self):
        
        # Create or overright settings file
        with open('data_bases/pifit_settings.txt', 'w') as w:
            
            # Save each setting from the gui into settings
            settings['Wave Start']   = str(self.wave_start_e.get())   
            settings['Wave Stop']    = str(self.wave_stop_e.get())      
            settings['int_time']     = str(self.int_time.get())         
            settings['coadds']       = str(self.coadds.get())
            settings['no_darks']     = str(self.no_darks.get())
            settings['ILS Width']    = str(self.ils_width_e.get())
            settings['Gauss Weight'] = str(self.gauss_weight.get())
            settings['Fit ILS']      = str(self.ils_width_b.get())
            settings['LDF']          = str(self.ldf_e.get())
            settings['Fit LDF']      = str(self.ldf_b.get())
            settings['poly_n']       = str(self.poly_n.get())
            settings['shift']        = str(self.shift.get())
            settings['stretch']      = str(self.stretch.get())
            settings['ring']         = str(self.ring_amt.get())
            settings['SO2']          = str(self.so2_amt.get())
            settings['NO2']          = str(self.no2_amt.get())
            settings['O3']           = str(self.o3_amt.get())
            settings['BrO']          = str(self.bro_amt.get())
            
            # Add all of the settings dictionary
            for s in settings:
                w.write(s + ';' + str(settings[s]) + '\n')
                
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
            settings['int_time'] = float(self.int_time.get())
            settings['spec'].integration_time_micros(float(settings['int_time'])*1000)
            
            self.print_output('Integration time updated to '+str(settings['int_time']) +\
                              '\nSpectrum number ' + str(settings['loop']))
            
        except KeyError:
            self.print_output('No spectrometer conected')

#========================================================================================
#===============================Read a single test spectrum==============================
#========================================================================================

    def test_spec(self):
        
        # Update status
        self.status.set('Acquiring')
        mygui.update(self)
        
        x, y, header, read_time = aquire_spectrum(self, 
                                                  settings['spec'],
                                                  settings['int_time'], 
                                                  int(self.coadds.get()))
        
        # Display the test spectrum
        self.line0.set_data(x, y)
        self.line1.set_data(x[0], y[0])
        
        # Rescale the axes
        self.ax.relim()
        self.ax.autoscale_view() 
        
        self.canvas.draw()
        
        # Update status
        self.status.set('Standby')

#========================================================================================
#========================================Read Darks======================================
#========================================================================================

    def read_darks(self):
        
        if tkMessageBox.askyesno('Record Darks', 'Ready to begin measuring darks?'):
        
            # Update status
            self.status.set('Acquiring darks')
            mygui.update(self)
            
            # Update notes
            self.print_output('Begin reading darks')
            
            # Reset progress bar
            self.progress['mode'] = 'determinate'
            self.progress['value'] = 0
            
            # Create zero array to hold dark spectra
            dark = np.zeros(2048)
            
            # Define dark filepath
            dark_fp = settings['folder'] + 'dark/'
            
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
                
                # Update the progress bar
                prog = ((i+1)/dark_n) * 100
                self.progress['value'] = prog
                
                # Force gui to update
                mygui.update(self)
            
            # Divide by number of darks to get average
            settings['dark_spec'] = np.divide(dark, dark_n)
            
            # Display the test spectrum
            self.line0.set_data(x, settings['dark_spec'])
            self.line1.set_data(x[0], settings['dark_spec'][0])
            
            # Rescale the axes
            self.ax.relim()
            self.ax.autoscale_view() 
            
            self.canvas.draw()
            
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
            
            # Update status
            self.status.set('Standby')
        
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
    def begin(self):
        
        # Turn off stopping flag
        settings['stop_flag'] = False
       
        # Create common dictionary
        common = {}
        
        # Populate common with other data from the GUI
        common['wave_start']       = float(self.wave_start_e.get())
        common['wave_stop']        = float(self.wave_stop_e.get())
        common['poly_n']           = int(self.poly_n.get()) + 1
        common['ils_width']        = float(self.ils_width_e.get())
        common['ils_gauss_weight'] = float(self.gauss_weight.get())
        common['ldf']              = float(self.ldf_e.get())
        common['ils_flag']         = self.ils_width_b.get()
        common['ldf_flag']         = self.ldf_b.get()
        common['dark_flag']        = int(settings['dark_flag'])
        common['flat_flag']        = int(settings['flat_flag'])
        common['solar_resid_flag'] = settings['solar_resid_flag']

        # Turn of dark flag if in real time and no darks have been taken
        if settings['rt_dark_flag'] == False:
            common['dark_flag'] = 0
            self.print_output('WARNING! No dark spectra aquired!')

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
        initial_params = common['params'].copy()

#========================================================================================
#=================================Read in xsecs and flat=================================
#========================================================================================

        # Update status
        self.status.set('Building Model')
        mygui.update(self)

        # Build filepath to flat spectrum from spectrometer serial number
        settings['flat_path'] = 'data_bases/Spectrometer/flat_' + \
                                 str(self.c_spec.get()) + '.txt'
        
        # Load fitting data files
        common = build_fwd_data(common, settings, self)

#========================================================================================
#===================================Build dark spectrum==================================
#========================================================================================

        if common['dark_flag'] == True:
            
            # Assign dark spectrum
            common['dark'] = settings['dark_spec']

#========================================================================================
#========================Read test spectrum to get wavelength grid=======================
#========================================================================================
 
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
            
        # If forming solar residual create empty array and counter
        if common['solar_resid_flag'] == 'Generate':
            common['solar_resid'] = np.zeros(len(grid))
            resid_count = 0
 
#========================================================================================
#===================================Create ouput folder==================================
#========================================================================================
            
        # Create new output folder
        if settings['create_out_flag'] == True:
           
            # Reset loop counter
            settings['loop'] = 0
            
            # Create filename for output file
            out_excel_fname = settings['folder'] + 'iFit_out.csv'
            
            # Open excel file and write header line
            with open(out_excel_fname, 'w') as writer:
                
                # Write header line
                writer.write('File,Number,Date,Time,so2 (ppm.m),so2 error,')
                
                for i in common['params'].keys():
                    writer.write(i + ',' + i + '_e,')
                    
                # Add column for fit info
                writer.write('Fit Quality,')
                    
                # Write other fit info
                writer.write('Fit window: ' + str(common['wave_start']) + ' - ' + \
                             str(common['wave_stop']) + ' nm,' + 'ILS width: '  + \
                             str(common['ils_width']) + ',ILS Gauss Weight: '   + \
                             str(common['ils_gauss_weight']) + '\n')
            
            # Create folder to hold spectra
            if not os.path.exists(settings['folder'] + 'spectra/'):
                    os.makedirs(settings['folder'] + 'spectra/')
            
            # Turn off flag to limit folders to one per program run
            settings['create_out_flag'] = False     
            
        else:
            out_excel_fname = settings['folder'] + 'iFit_out.csv'
            
#========================================================================================
#===========================Set Progress bar to correct format===========================
#========================================================================================
            
        self.progress['mode'] = 'indeterminate'
        self.progress['value'] = 0
             
#========================================================================================
#===================================Start Analysis Loop==================================
#========================================================================================

        # Open excel file and write header line
        with open(out_excel_fname, 'a') as writer:

            # Print output
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
            self.status.set('Acquiring')
            mygui.update(self)

            # Begin analysis loop
            while True:
                            
                # End loop if finished
                if settings['stop_flag'] == True:
                    break 
                        
                
                # Read spectrum from spectrometer and fit
                if self.toggle_button.config('text')[-1]=='FITTING ON' and \
                common['last_spec'].all() != 0:

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
                    fit_p, err_dict, y_data, gas_T, fit_flag = thread_out['fit']
                   
                    # Get spectrum
                    x, y, header, t = thread_out['spectrum']
                    read_date, read_time = t.split(' ')
                    
                    # Build file name
                    n = str('{num:05d}'.format(num=settings['loop']))
                    fname = settings['folder'] + 'spectra/spectrum_' + n + '.txt'
                    
                    # Save
                    np.savetxt(fname, np.column_stack((x, y)), header = header)
                    
                    # Update last spec variable and spec number
                    common['last_spec'] = y
                    spec_no = settings['loop']  
                    
                    # Update progress bar
                    prog = settings['loop']
                    self.progress['value'] = prog
                    
                    wait_flag = False
                
#========================================================================================
#====================================Just read spectrum==================================
#========================================================================================                
                                
                # Read spectrum from spectrometer but do not fit
                else:

                    # Read spectrum
                    x, y, header, t = aquire_spectrum(self, settings['spec'], 
                                                      settings['int_time'],
                                                      int(self.coadds.get()))
                    read_date, read_time = t.split(' ')
                    
                    # Build file name
                    n = str('{num:05d}'.format(num=settings['loop']))
                    fname = settings['folder'] + 'spectra/spectrum_' + n + '.txt'
                    
                    # Save
                    np.savetxt(fname, np.column_stack((x,y)), header = header)
                    
                    # Update last spec variable
                    common['last_spec'] = y
                    
                    # Update progress bar
                    prog = settings['loop']
                    self.progress['value'] = prog
                    
                    wait_flag = True
                
#========================================================================================
#========================Unpack fit results and add to output file=======================
#========================================================================================
                 
                if self.toggle_button.config('text')[-1]=='FITTING ON' and \
                wait_flag == False:
                    
                    # Unpack fit results
                    fit_dict  = {}
                    for m, l in enumerate(common['params'].keys()):
                        fit_dict[l] = fit_p[m]
    
                     # Feed fit params into forward
                    fit = ifit_fwd(grid, *fit_p)
    
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
                            for i in fit_dict:
                                common['params'][i] = fit_dict[i]
                            fit_msg = 'Good'
    
                    # Write results to excel file, starting with spectrum info
                    writer.write(str(fname)          + ',' + \
                                 str(spec_no)        + ',' + \
                                 str(str(read_date)) + ',' + \
                                 str(str(read_time)))
                    
                    # Print so2 amount and error in ppm.m for ease
                    writer.write(',' + str(fit_dict['so2_amt']/2.463e15) + ',' + \
                                 str(err_dict['so2_amt']/2.463e15))            
        
                    # Print fit results and error for each parameter          
                    for l in common['params'].keys():
                        writer.write(','+str(fit_dict[l])+','+str(err_dict[l]))
                        
                    # Write fit quality and start new line
                    writer.write(',' + fit_msg + '\n')
                
                    # Add values to array for plotting
                    spec_nos.append(spec_no)
                    gas['SO2_amts'].append(fit_dict['so2_amt']/2.463e15)
                    gas['SO2_errs'].append(err_dict['so2_amt']/2.463e15)
                    gas['O3_amts'].append(fit_dict['o3_amt']/2.463e15)
                    gas['O3_errs'].append(err_dict['o3_amt']/2.463e15)
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
                if int(settings['Show Graphs']) == 1 and wait_flag == False:            
                    
                    if self.toggle_button.config('text')[-1]=='FITTING ON':
                        
                        # Get selected transmittance data
                        gas_tran = gas_T[settings['analysis_gas'] + '_tran']
                        gas_spec = gas_T[settings['analysis_gas'] + '_spec']
                        gas_amts = gas[settings['analysis_gas'] + '_amts']
                        
                        graph_view = self.graph_view.get()
                        
                        if graph_view == 'Spectrum':
                            self.line0.set_data(x, y)
                            self.line1.set_data(x[0], y[0])
                            
                        if graph_view == 'Fit':
                            self.line0.set_data(grid, y_data)
                            self.line1.set_data(grid, fit)
                            
                        if graph_view == 'Residual':
                            self.line0.set_data(grid, resid)
                            self.line1.set_data(grid[0], resid[0])
                            
                        if graph_view == 'Absorbance':
                            self.line0.set_data(grid, gas_tran)
                            self.line1.set_data(grid, gas_spec)
                            
                        if graph_view == 'Gas amount':
                            self.line0.set_data(spec_nos, gas_amts)
                            self.line1.set_data(spec_nos[0], gas_amts[0])
                             
                    
                    else:
                        
                        # Just display the spectrum
                        self.line0.set_data(x, y)
                        self.line1.set_data(x[0], y[0])
                    
                    # Rescale the axes
                    self.ax.relim()
                    self.ax.autoscale_view() 
                    
                    # Update canvas
                    self.canvas.draw()
                               
                # Add to the count cycle
                settings['loop'] += 1
                
                # Force gui to update
                mygui.update(self)
                
            # Update solar residual
            if settings['solar_resid_flag'] == 'Generate':
                
                # Find average residual
                common['solar_resid'] = np.divide(common['solar_resid'], resid_count)
                
                # Save
                np.savetxt('data_bases/gas data/solar_resid.txt', common['solar_resid'])
                
                self.print_output('Solar residual spectrum updated')
                
            # Update status
            self.status.set('Standby')






#========================================================================================
#========================================================================================
#====================================Advanced Setings====================================
#========================================================================================
#========================================================================================


def adv_settings(self):
    
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
        settings['Show Error Bars']  = err_b.get()
        settings['scroll_flag']      = scroll_b.get()
        settings['scroll_spec_no']   = int(spec_no_e.get())
        settings['resid_type']       = resid_type.get()
        settings['analysis_gas']     = gas.get()
        
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
 
    # Create frames for diferent sections for post analysis and buttons
    model_frame = tk.LabelFrame(popup, text = 'Model Settings', font = LARG_FONT)
    model_frame.grid(row=0, column=0, padx = 10, pady = 10, sticky="NW")
    
    datab_frame = tk.LabelFrame(popup, text='Data Base Settings', font = LARG_FONT)
    datab_frame.grid(row=0, column=1, padx = 10, pady = 10, sticky="NW", rowspan=2)
    
    graph_frame = tk.LabelFrame(popup, text='Graph Settings', font = LARG_FONT)
    graph_frame.grid(row=1, column=0, padx = 10, pady = 10, sticky="NW")
    
    button_frame = ttk.Frame(popup)
    button_frame.grid(row=2, column=0, padx = 10, pady = 10, sticky="EW", columnspan = 2)

    # Button to apply changes and close
    b1 = ttk.Button(button_frame, text='Ok', command=lambda: update_settings(settings, 
                                                                             True))
    b1.grid(row = 0, column = 0, padx = 100, pady = 5, sticky = "EW")
    
    # Button to just apply changes
    b2 = ttk.Button(button_frame, text='Apply', command=lambda: update_settings(settings,
                                                                                False))
    b2.grid(row = 0, column = 1, padx = 100, pady = 5, sticky = "EW")
    
    # Buttong to cancel
    b3 = ttk.Button(button_frame, text='Cancel', command=lambda: popup.destroy())
    b3.grid(row = 0, column = 2, padx = 100, pady = 5, sticky = "EW")

#========================================================================================
#=====================================Model Settings=====================================
#========================================================================================

    # Create row number counter
    row_n = 0
    
    # Set resolution of model grid
    popup.model_res = tk.DoubleVar(model_frame, value = settings['model_resolution'])
    model_res_l = tk.Label(model_frame, text = 'Model Grid spacing:', font = NORM_FONT)
    model_res_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    model_res_e = ttk.Entry(model_frame, textvariable = popup.model_res)
    model_res_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove dark spectra
    dark_b = tk.IntVar(model_frame, value = settings['dark_flag'])
    dark_l = tk.Label(model_frame, text = 'Remove Dark Spectra?', font = NORM_FONT)
    dark_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    dark_c = ttk.Checkbutton(model_frame, variable = dark_b)
    dark_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove flat spectra
    flat_b = tk.IntVar(model_frame, value = settings['flat_flag'])
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
    update_b = tk.StringVar(model_frame, value = settings['update_params'])
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
    
    # Turn on/off error bars
    err_b = tk.IntVar(graph_frame, value = settings['Show Error Bars'])
    err_l = tk.Label(graph_frame, text = 'Show Error Bars?', font = NORM_FONT)
    err_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    err_c = ttk.Checkbutton(graph_frame, variable = err_b)
    err_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
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
    
# Run the App!
if __name__ == '__main__':    
    mygui().mainloop()
