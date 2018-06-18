# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 09:24:05 2018

@author: mqbpwbe2
"""

# Import required libraries
import matplotlib
matplotlib.use('TkAgg')
import os
import numpy as np
from tkinter import ttk
import tkinter as tk
from tkinter import filedialog as fd
import tkinter.scrolledtext as tkst
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import OrderedDict
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

from ifit_lib.build_fwd_data import build_fwd_data
from ifit_lib.fit import fit_spec, ifit_fwd
from ifit_lib.read_binary_block import read_binary_block
from ifit_lib.find_nearest import extract_window
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
        self.protocol("WM_DELETE_WINDOW",self.handler)
        
        # Button Style
        ttk.Style().configure('TButton', width = 20, height = 20, relief="flat") 
        
        # Add a title and icon
        tk.Tk.wm_title(self, 'iFit-2-2 Scanner')
        tk.Tk.iconbitmap(self, default = 'data_bases/icon.ico')

        # Build a menubar to hold options for the user
        menubar = tk.Menu(self)
        filemenu = tk.Menu(menubar, tearoff = 0)
        filemenu.add_command(label = 'Model Settings', command = model_settings)
        filemenu.add_command(label='Graph Settings',command=lambda: graph_settings(self))
        filemenu.add_command(label = 'Data Bases', command = data_bases)
        filemenu.add_separator()
        menubar.add_cascade(label = 'Options', menu = filemenu)
        tk.Tk.config(self, menu = menubar)

        # Create Frame to hold setup
        setup_frame = ttk.Frame(self, relief = 'groove')
        setup_frame.grid(row=0, column=0, padx=10, pady=10, sticky="NW")
        
        # Create frame to hold graphs
        graph_frame = ttk.Frame(self, relief = 'groove')
        graph_frame.grid(row=0, column=1, padx=10, pady=10, rowspan=10, sticky="NW")
        graph_frame.columnconfigure(index=0, weight=1)
        graph_frame.rowconfigure(index = 0, weight = 1)
        
        # Create frame to hold text output
        text_frame = ttk.Frame(self, relief = 'groove')
        text_frame.grid(row=2, column=0, padx=10, pady=10, sticky="NW")
        
        mygui.columnconfigure(index=1, weight=1, self = self)
        mygui.rowconfigure(index = 5, weight = 1, self = self)

#========================================================================================
#====================================Create text output==================================
#========================================================================================        
                 
         # Build text box
        self.text_box = tkst.ScrolledText(text_frame, width = 60, height = 10)
        self.text_box.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W',
                           columnspan = 4)
        self.text_box.insert('1.0', 'Welcome to iFit! Written by Ben Esse\n\n')    
        
        # Create button to start
        start_b = ttk.Button(text_frame, text = 'Begin!', command = self.begin)
        start_b.grid(row = 0, column = 0, padx = 5, pady = 5)
        
        # Create button to stop
        stop_b = ttk.Button(text_frame, text = 'Stop', command = self.stop)
        stop_b.grid(row = 0, column = 1, padx = 5, pady = 5)
        
        # Create button to save settings
        save_b = ttk.Button(text_frame, text = 'Save Settings', command = self.save)
        save_b.grid(row = 0, column = 2,  padx = 5, pady = 5)
        
        # Create a button to exit
        exit_b = ttk.Button(text_frame, text = 'Exit', command = self.quit)
        exit_b.grid(row = 0, column = 3, padx = 5, pady = 5)
        
#========================================================================================
#===================================Set program settings=================================
#========================================================================================

        # Create settings dictionary
        global settings
        settings = {}
        
        # Read in settings file
        try:
            with open('data_bases/ifit_scan_settings.txt', 'r') as r:
                
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
            settings['sol_path']          = 'data_bases/sao2010_full.txt'
            settings['ring_path']         = 'data_bases/ring.dat'
            settings['so2_path']          = 'data_bases/SO2_293K.dat'
            settings['no2_path']          = 'data_bases/No2_223l.dat'
            settings['o3_path']           = 'data_bases/O3_293K.dat'
            settings['bro_path']          = 'data_bases/BrO_Cross_298K.txt'
            settings['solar_resid_path']  = 'data_bases/solar_resid.txt'
            settings['Scan Filepaths'] = ''
 
        # Create loop counter to keep track of the analysis
        settings['loop'] = 0
        
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
        
        self.ax3.set_ylabel(settings['analysis_gas'] + ' Transmittance', fontsize = 10)
        self.ax3.set_xlabel('Wavelength (nm)', fontsize=10)
        
        self.ax4.set_ylabel(settings['analysis_gas'] + ' amt (ppm.m)', fontsize = 10)
        self.ax4.set_xlabel('Spcan Angle', fontsize=10)
        
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
        self.line5, = self.ax3.plot(0, 0, 'r', label = 'Gas trans')
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
#====================================Create GUI frames===================================
#========================================================================================

        # Frame for model settings
        model_frame = tk.LabelFrame(setup_frame, text = 'Model Setup', font = LARG_FONT)
        model_frame.grid(row=1, column=0, padx=10, pady=6, sticky="NW")

        # Frame for parameters
        param_frame = tk.LabelFrame(setup_frame, text = 'Fit Paramters',font = LARG_FONT)
        param_frame.grid(row=2, column=0, padx=10, pady=5, sticky="NW")
        
        # Frame for quick analysis
        quick_frame = tk.LabelFrame(setup_frame, text = 'Quick Analysis', font = LARG_FONT)
        quick_frame.grid(row=3, column=0, padx=10, pady=10, sticky="NW")
        
#========================================================================================
#==================================Create file dialouges=================================
#========================================================================================        
        
        # Reformat filepath strings
        self.scan_fpaths = []
        for i in settings['Scan Filepaths'][1:-1].split(', '):
            self.scan_fpaths.append(str(i[1:-1]))
        
        # File dialouge for spectra
        if self.scan_fpaths == ['']:
            message = 'No scans selected'
        else:
            message = str(len(self.scan_fpaths)) + ' scans selected'
        self.scan_ent = tk.StringVar(value = message)
        self.scanfp_l = tk.Entry(model_frame, font = NORM_FONT, width = 40, 
                                 text = self.scan_ent)
        self.scanfp_l.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = 'W', 
                      columnspan = 3)
        scanfp_b = ttk.Button(model_frame, text="Select Scans", command = self.scan_fp)
        scanfp_b.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')

#========================================================================================
#=======================================Model Setup======================================
#========================================================================================

        # Create entry for start and stop wavelengths
        self.wave_start = tk.DoubleVar(model_frame, value = settings['Wave Start'])
        self.wave_start_l = tk.Label(model_frame, text = 'Wave Start:', font = NORM_FONT)
        self.wave_start_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.wave_start_e = ttk.Entry(model_frame, textvariable = self.wave_start)
        self.wave_start_e.grid(row = 1, column = 1, padx = 5, pady = 5)
        
        self.wave_stop = tk.DoubleVar(model_frame, value = settings['Wave Stop'])
        self.wave_stop_l = tk.Label(model_frame, text = 'Wave Stop:', font = NORM_FONT)
        self.wave_stop_l.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.wave_stop_e = ttk.Entry(model_frame, textvariable = self.wave_stop)
        self.wave_stop_e.grid(row = 1, column = 3, padx = 5, pady = 5)
        
        # Instrument lineshape width
        self.ils_width = tk.DoubleVar(model_frame, value = settings['ILS Width'])
        self.ils_width_b = tk.BooleanVar(model_frame, value = settings['Fit ILS'])
        self.ils_width_l = tk.Label(model_frame, text = 'ILS Width:', font = NORM_FONT)
        self.ils_width_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.ils_width_e = ttk.Entry(model_frame, textvariable = self.ils_width)
        self.ils_width_e.grid(row = 2, column = 1, padx = 5, pady = 5)
        self.ils_width_l2 = tk.Label(model_frame, text = 'Fit ILS?', font = NORM_FONT)
        self.ils_width_l2.grid(row = 2, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.ils_width_c = ttk.Checkbutton(model_frame, variable = self.ils_width_b)
        self.ils_width_c.grid(row = 2, column = 3, padx = 5, pady = 5)
        
        # ILS Gaussian weighting
        self.gauss_weight = tk.DoubleVar(model_frame, value = settings['Gauss Weight'])
        gauss_weight_l = tk.Label(model_frame, text = 'ILS Gauss Weight:', 
                                  font = NORM_FONT)
        gauss_weight_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
        gauss_weight_e = ttk.Entry(model_frame, textvariable = self.gauss_weight)
        gauss_weight_e.grid(row = 3, column = 1, padx = 5, pady = 5)
        
        # Light dilution factor
        self.ldf = tk.DoubleVar(model_frame, value = settings['LDF'])
        self.ldf_b = tk.BooleanVar(model_frame, value = settings['Fit LDF'])
        self.ldf_l = tk.Label(model_frame, text = 'LDF:', font = NORM_FONT)
        self.ldf_l.grid(row = 4, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.ldf_e = ttk.Entry(model_frame, textvariable = self.ldf)
        self.ldf_e.grid(row = 4, column = 1, padx = 5, pady = 5)
        self.ldf_l2 = tk.Label(model_frame, text = 'Fit LDF?', font = NORM_FONT)
        self.ldf_l2.grid(row = 4, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.ldf_c = ttk.Checkbutton(model_frame, variable = self.ldf_b)
        self.ldf_c.grid(row = 4, column = 3, padx = 5, pady = 5)

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
#==============================Create quick analysis outputs=============================
#========================================================================================
        
        # Create loop counter for spectra
        self.spec_count = tk.StringVar(self, value = '1 / 1')
        spec_count_l = tk.Label(quick_frame, text = 'Spectrum:', font = NORM_FONT)
        spec_count_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        spec_count_e = ttk.Entry(quick_frame, textvariable = self.spec_count)
        spec_count_e.grid(row = 0, column = 1, padx = 5, pady = 5)
        
        # Create loop counter for scans
        self.scan_count = tk.StringVar(self, value = '1 / 1')
        scan_count_l = tk.Label(quick_frame, text = 'Scan:', font = NORM_FONT)
        scan_count_l.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = 'W')
        scan_count_e = ttk.Entry(quick_frame, textvariable = self.scan_count)
        scan_count_e.grid(row = 0, column = 3, padx = 5, pady = 5)
        
        # Create ouput for last so2 amount
        self.last_so2_amt = tk.DoubleVar(self, value = 0)
        last_so2_amt_l = tk.Label(quick_frame, text = 'Last amt:', 
                                  font = NORM_FONT)
        last_so2_amt_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        last_so2_amt_e = ttk.Entry(quick_frame, textvariable = self.last_so2_amt)
        last_so2_amt_e.grid(row = 1, column = 1, padx = 5, pady = 5)
        
        # Create ouput for last so2 error
        self.last_so2_err = tk.DoubleVar(self, value = 0)
        last_so2_err_l = tk.Label(quick_frame, text = 'Last error:', 
                                  font = NORM_FONT)
        last_so2_err_l.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
        last_so2_err_e = ttk.Entry(quick_frame, textvariable = self.last_so2_err)
        last_so2_err_e.grid(row = 1, column = 3, padx = 5, pady = 5)



        
#========================================================================================         
#========================================================================================
#=====================================Button Functions===================================
#======================================================================================== 
#======================================================================================== 

    def handler(self):

        self.quit()

#========================================================================================
#=====================================Filepath Buttons===================================
#========================================================================================
   
    # Function to select spectra to analyse
    def scan_fp(self):
            
        # Open dialouge to get files
        fpaths = fd.askopenfilenames()

        if fpaths != '':
            self.scan_fpaths = []
            for i in fpaths:
                self.scan_fpaths.append(str(i))
            
            # Save output to input 
            self.scan_ent.set(str(len(self.scan_fpaths)) + ' scans selected')
             
#========================================================================================
#=======================================Save Settings====================================
#========================================================================================
    
    # Function to save setting to the ifit_settings.txt file        
    def save(self):
        
        # Create or overright settings file
        with open('data_bases/ifit_scan_settings.txt', 'w') as w:
            
            # Save each setting into the file
            w.write('Wave Start;'       + str(self.wave_start_e.get())      + '\n')
            w.write('Wave Stop;'        + str(self.wave_stop_e.get())       + '\n')
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
            w.write('analysis_gas;'     + str(settings['analysis_gas'])     + '\n')
            w.write('resid_type;'       + str(settings['resid_type'])       + '\n')
            w.write('solar_resid_flag;' + str(settings['solar_resid_flag']) + '\n')
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
            w.write('solar_resid_path;' + str(settings['solar_resid_path']) + '\n')
            
            try:
                w.write('Scan Filepaths;' + str(self.scan_fpaths) + '\n')
            except AttributeError:
                w.write('Scan Filepaths; \n') 
                
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
        
        # Force gui to update
        mygui.update(self)    
        
#========================================================================================
#===========================================Stop=========================================
#========================================================================================
        
    def stop(self):
        settings['stop_flag'] = True
        self.print_output('Loop Stopped\nSpectrum number ' + str(settings['loop']-1))
        
#========================================================================================
#========================================Begin iFit======================================
#========================================================================================
   
    # Function to begin analysis loop
    def begin(self):
        
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
        common['scan_files']       = self.scan_fpaths
        common['dark_flag']        = bool(settings['dark_flag'])
        common['flat_flag']        = bool(settings['flat_flag'])
        common['solar_resid_flag'] = settings['solar_resid_flag']

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
#========================Read test spectrum to get wavelength grid=======================
#========================================================================================
        
        # Read in first scan in array
        err,x,header,info_block,spec_block = read_binary_block(common['scan_files'][0])
        
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
#=================================Read in xsecs and flat=================================
#========================================================================================

        # Build filepath to flat spectrum from spectrometer serial number
        settings['flat_path'] = 'data_bases/flat_'#+str(self.c_spec.get())+'.txt'
        
        # Load fitting data files
        common = build_fwd_data(common, settings, self)
 
#========================================================================================
#===================================Create ouput folder==================================
#========================================================================================
            
        # Reset loop counter
        settings['loop'] = 1
        
        # Get date
        meas_date = common['scan_files'][0].split('/')[-1][:6]

        # Create filepath to directory to hold program outputs
        results_folder = 'Results/iFit_scan/' + str(meas_date) + '/'
        
        # Create folder if it doesn't exist
        if not os.path.exists(results_folder):
                os.makedirs(results_folder)                     

#========================================================================================
#===================================Start Analysis Loop==================================
#========================================================================================

        # Initiate loop counter display
        msg = '1 / ' + str(len(common['scan_files']))
        self.scan_count.set(msg)

        # Output initiation message
        self.print_output('Loop Started\nSpectrum number ' + str(settings['loop']))

        for fpath in common['scan_files']:
            
            # Create empty arrays to hold the loop number and so2_amt values
            gas = {}
            spec_nos = []
            gas['SO2_amts'] = []
            gas['SO2_errs'] = []
            gas['O3_amts']  = []
            gas['O3_errs']  = []
            gas['BrO_amts'] = []
            gas['BrO_errs'] = []
            
            # End loop if finished
            if settings['stop_flag'] == True:
                break

            # Read in scan block
            err, x, header, info_block, spec_block = read_binary_block(fpath)
            
            if err == 0:
            
                # First spectrum is the dark
                common['dark'] = spec_block[:,0]
                
                # Create output file
                out_file = results_folder + 'scan' + str(settings['loop']) + '.csv'
                with open(out_file, 'w') as w:
                    w.write('Scan no,Time,Motor Position,Angle,so2 (ppm.m),so2 error,')
                            
                    for i in common['params'].keys():
                        w.write(i + ',' + i + '_e,')
                        
                    # Starrt new line 
                    w.write('\n')
                
                    # Analyse rest (106 spectra to a scan, first is dark, last is crap)
                    for n in range(1, len(spec_block.T)-1):
                        
                        # End loop if finished
                        if settings['stop_flag'] == True:
                            break
    
                        # Load spectrum
                        y = spec_block[:,n]
    
                        # Fit the spectrum
                        fit_p,err_dict,y_data,gas_T,fit_flag = fit_spec(common, y, grid)
                            
                        # Unpack fit results
                        fit_dict  = {}
                        for m, l in enumerate(common['params'].keys()):
                            fit_dict[l] = fit_p[m]
                            
                        # Unpack spec no, timestamp and motor position
                        spec_no = str(info_block[0][n])
                        timestamp = str(int(info_block[1][n])) + ':' + \
                                    str(int(info_block[2][n])) + ':' + \
                                    str(int(info_block[3][n]))
                        motor_pos = str(info_block[4][n])
                        view_ang  = str(float(motor_pos) * 0.06 - 102)
                            
                        # Add the data to the results file
                        w.write(spec_no + ','+timestamp+','+motor_pos+','+view_ang)
                        
                        # Print so2 amount and error in ppm.m for ease
                        w.write(',' + str(fit_dict['so2_amt']/2.463e15) + ',' + \
                                str(err_dict['so2_amt']/2.463e15))
                        
                        # Print fit results and error for each parameter          
                        for l in common['params'].keys():
                            w.write(',' + str(fit_dict[l]) + ',' + str(err_dict[l]))
                            
                        # Start a new line
                        w.write('\n')
                        
                        # Add values to array for plotting
                        spec_nos.append(float(view_ang))
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
                            common['solar_resid'] = np.add(common['solar_resid'], resid)
                            resid_count += 1

                        if bool(settings['update_params']) == True:
                            if fit_flag == False:
                                common['params'] = initial_params.copy()
                                self.print_output('Fitting for spectrum '+str(spec_no)+\
                                                  ' failed, resetting parameters')
                               
                            elif max_resid > 0.15:
                                common['params'] = initial_params.copy()
                                self.print_output('Fitting for spectrum '+str(spec_no)+\
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
                            
                            
                            # Get selected transmittance data
                            gas_tran = gas_T[settings['analysis_gas'] + '_tran']
                            gas_spec = gas_T[settings['analysis_gas'] + '_spec']
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
                            #                 x data    y data    x limits   y limits
                            data = np.array(([grid,     y_data,   'auto',  [y_lo,y_hi]],
                                             [grid,     fit,      'auto',  [y_lo,y_hi]],
                                             [x   ,     y,        'auto',  'auto'     ],
                                             [grid,     resid,    'auto',  'auto'     ],
                                             [grid,     gas_tran, 'auto',  [t_lo,t_hi]],
                                             [grid,     gas_spec, 'auto',  [t_lo,t_hi]],
                                             [spec_nos, gas_amts, 'auto',  'auto'     ]))
   
                            # Update graph
                            update_graph(lines, axes, self.canvas, data)
                        
                        # Update loop counter display
                        msg = str(n) + ' / ' + str(len(spec_block.T))
                        self.spec_count.set(msg)
                        
                        # Force gui to update
                        mygui.update(self)
                        
            else:
                self.print_output('Error in file ' + str(settings['loop']))  
           
            # Update loop counter display
            msg = str(settings['loop']) + ' / ' + str(len(common['scan_files']))
            self.scan_count.set(msg)
                     
            # Add to the count cycle
            settings['loop'] += 1   
                
        self.print_output('Analysis Complete!')        
                
        # Update solar residual
        if settings['solar_resid_flag'] == 'Generate':
            
            # Find average residual
            common['solar_resid'] = np.divide(common['solar_resid'], resid_count)
            
            # Save
            np.savetxt('data_bases/solar_resid.txt', common['solar_resid'])
            
            self.print_output('Solar residual spectrum updated')
















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
        settings['solar_resid_flag'] = resid_b.get()
        
        # Close the window
        popup.destroy()
    
    # Create row number counter
    row_n = 0
    
    # Set resolution of model grid
    popup.model_res = tk.DoubleVar(popup, value = settings['model_resolution'])
    model_res_l = tk.Label(popup, text = 'Model Grid spacing:', font = NORM_FONT)
    model_res_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    model_res_e = ttk.Entry(popup, textvariable = popup.model_res)
    model_res_e.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove dark spectra
    dark_b = tk.BooleanVar(popup, value = settings['dark_flag'])
    dark_l = tk.Label(popup, text = 'Remove Dark Spectra?', font = NORM_FONT)
    dark_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    dark_c = ttk.Checkbutton(popup, variable = dark_b)
    dark_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove flat spectra
    flat_b = tk.BooleanVar(popup, value = settings['flat_flag'])
    flat_l = tk.Label(popup, text = 'Remove Flat Spectra?', font = NORM_FONT)
    flat_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    flat_c = ttk.Checkbutton(popup, variable = flat_b)
    flat_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to remove, form or ignore the solar residual spectrum
    resid_options = [settings['solar_resid_flag'],
                     'Ignore',
                     'Generate',
                     'Remove']
    resid_b = tk.StringVar(popup, value = settings['solar_resid_flag'])
    resid_l = tk.Label(popup, text = 'Solar residual:', font = NORM_FONT)
    resid_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    resid_c = ttk.OptionMenu(popup, resid_b, *resid_options)
    resid_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
    row_n += 1
    
    # Control whether or not to update fit parameter guesses with the last fit values
    update_b = tk.BooleanVar(popup, value = settings['update_params'])
    update_l = tk.Label(popup, text = 'Auto-update fit parameters?', font = NORM_FONT)
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
def graph_settings(self):
    popup = tk.Tk()
    tk.Tk.wm_title(popup, 'Graph Settings')
    
    def update_graph_settings(settings):
        
        # Update graph settings in common
        settings['Show Graphs']     = graph_b.get()
        settings['Show Error Bars'] = err_b.get()
        settings['scroll_flag']     = scroll_b.get()
        settings['scroll_spec_no']  = int(spec_no_e.get())
        settings['resid_type']      = resid_type.get()
        settings['analysis_gas']    = gas.get()
        
        # Update graph
        if settings['resid_type'] == 'Percentage':
            self.ax2.set_ylabel('Fit residual (%)', fontsize=10)
        if settings['resid_type'] == 'Absolute':
            self.ax2.set_ylabel('Fit residual (Abs)', fontsize=10)
        if settings['resid_type'] == 'Spec/Fit':
            self.ax2.set_ylabel('Fit residual (Spec/Fit)', fontsize=10)
            
        self.ax3.set_ylabel(settings['analysis_gas'] + ' Transmittance', fontsize = 10)
        
        self.ax4.set_ylabel(settings['analysis_gas'] + ' amt (ppm.m)', fontsize = 10)
            
        self.canvas.draw()
        
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
    
    # Select which gas to analyse
    gas_options = [settings['analysis_gas'],
                   'SO2',
                   'O3',
                   'BrO']
    gas = tk.StringVar(popup, value = settings['scroll_flag'])
    gas_l = tk.Label(popup, text = 'Gas to analyse:', font = NORM_FONT)
    gas_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    gas_c = ttk.OptionMenu(popup, gas, *gas_options)
    gas_c.grid(row = row_n, column = 1, padx = 5, pady = 5)
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
    
    # Set format of residual
    resid_options = [settings['resid_type'],
                     'Absolute',
                     'Percentage',
                     'Spec/Fit']
    resid_type = tk.StringVar(popup, value = settings['resid_type'])
    resid_type_l = tk.Label(popup, text = 'Residual Display', font = NORM_FONT)
    resid_type_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    resid_type_m = ttk.OptionMenu(popup, resid_type, *resid_options)
    resid_type_m.grid(row = row_n, column = 1, padx = 5, pady = 5)
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
        settings['sol_path']         = sol_path_e.get()
        settings['ring_path']        = ring_path_e.get()
        settings['so2_path']         = so2_path_e.get()
        settings['no2_path']         = no2_path_e.get()
        settings['o3_path']          = o3_path_e.get()
        settings['bro_path']         = bro_path_e.get()
        settings['solar_resid_path'] = solar_resid_path_e.get()
        
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
    
    # Solar residual
    popup.solar_resid_path = tk.StringVar(popup, value = settings['solar_resid_path'])
    solar_resid_path_l = tk.Label(popup, text = 'Solar Residual:', font = NORM_FONT)
    solar_resid_path_l.grid(row = row_n, column = 0, padx = 5, pady = 5, sticky = 'W')
    solar_resid_path_e = ttk.Entry(popup, textvariable = popup.solar_resid_path,
                                   font = NORM_FONT, width = 40)
    solar_resid_path_e.grid(row = row_n, column = 1, padx = 5, pady = 5) 
    solar_resid_path_b = ttk.Button(popup, text = "Select", 
                                    command = lambda: update_fp(popup.solar_resid_path))
    solar_resid_path_b.grid(row = row_n, column = 2, padx = 5, pady = 5, sticky = 'W')
    row_n += 1
    
    # Button to apply changes and close
    b1 = ttk.Button(popup, text='Apply', command=lambda: update_data_bases(settings))
    b1.grid(row = row_n, column = 0, padx = 5, pady = 5, columnspan=2)
    b2 = ttk.Button(popup, text='Cancel', command=lambda: popup.destroy())
    b2.grid(row = row_n, column = 1, padx = 5, pady = 5, columnspan=2) 
    row_n += 1
          
# Tkinter stuff 
if __name__ == '__main__':    
    mygui().mainloop()
