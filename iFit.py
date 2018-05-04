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
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

from ifit_lib.build_fwd_data import build_fwd_data
from ifit_lib.read_spectrum import read_spectrum, average_spectra
from ifit_lib.fit import fit_spec, gen_fit_norm, ifit_fwd, gen_fit_ils, ifit_ils_fwd, ifit_ldf_fwd, gen_fit_ldf
from ifit_lib.find_nearest import extract_window
from ifit_lib.make_ils import make_ils

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
        
        # Add a title and icon
        tk.Tk.wm_title(self, 'iFit-2-1')
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
       
        # Create frames for diferent sections
        break_frame0 = tk.Frame(self, width = 10, pady=3)
        break_frame0.grid(row=0, column=0, rowspan=10, sticky="ew")
        
        break_frame1 = tk.Frame(self, height = 10, pady=3)
        break_frame1.grid(row=0, column=1, columnspan = 2, sticky="ew")        
        
        setup_frame = ttk.Frame(self, padding=3, borderwidth=5, relief = tk.GROOVE)
        setup_frame.grid(row=1, column=1, sticky="ew")
        
        break_frame2 = tk.Frame(self, height = 10, pady=3)
        break_frame2.grid(row=2, column=1, sticky="ew")
        
        param_frame = ttk.Frame(self, padding=3, borderwidth=5, relief = tk.GROOVE)
        param_frame.grid(row=3, column=1, sticky="ew")
        
        break_frame3 = tk.Frame(self, height = 10, pady=3)
        break_frame3.grid(row=4, column=1, sticky="ew")
        
        button_frame = ttk.Frame(self, padding=3, borderwidth=5, relief = tk.GROOVE)
        button_frame.grid(row=5, column=1, sticky="ew")
        
        break_frame4 = tk.Frame(self, height = 10, pady=3)
        break_frame4.grid(row=6, column=1, sticky="ew")
        break_frame6 = tk.Frame(self, width = 10, pady=3)
        break_frame6.grid(row=0, column=2,rowspan = 10, sticky="ew")
        
        graph_frame = ttk.Frame(self, padding=3, borderwidth=5, relief = tk.GROOVE)
        graph_frame.grid(row=1, column=3, rowspan = 10, sticky="NW")
        
        break_frame7 = tk.Frame(self, width = 10, pady=3)
        break_frame7.grid(row=0, column=4,rowspan = 10, columnspan = 10, sticky="ew")
        
        mygui.columnconfigure(index=3, weight=1, self = self)
        mygui.rowconfigure(index = 5, weight = 1, self = self)
        graph_frame.columnconfigure(index=0, weight=1)
        graph_frame.rowconfigure(index = 0, weight = 1)
        button_frame.rowconfigure(index=1, weight = 1)
        
#========================================================================================
#===================================Define ttk styles====================================
#========================================================================================
        
        # Button Style
        ttk.Style().configure('TButton', width = 20, height = 20, relief="flat")     
        
#========================================================================================
#====================================Create text output==================================
#======================================================================================== 
       
        # Create a scroll bar
        scrollbar = ttk.Scrollbar(button_frame)
        scrollbar.grid(row=1, column=3, sticky='NSE')
        
        # Build text box
        self.text_box = tk.Text(button_frame, width = 60, height = 10, 
                                yscrollcommand = scrollbar.set)
        self.text_box.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W',
                           columnspan = 3)
        self.text_box.insert('1.0', 'Welcome to iFit! Written by Ben Esse\n\n')
        
        # Connect the scrollbar to the textbox
        scrollbar.config(command = self.text_box.yview)
        
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
            settings['Wave Start']        = 310
            settings['Wave Stop']         = 319.9
            settings['Spectrometer']      = '-select-'
            settings['Spectra Type']      = '-select-'
            settings['ILS Width']         = 0.25
            settings['Fit ILS']           = 0
            settings['LDF']               = 0.0
            settings['Fit LDF']           = 0
            settings['Show Graphs']       = 1
            settings['Show Error Bars']   = 1
            settings['scroll_flag']       = 1
            settings['scroll_spec_no']    = 200
            settings['a']                 = 1.0
            settings['b']                 = 1.0
            settings['c']                 = 1.0
            settings['d']                 = 1.0
            settings['shift']             = -0.2
            settings['stretch']           = 0.05
            settings['ring']              = 1.0
            settings['SO2']               = 1e+16
            settings['NO2']               = 1e+17
            settings['O3']                = 1e+19
            settings['BrO']               = 1e+15
            settings['Spectra Filepaths'] = ''
            settings['Dark Filepaths']    = ''
                
        setup_title_l = tk.Label(setup_frame, text = 'Program Setup:', font = LARG_FONT)
        setup_title_l.grid(row = 0, column = 0, padx = 5, pady = 5, columnspan = 2, 
                           stick = "W")
        
        # Create entry for start and stop wavelengths
        self.wave_start = tk.IntVar(setup_frame, value = settings['Wave Start'])
        self.wave_start_l = tk.Label(setup_frame, text = 'Wave Start:', font = NORM_FONT)
        self.wave_start_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.wave_start_e = ttk.Entry(setup_frame, textvariable = self.wave_start)
        self.wave_start_e.grid(row = 1, column = 1, padx = 5, pady = 5)
        
        self.wave_stop = tk.IntVar(setup_frame, value = settings['Wave Stop'])
        self.wave_stop_l = tk.Label(setup_frame, text = 'Wave Stop:', font = NORM_FONT)
        self.wave_stop_l.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.wave_stop_e = ttk.Entry(setup_frame, textvariable = self.wave_stop)
        self.wave_stop_e.grid(row = 1, column = 3, padx = 5, pady = 5)
        
        # Create entry to select spectrometer
        options = [settings['Spectrometer'],
                   'FLMS02101',
                   'FLMS02929',
                   'USB2+H04173',
                   'USB2+H15972',
                   'USB2+F02057']
        self.spec_name = tk.StringVar(setup_frame, value = options[0])
        spectro_l = tk.Label(setup_frame, text = 'Spectrometer:', font = NORM_FONT)
        spectro_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        spectro_c = ttk.OptionMenu(setup_frame, self.spec_name, *options)
        spectro_c.grid(row = 2, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        # Create entry to select spectra type
        spec_options = [settings['Spectra Type'],
                        'IFRiT',
                        'Master.Scope',
                        'Jai Spec',
                        'Spectrasuite']
        self.spec_type = tk.StringVar(setup_frame, value = spec_options[0])
        self.spec_l = tk.Label(setup_frame, text = 'Spectra Type:', font = NORM_FONT)
        self.spec_l.grid(row = 2, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.spec_c = ttk.OptionMenu(setup_frame, self.spec_type, *spec_options)
        self.spec_c.grid(row = 2, column = 3, padx = 5, pady = 5, sticky = 'W')
        
        # Instrument lineshape width
        self.ils_width = tk.IntVar(setup_frame, value = settings['ILS Width'])
        self.ils_width_b = tk.IntVar(setup_frame, value = settings['Fit ILS'])
        self.ils_width_l = tk.Label(setup_frame, text = 'ILS Width:', font = NORM_FONT)
        self.ils_width_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.ils_width_e = ttk.Entry(setup_frame, textvariable = self.ils_width)
        self.ils_width_e.grid(row = 3, column = 1, padx = 5, pady = 5)
        self.ils_width_l2 = tk.Label(setup_frame, text = 'Fit ILS?', font = NORM_FONT)
        self.ils_width_l2.grid(row = 3, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.ils_width_c = ttk.Checkbutton(setup_frame, variable = self.ils_width_b)
        self.ils_width_c.grid(row = 3, column = 3, padx = 5, pady = 5)
        
        # Light dilution factor
        self.ldf = tk.IntVar(setup_frame, value = settings['LDF'])
        self.ldf_b = tk.IntVar(setup_frame, value = settings['Fit LDF'])
        self.ldf_l = tk.Label(setup_frame, text = 'LDF:', font = NORM_FONT)
        self.ldf_l.grid(row = 4, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.ldf_e = ttk.Entry(setup_frame, textvariable = self.ldf)
        self.ldf_e.grid(row = 4, column = 1, padx = 5, pady = 5)
        self.ldf_l2 = tk.Label(setup_frame, text = 'Fit LDF?', font = NORM_FONT)
        self.ldf_l2.grid(row = 4, column = 2, padx = 5, pady = 5, sticky = 'W')
        self.ldf_c = ttk.Checkbutton(setup_frame, variable = self.ldf_b)
        self.ldf_c.grid(row = 4, column = 3, padx = 5, pady = 5)
        
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
        self.spec_ent=tk.StringVar(value=str(len(self.spec_fpaths))+' spectra selected')
        self.specfp_l = tk.Entry(setup_frame, font = NORM_FONT, width = 40, 
                                 text = self.spec_ent)
        self.specfp_l.grid(row = 5, column = 1, padx = 5, pady = 5, sticky = 'W', 
                      columnspan = 3)
        specfp_b = ttk.Button(setup_frame, text="Select Spectra", command = self.spec_fp)
        specfp_b.grid(row = 5, column = 0, padx = 5, pady = 5, sticky = 'W')
        
        # File dialouge for darks
        self.dark_ent=tk.StringVar(value=str(len(self.dark_fpaths))+' spectra selected')
        self.darkfp_l = tk.Entry(setup_frame, font = NORM_FONT, width = 40, 
                                 text = self.dark_ent)
        self.darkfp_l.grid(row = 6, column = 1, padx = 5, pady = 5, sticky = 'W', 
                      columnspan = 3)
        darkfp_b = ttk.Button(setup_frame, text = "Select Darks", command = self.dark_fp)
        darkfp_b.grid(row = 6, column = 0, padx = 5, pady = 5, sticky = 'W')
               
#========================================================================================
#================================Set initial fit parameters==============================
#========================================================================================
        
        # Heading    
        params_l = tk.Label(param_frame, text = 'Fit parameters:', font = LARG_FONT)
        params_l.grid(row = 0, column = 0, padx = 5, pady = 5, columnspan = 2, 
                           stick = "W")
        
        # Polynomial coefficents
        self.a = tk.DoubleVar(self, value = settings['a'])
        a_l = tk.Label(param_frame, text = 'a:', font = NORM_FONT)
        a_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        a_e = ttk.Entry(param_frame, textvariable = self.a)
        a_e.grid(row = 1, column = 1, padx = 5, pady = 5)

        self.b = tk.DoubleVar(param_frame, value = settings['b'])
        b_l = tk.Label(param_frame, text = 'b:', font = NORM_FONT)
        b_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        b_e = ttk.Entry(param_frame, textvariable = self.b)
        b_e.grid(row = 2, column = 1, padx = 5, pady = 5)
        
        self.c = tk.DoubleVar(param_frame, value = settings['c'])
        c_l = tk.Label(param_frame, text = 'c:', font = NORM_FONT)
        c_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
        c_e = ttk.Entry(param_frame, textvariable = self.c)
        c_e.grid(row = 3, column = 1, padx = 5, pady = 5)
        
        self.d = tk.DoubleVar(param_frame, value = settings['d'])
        d_l = tk.Label(param_frame, text = 'd:', font = NORM_FONT)
        d_l.grid(row = 4, column = 0, padx = 5, pady = 5, sticky = 'W')
        d_e = ttk.Entry(param_frame, textvariable = self.d)
        d_e.grid(row = 4, column = 1, padx = 5, pady = 5)
        
        # Spectrometer wavelength shift parameters
        self.shift = tk.DoubleVar(param_frame, value = settings['shift'])
        shift_l = tk.Label(param_frame, text = 'shift:', font = NORM_FONT)
        shift_l.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
        shift_e = ttk.Entry(param_frame, textvariable = self.shift)
        shift_e.grid(row = 1, column = 3, padx = 5, pady = 5)
        
        self.stretch = tk.DoubleVar(param_frame, value = settings['stretch'])
        stretch_l = tk.Label(param_frame, text = 'stretch:', font = NORM_FONT)
        stretch_l.grid(row = 2, column = 2, padx = 5, pady = 5, sticky = 'W')
        stretch_e = ttk.Entry(param_frame, textvariable = self.stretch)
        stretch_e.grid(row = 2, column = 3, padx = 5, pady = 5)
        
        # Ring effect
        self.ring_amt = tk.DoubleVar(param_frame, value = settings['ring'])
        ring_amt_l = tk.Label(param_frame, text = 'Ring:', font = NORM_FONT)
        ring_amt_l.grid(row = 3, column = 2, padx = 5, pady = 5, sticky = 'W')
        ring_amt_e = ttk.Entry(param_frame, textvariable = self.ring_amt)
        ring_amt_e.grid(row = 3, column = 3, padx = 5, pady = 5)
        
        # Gas amounts
        self.so2_amt = tk.DoubleVar(param_frame, value = settings['SO2'])
        so2_amt_l = tk.Label(param_frame, text = 'SO2:', font = NORM_FONT)
        so2_amt_l.grid(row = 1, column = 4, padx = 5, pady = 5, sticky = 'W')
        so2_amt_e = ttk.Entry(param_frame, textvariable = self.so2_amt)
        so2_amt_e.grid(row = 1, column = 5, padx = 5, pady = 5)
        
        self.no2_amt = tk.DoubleVar(param_frame, value = settings['NO2'])
        no2_amt_l = tk.Label(param_frame, text = 'NO2:', font = NORM_FONT)
        no2_amt_l.grid(row = 2, column = 4, padx = 5, pady = 5, sticky = 'W')
        no2_amt_e = ttk.Entry(param_frame, textvariable = self.no2_amt)
        no2_amt_e.grid(row = 2, column = 5, padx = 5, pady = 5)
        
        self.o3_amt = tk.DoubleVar(param_frame, value = settings['O3'])
        o3_amt_l = tk.Label(param_frame, text = 'O3:', font = NORM_FONT)
        o3_amt_l.grid(row = 3, column = 4, padx = 5, pady = 5, sticky = 'W')
        o3_amt_e = ttk.Entry(param_frame, textvariable = self.o3_amt)
        o3_amt_e.grid(row = 3, column = 5, padx = 5, pady = 5)
       
        self.bro_amt = tk.DoubleVar(param_frame, value = settings['BrO'])
        bro_amt_l = tk.Label(param_frame, text = 'BrO:', font = NORM_FONT)
        bro_amt_l.grid(row = 4, column = 4, padx = 5, pady = 5, sticky = 'W')
        bro_amt_e = ttk.Entry(param_frame, textvariable = self.bro_amt)
        bro_amt_e.grid(row = 4, column = 5, padx = 5, pady = 5)
          
#========================================================================================
#==============================Create start and exit buttons=============================
#========================================================================================         
        
        # Create button to start
        start_b = ttk.Button(button_frame, text = 'Begin!', command = self.begin)
        start_b.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        
        # Create button to save settings
        save_b = ttk.Button(button_frame, text = 'Save Settings', command = self.save)
        save_b.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        # Create a button to exit
        exit_b = ttk.Button(button_frame, text = 'Exit', command = self.quit)
        exit_b.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = 'W')
        
#========================================================================================
#====================================Create plot canvas==================================
#========================================================================================        
                    
        # Create figure to hold the graphs
        plt.rcParams.update({'font.size': 8} )
        self.fig = plt.figure(figsize = (6,6))
        gs = gridspec.GridSpec(3,1, height_ratios = (3,2,3))
        
        # Create plot axes
        self.ax0 = self.fig.add_subplot(gs[0])
        self.ax1 = self.fig.add_subplot(gs[1])
        self.ax2 = self.fig.add_subplot(gs[2])
        
        # Five axes: 1) Spectrum and fit
        #            2) Residual
        #            3) SO2 amount time series
        
        # Create axis title
        self.ax0.set_title('Spectrum 0 / 0', fontsize=10)
        
        # Set axis labels
        self.ax0.set_ylabel('Intensity (arb)', fontsize=10)
        
        self.ax1.set_ylabel('Fit residual (%)', fontsize=10)
        self.ax1.set_xlabel('Wavelength (nm)', fontsize=10)
        
        self.ax2.set_ylabel('SO2 amt (ppm.m)', fontsize=10)
        self.ax2.set_xlabel('Spectrum number', fontsize=10)
        
        # Create lines to plot data series
        
        # Spectral data
        self.line1, = self.ax0.plot(0, 0, 'b', label = 'Spectrum')
        self.line2, = self.ax0.plot(0, 0, 'darkorange', label = 'Fit')
        self.ax0.legend(loc = 0)
        
        # Fit residual
        self.line3, = self.ax1.plot(0, 0, 'r')
        
        # Full spectrum
        self.line4, = self.ax2.plot(0, 0, 'g')
        
        # Pack into arrays to pass to plotting function
        self.lines = [self.line1, self.line2, self.line3, self.line4]
        self.axes  = [self.ax0,   self.ax0,   self.ax1,   self.ax2  ]
        
        # Make it look nice
        plt.tight_layout()
        
        # Create the canvas to hold the graph in the GUI
        self.canvas = FigureCanvasTkAgg(self.fig, graph_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=0, column=0, padx=10)
        
        # Add matplotlib toolbar above the plot canvas
        toolbar_frame = tk.Frame(graph_frame, bg = 'black')  
        toolbar_frame.grid(row=1,column=0, sticky = 'W')                             
        toolbar = NavigationToolbar2TkAgg(self.canvas, toolbar_frame)
        toolbar.update()     
        
#========================================================================================         
#========================================================================================
#=====================================Button Functions===================================
#======================================================================================== 
#======================================================================================== 

#========================================================================================
#=====================================Filepath Buttons===================================
#========================================================================================
    
    def spec_fp(self):
            
        # Open dialouge to get files
        fpaths = fd.askopenfilenames()

        if fpaths != None:
            self.spec_fpaths = []
            for i in fpaths:
                self.spec_fpaths.append(str(i))
            
            # Save output to input 
            self.spec_ent.set(str(len(self.spec_fpaths)) + ' spectra selected')
    
    def dark_fp(self):
            
        # Open dialouge to get files
        fpaths = fd.askopenfilenames()
        
        if fpaths != None:
            self.dark_fpaths = []
            for i in fpaths:
                self.dark_fpaths.append(str(i))
            
            # Save output to input
            self.dark_ent.set(str(len(self.dark_fpaths)) + ' spectra selected')
             
#========================================================================================
#=======================================Save Settings====================================
#========================================================================================
            
    def save(self):
        
        # Create or overright settings file
        with open('data_bases/ifit_settings.txt', 'w') as w:
            
            # Save each setting into the file
            w.write('Wave Start;'       + str(self.wave_start_e.get())      + '\n')
            w.write('Wave Stop;'        + str(self.wave_stop_e.get())       + '\n')
            w.write('Spectrometer;'     + str(self.spec_name.get())         + '\n')
            w.write('Spectra Type;'     + str(self.spec_type.get())         + '\n')
            w.write('ILS Width;'        + str(self.ils_width_e.get())       + '\n')
            w.write('Fit ILS;'          + str(self.ils_width_b.get())       + '\n')
            w.write('LDF;'              + str(self.ldf_e.get())             + '\n')
            w.write('Fit LDF;'          + str(self.ldf_b.get())             + '\n')
            w.write('Show Graphs;'      + str(settings['Show Graphs'])      + '\n')
            w.write('Show Error Bars;'  + str(settings['Show Error Bars'])  + '\n')
            w.write('scroll_flag;'      + str(settings['scroll_flag'])      + '\n')
            w.write('scroll_spec_no;'   + str(settings['scroll_spec_no'])   + '\n')
            w.write('a;'                + str(self.a.get())                 + '\n')
            w.write('b;'                + str(self.b.get())                 + '\n')
            w.write('c;'                + str(self.c.get())                 + '\n')
            w.write('d;'                + str(self.d.get())                 + '\n')
            w.write('shift;'            + str(self.shift.get())             + '\n')
            w.write('stretch;'          + str(self.stretch.get())           + '\n')
            w.write('ring;'             + str(self.ring_amt.get())          + '\n')
            w.write('SO2;'              + str(self.so2_amt.get())           + '\n')
            w.write('NO2;'              + str(self.no2_amt.get())           + '\n')
            w.write('O3;'               + str(self.o3_amt.get())            + '\n')
            w.write('BrO;'              + str(self.bro_amt.get())           + '\n')
            w.write('model_resolution;' + str(settings['model_resolution']) + '\n')
            w.write('sol_path;'         + str(settings['sol_path'])         + '\n')
            w.write('resid_path;'       + str(settings['resid_path'])       + '\n')
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
#========================================Begin iFit======================================
#========================================================================================
            
    def begin(self):
        
        # Create common dictionary
        common = {}
        
        # Build parameter array
        common['params'] = [['a',        float(self.a.get())       ], 
                            ['b',        float(self.b.get())       ],
                            ['c',        float(self.c.get())       ],
                            ['d',        float(self.d.get())       ],
                            ['shift',    float(self.shift.get())   ], 
                            ['stretch',  float(self.stretch.get()) ], 
                            ['ring_amt', float(self.ring_amt.get())],
                            ['so2_amt',  float(self.so2_amt.get()) ],
                            ['no2_amt',  float(self.no2_amt.get()) ],       
                            ['o3_amt',   float(self.o3_amt.get())  ]] 

        # Populate common with other data from the GUI
        common['wave_start']       = float(self.wave_start_e.get())
        common['wave_stop']        = float(self.wave_stop_e.get())
        common['ils_width']        = float(self.ils_width_e.get())
        common['ils_gauss_weight'] = 1.0
        common['ldf']              = float(self.ldf_e.get())
        common['solar_resid_flag'] = 'Ignore'
        common['ils_flag']         = self.ils_width_b.get()
        common['ldf_flag']         = self.ldf_b.get()
        common['spectra_files']    = self.spec_fpaths
        common['dark_files']       = self.dark_fpaths
        common['spec_name']        = self.spec_name.get()
        
        # Read format of spectra
        spec_type = self.spec_type.get()

        # Filepaths to solar, ring, flat and gas xsec spectra
        settings['flat_path']  = 'data_bases/flat_' + str(self.spec_name.get()) + '.txt'
        
        # Load fitting data files
        common = build_fwd_data(self, common, settings)
        
#========================================================================================
#====================Adjust initial parameters if fitting extra params===================
#========================================================================================

        # If fitting ILS, select correct forward model and add ils to params
        if common['ils_flag'] == True:

            fwd_model = ifit_ils_fwd
            gen_fit   = gen_fit_ils
            
            common['params'].append(['ils_width', common['ils_width']])
            
            # Create empty ils parameter
            common['ils'] = None
            
        # If fitting LDF, select correct forward model and add ils to params
        elif common['ldf_flag'] == True:

            fwd_model = ifit_ldf_fwd
            gen_fit   = gen_fit_ldf
            
            common['params'].append(['ldf', common['ldf']])
            
            # Generate spectrometer ILS, a smoothing function applied to the fwdmodel
            common['ils'] = make_ils(common['model_grid'], common['ils_width'], 
                                     common['ils_gauss_weight'])
            
        else:
            fwd_model = ifit_fwd
            gen_fit = gen_fit_norm
            
            # Generate spectrometer ILS, a smoothing function applied to the fwdmodel
            common['ils'] = make_ils(common['model_grid'], common['ils_width'], 
                                     common['ils_gauss_weight'])

        # Save initial guess parameters
        common['initial_params'] = common['params']

#========================================================================================
#===================================Build dark spectrum==================================
#========================================================================================

        # Read dark spectrum
        x, common['dark'] = average_spectra(common['dark_files'], spec_type)
        
        # Find indices of desired wavelength window and add to common
        grid, common['ind1'], common['ind2'] = extract_window(x, common['wave_start'], 
                                                              common['wave_stop'])
        
        # Find stray light window
        stray_grid, common['stray_i1'], common['stray_i2'] = extract_window(x, 282, 290)
        
        # Create flag to control whether stray light is removed
        common['stray_flag'] = True
        
        # If no stray light pixels available, turn off the flag
        if common['stray_i1'] == common['stray_i2']:
            common['stray_flag'] = False 

#========================================================================================
#====================================Open output files===================================
#========================================================================================
    
        # Read first dark spectrum to get date of data
        x, y, data_date, data_time, spec_no = read_spectrum(common['dark_files'][0],
                                                            spec_type)
        
        # Create directory to hold program outputs
        results_folder = 'Results/iFit/' + data_date + '/'
        
        # Create folder if it doesn't exist
        if not os.path.exists(results_folder):
                # Make the directory
                os.makedirs(results_folder)
        
        # Create filename for output file
        out_excel_fname = results_folder + 'iFit_out.csv'
        
        try:
            
            # Open excel file and write header line
            with open(out_excel_fname, 'w') as writer:
                
                # Write header line
                writer.write('File,Number,Date,Time,so2 (ppm.m),so2 error,a,a_e,b,b_e,c,'+\
                             'c_e,d,d_e,shift,shift_e,stretch,stretch_e,ring,ring_e,so2,'+\
                             'so2_e,no2,no2_e,o3,o3_e,')
                
                if common['ils_flag'] == True:
                    writer.write('ILS width,ILS width_e,')  
                    
                if common['ldf_flag'] == True:
                    writer.write('LDF,LDF_e,')
                    
                # Write other fit info
                writer.write('Fit window: '+str(common['wave_start'])+' - '+str(common['wave_stop'])+\
                             ' nm,' + 'ILS width: ' + str(common['ils_width']) + \
                             'ILS Gauss Weight: ' + str(common['ils_gauss_weight']) + '\n')
                
        except PermissionError:
            self.print_output('Please close iFit output file to continue')
            return

#========================================================================================
#===================================Start Analysis Loop==================================
#========================================================================================

        # Open excel file and write header line
        with open(out_excel_fname, 'a') as writer:

            self.print_output('Begin fitting')
        
            count = 0
            
            # Create empty arrays to hold the loop number and so2_amt values
            spec_nos    = []
            so2_amts    = []
            amt_errs    = []
            solar_resid = np.zeros(len(grid))

            for fname in common['spectra_files']:
                
                # Read in spectrum file
                try:
                    x, y, data_date, data_time, spec_no = read_spectrum(fname,
                                                                        spec_type)
                    
                except FileNotFoundError:
                    self.print_output('File number ' + str(count) + ' not found')
                    
                except OSError:
                    self.print_output('File number ' + str(count) + ' not found')
                    
                else:

#========================================================================================
#=========================Analyse spectrum and add to output file========================
#========================================================================================
                    
                    # Fit
                    results = fit_spec(common, y, grid, fwd_model)
                    
                    fit_params, err_dict, y_data, fit_flag = results
                                                                             
                    # Unpack fit results
                    fit_dict  = {}
                    for m, l in enumerate(common['params']):
                        fit_dict[l[0]] = fit_params[m]
                                
                    # Write results to excel file, starting with spectrum info
                    writer.write(str(fname)     + ',' + \
                                 str(spec_no)   + ',' + \
                                 str(data_date) + ',' + \
                                 str(data_time))
                    
                    # Print so2 amount and error in ppm.m for ease of access
                    writer.write(',' + str(fit_dict['so2_amt']/2.463e15) + ',' + \
                                 str(err_dict['so2_amt']/2.463e15))            
        
                    # Print fit results and error for each parameter          
                    for l in common['params']:
                        writer.write(',' + str(fit_dict[l[0]]) + ',' + str(err_dict[l[0]]))
                        
                    # Start new line
                    writer.write('\n')
                
                    # Add values to array for last 200 spectra counts and so2_amt for plotting
                    spec_nos.append(spec_no)
                    so2_amts.append(fit_dict['so2_amt']/2.463e15)
                    amt_errs.append(err_dict['so2_amt']/2.463e15)
                       
                    # Cut if too long
                    if settings['scroll_flag'] == True:
                        if len(spec_nos) > settings['scroll_spec_no']:
                            spec_nos = spec_nos[1:]
                            so2_amts = so2_amts[1:]
                            amt_errs = amt_errs[1:]
            
#========================================================================================
#=======================================Update plot======================================
#========================================================================================
        
                    # Feed fit params into forward
                    fit = gen_fit(grid, fit_params)
                    
                    # Calculate the residual of the fit
                    resid = np.multiply(np.divide(np.subtract(y_data, fit), y_data), 100)
                    
                    # Add to solar residual
                    if common['solar_resid_flag'] == 'Save':
                        solar_resid = np.add(solar_resid, resid)

                    # Replot data
                    if int(settings['Show Graphs']) == 1:            
                        
                        # Set graph title as spectrum number
                        self.ax0.set_title('Spectrum ' + str(count) + ' / ' + \
                                           str(len(common['spectra_files']) - 1))
                        
                        # Graph 1: Spectrum and fit
                        y_low  = min(y_data) - abs((0.1*max(y_data)))
                        y_high = max(y_data) + abs((0.1*max(y_data)))
                        f_low  = min(fit) - abs((0.1*max(fit)))
                        f_high = max(fit) + abs((0.1*max(fit)))
                        low, high = min([f_low, y_low]), max([f_high, y_high])
                        self.line1.set_data(grid, y_data)
                        self.line2.set_data(grid, fit)
                        self.ax0.set_xlim(grid.min() - 1, grid.max() + 1)
                        self.ax0.set_ylim(low, high)
                        
                        # Graph 2: Residual
                        r_low  = min(resid) - abs((0.1*max(resid)))
                        r_high = max(resid) + abs((0.1*max(resid)))
                        self.line3.set_data(grid, resid)
                        self.ax1.set_xlim(grid.min() - 1, grid.max() + 1)
                        if max(np.abs(resid)) < 5:
                            self.ax1.set_ylim(-5, 5)
                        else:
                            self.ax1.set_ylim(r_low, r_high)
                        
                        # Graph 3: SO2 time series
                        s_low  = min(so2_amts) - abs((0.1*max(so2_amts)))
                        s_high = max(so2_amts) + abs((0.1*max(so2_amts)))
                        self.line4.set_data(spec_nos, so2_amts)
                        self.ax2.set_xlim(min(spec_nos) - 1, max(spec_nos) + 1)
                        self.ax2.set_ylim(s_low, s_high)

                        # Add error bars to SO2 output
                        if int(settings['Show Error Bars']) == True:
                            error = [np.subtract(so2_amts, amt_errs), 
                                     np.add(so2_amts, amt_errs)]
                            self.ax2.fill_between(spec_nos, error[0], error[1], 
                                                  color = 'lightgreen')
                        
                        self.canvas.draw()
                        
#========================================================================================
#==================================Update fit parameters=================================
#========================================================================================
                    
                    # If fitting fails or max resid > 10% revert to initial fit parameters
                    if fit_flag == False:
                        common['params'] = common['initial_params']
                        self.print_output('Fitting for spectrum ' + str(count) + \
                                          ' failed, resetting parameters')
                        
                    elif max((resid)**2)**0.5 > 10:
                        common['params'] = common['initial_params']
                        self.print_output('Fitting for spectrum ' + str(count) + \
                                          ' bad, resetting parameters')
                    
                    else:
                        
                        # Update parameters with those from the last fit
                        common['params'] = [['a',         fit_dict['a']         ],
                                            ['b',         fit_dict['b']         ],
                                            ['c',         fit_dict['c']         ],
                                            ['d',         fit_dict['d']         ],
                                            ['shift',     fit_dict['shift']     ],
                                            ['stretch',   fit_dict['stretch']   ],
                                            ['ring_amt',  fit_dict['ring_amt']  ],
                                            ['so2_amt',   fit_dict['so2_amt']   ],
                                            ['no2_amt',   fit_dict['no2_amt']   ],
                                            ['o3_amt',    fit_dict['o3_amt']    ]]
                        
                        if common['ils_flag'] == True:
                            common['params'].append(['ils_width', common['ils_width']])
                            
                        if common['ldf_flag'] == True:
                            common['params'].append(['ldf', common['ldf']])
                             
                    # Add to the count cycle
                    count += 1
                    
                    # Force gui to update
                    mygui.update(self)
                    
            self.print_output('Fitting complete')

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
        
        # Close the window
        popup.destroy()
        
    # Set resolution of model grid
    popup.model_res = tk.IntVar(popup, value = settings['model_resolution'])
    model_res_l = tk.Label(popup, text = 'Model Grid spacing:', font = NORM_FONT)
    model_res_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
    model_res_e = ttk.Entry(popup, textvariable = popup.model_res)
    model_res_e.grid(row = 0, column = 1, padx = 5, pady = 5)
    
    # Button to apply changes and close
    b1 = ttk.Button(popup, text='Apply', command=lambda: update_model_settings(settings))
    b1.grid(row = 1, column = 0, padx = 5, pady = 5)
    
    
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
        
    # Control graph display settings
    graph_b = tk.IntVar(popup, value = settings['Show Graphs'])
    graph_l = tk.Label(popup, text = 'Show Graphs?', font = NORM_FONT)
    graph_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
    graph_c = ttk.Checkbutton(popup, variable = graph_b)
    graph_c.grid(row = 0, column = 1, padx = 5, pady = 5)
    
    # Turn on/off error bars
    err_b = tk.IntVar(popup, value = settings['Show Error Bars'])
    err_l = tk.Label(popup, text = 'Show Error Bars?', font = NORM_FONT)
    err_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
    err_c = ttk.Checkbutton(popup, variable = err_b)
    err_c.grid(row = 1, column = 1, padx = 5, pady = 5)
    
    # Turn on/off graph scrolling
    scroll_b = tk.IntVar(popup, value = settings['scroll_flag'])
    scroll_l = tk.Label(popup, text = 'Scroll Graph?', font = NORM_FONT)
    scroll_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
    scroll_c = ttk.Checkbutton(popup, variable = scroll_b)
    scroll_c.grid(row = 2, column = 1, padx = 5, pady = 5)
    
    # Set number of spectra to display on graph
    popup.spec_no = tk.IntVar(popup, value = settings['scroll_spec_no'])
    spec_no_l = tk.Label(popup, text = 'No. Spectra to display', font = NORM_FONT)
    spec_no_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
    spec_no_e = ttk.Entry(popup, textvariable = popup.spec_no, font = NORM_FONT)
    spec_no_e.grid(row = 3, column = 1, padx = 5, pady = 5)
    
    # Button to apply changes and close
    b1 = ttk.Button(popup, text='Apply', command=lambda: update_graph_settings(settings))
    b1.grid(row = 4, column = 0, padx = 5, pady = 5)


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
        settings['resid_path'] = resid_path_e.get()
        settings['ring_path']  = ring_path_e.get()
        settings['so2_path']   = so2_path_e.get()
        settings['no2_path']   = no2_path_e.get()
        settings['o3_path']    = o3_path_e.get()
        settings['bro_path']   = bro_path_e.get()
        
        # Close the window
        popup.destroy()
    
    # Create inputs for data base filepaths

    # Solar spectrum
    popup.sol_path = tk.StringVar(popup, value = settings['sol_path'])
    sol_path_l = tk.Label(popup, text = 'Solar spectrum:', font = NORM_FONT)
    sol_path_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
    sol_path_e = ttk.Entry(popup, text = popup.sol_path, font = NORM_FONT, width = 40)
    sol_path_e.grid(row = 0, column = 1, padx = 5, pady = 5)
    sol_path_b = ttk.Button(popup, text = "Select", command = lambda: update_fp(popup.sol_path))
    sol_path_b.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = 'W')
    
    # Solar residual
    popup.resid_path = tk.StringVar(popup, value = settings['resid_path'])
    resid_path_l = tk.Label(popup, text = 'Solar residual:', font = NORM_FONT)
    resid_path_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
    resid_path_e = ttk.Entry(popup, textvariable = popup.resid_path, font = NORM_FONT,
                             width = 40)
    resid_path_e.grid(row = 1, column = 1, padx = 5, pady = 5)
    resid_path_b = ttk.Button(popup, text = "Select", command = lambda: update_fp(popup.resid_path))
    resid_path_b.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
    
    # Ring spectrum
    popup.ring_path = tk.StringVar(popup, value = settings['ring_path'])
    ring_path_l = tk.Label(popup, text = 'Ring Spectrum:', font = NORM_FONT)
    ring_path_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
    ring_path_e = ttk.Entry(popup, textvariable = popup.ring_path, font = NORM_FONT, 
                            width = 40)
    ring_path_e.grid(row = 2, column = 1, padx = 5, pady = 5)
    ring_path_b = ttk.Button(popup, text = "Select", command = lambda: update_fp(popup.ring_path))
    ring_path_b.grid(row = 2, column = 2, padx = 5, pady = 5, sticky = 'W')
    
    # SO2 xsec
    popup.so2_path = tk.StringVar(popup, value = settings['so2_path'])
    so2_path_l = tk.Label(popup, text = 'SO2 xsec:', font = NORM_FONT)
    so2_path_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
    so2_path_e = ttk.Entry(popup, textvariable = popup.so2_path, font = NORM_FONT,
                           width = 40)
    so2_path_e.grid(row = 3, column = 1, padx = 5, pady = 5)
    so2_path_b = ttk.Button(popup, text = "Select", command = lambda: update_fp(popup.so2_path))
    so2_path_b.grid(row = 3, column = 2, padx = 5, pady = 5, sticky = 'W')
    
    # NO2 xsec
    popup.no2_path = tk.StringVar(popup, value = settings['no2_path'])
    no2_path_l = tk.Label(popup, text = 'NO2 xsec:', font = NORM_FONT)
    no2_path_l.grid(row = 4, column = 0, padx = 5, pady = 5, sticky = 'W')
    no2_path_e = ttk.Entry(popup, textvariable = popup.no2_path, font = NORM_FONT,
                           width = 40)
    no2_path_e.grid(row = 4, column = 1, padx = 5, pady = 5)
    no2_path_b = ttk.Button(popup, text = "Select", command = lambda: update_fp(popup.no2_path))
    no2_path_b.grid(row = 4, column = 2, padx = 5, pady = 5, sticky = 'W')
    
    # O3 xsec
    popup.o3_path = tk.StringVar(popup, value = settings['o3_path'])
    o3_path_l = tk.Label(popup, text = 'O3 xsec:', font = NORM_FONT)
    o3_path_l.grid(row = 5, column = 0, padx = 5, pady = 5, sticky = 'W')
    o3_path_e = ttk.Entry(popup, textvariable = popup.o3_path, font = NORM_FONT,
                          width = 40)
    o3_path_e.grid(row = 5, column = 1, padx = 5, pady = 5)
    o3_path_b = ttk.Button(popup, text = "Select", command = lambda: update_fp(popup.o3_path))
    o3_path_b.grid(row = 5, column = 2, padx = 5, pady = 5, sticky = 'W')
    
    # Bro xsec
    popup.bro_path = tk.StringVar(popup, value = settings['bro_path'])
    bro_path_l = tk.Label(popup, text = 'BrO xsec:', font = NORM_FONT)
    bro_path_l.grid(row = 6, column = 0, padx = 5, pady = 5, sticky = 'W')
    bro_path_e = ttk.Entry(popup, textvariable = popup.bro_path, font = NORM_FONT,
                           width = 40)
    bro_path_e.grid(row = 6, column = 1, padx = 5, pady = 5) 
    bro_path_b = ttk.Button(popup, text = "Select", command = lambda: update_fp(popup.bro_path))
    bro_path_b.grid(row = 6, column = 2, padx = 5, pady = 5, sticky = 'W')
    
    # Button to apply changes and close
    b1 = ttk.Button(popup, text='Apply', command=lambda: update_data_bases(settings))
    b1.grid(row = 7, column = 0, padx = 5, pady = 5)
            
# Tkinter stuff       
app = mygui()
app.mainloop()