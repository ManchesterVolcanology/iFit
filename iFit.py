# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 09:24:05 2018

@author: mqbpwbe2
"""

# Import required libraries
import matplotlib
matplotlib.use('TkAgg')
import traceback
import tkinter.messagebox as tkMessageBox
import tkinter.scrolledtext as tkst
import glob
import numpy as np
from tkinter import ttk
import tkinter as tk
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import OrderedDict
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

from ifit_lib.build_fwd_data import build_fwd_data
from ifit_lib.build_gui import make_input
from ifit_lib.fit import fit_spec
from ifit_lib.read_spectrum import read_spectrum
from ifit_lib.julian_time import hms_to_julian
from ifit_lib.update_graph import update_graph
from ifit_lib.control_loop import rt_setup, post_setup, rt_analyse
from ifit_lib.gui_funcs import adv_settings, fit_toggle, spec_fp, dark_fp, stop, \
                               connect_spec, update_int_time, test_spec, read_darks, \
                               read_settings

# Define some fonts to use in the program
NORM_FONT = ('Verdana', 8)
MED_FONT  = ('Veranda', 11)
LARG_FONT = ('Verdana', 12, 'bold')

class mygui(tk.Tk):
    
    def __init__(self, *args, **kwargs):
        
#========================================================================================
#================================ Build GUI Containers ==================================
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
        tk.Tk.wm_title(self, 'iFit-2-4')
        tk.Tk.iconbitmap(self, default = 'data_bases/icon.ico')
        
        # Create notebook to hold different frames
        self.nb = ttk.Notebook(self)
        page1 = ttk.Frame(self.nb)
        page2 = ttk.Frame(self.nb)
        
        # Create two frames, one for post analysis and one for real time acquisition
        self.nb.add(page2, text = 'Real Time Acquisition')
        self.nb.add(page1, text = 'Post Analysis')
        
        self.nb.grid(column=0, padx=10, pady=10, sticky = 'NW')
        
        # Create frame to hold graphs
        graph_frame = ttk.Frame(self, relief = 'groove')
        graph_frame.grid(row=0, column=1, padx=10, pady=10, rowspan=10, sticky="NW")
        graph_frame.columnconfigure(index=0, weight=1)
        graph_frame.rowconfigure(index = 0, weight = 1)
        
        # Frame for quick analysis
        quick_frame = tk.LabelFrame(self, text = 'Control', font = LARG_FONT)
        quick_frame.grid(row=1, column=0, padx=10, pady=10, sticky="NW")
        
        mygui.columnconfigure(index=1, weight=1, self = self)
        mygui.rowconfigure(index = 5, weight = 1, self = self)
        
#========================================================================================
#=================================== Create text output =================================
#========================================================================================  
        
        # Create frame to hold text output
        text_frame = ttk.Frame(quick_frame)
        text_frame.grid(row=3, column=0, padx=10, pady=10, columnspan=5, sticky="NW")      
                 
        # Build text box
        self.text_box = tkst.ScrolledText(text_frame, width = 42, height = 8)
        self.text_box.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W',
                           columnspan = 2)
        self.text_box.insert('1.0', 'Welcome to iFit! Written by Ben Esse\n\n')  
        
#========================================================================================
#================================== Set program settings ================================
#========================================================================================

        # Create settings dictionary
        global settings
        settings = {}
        
        # Read in settings file
        try:
            settings = read_settings('data_bases/ifit_settings.txt', settings)
   
        except FileNotFoundError:
            self.print_output('No settings file found\nReverting to origional')
            settings['wave_start']        = 305
            settings['wave_stop']         = 318
            settings['Spectrometer']      = '-select-'
            settings['Spectra Type']      = '-select-'
            settings['int_time']          = 100
            settings['coadds']            = 10
            settings['no_darks']          = 10
            settings['ils_width']         = 0.52
            settings['gauss_weight']      = 1.0
            settings['Fit ILS']           = 'File'
            settings['ldf']               = 0.0
            settings['Fit LDF']           = 'N/A'
            settings['dark_flag']         = True
            settings['flat_flag']         = True
            settings['update_params']     = True
            settings['good_fit_bound']    = 10
            settings['fit_weight']        = 'None'
            settings['Show Graphs']       = True
            settings['analysis_gas']      = 'so2'
            settings['scroll_flag']       = True
            settings['scroll_spec_no']    = 200
            settings['x_plot']            = 'Time'
            settings['resid_type']        = 'Spec/Fit'
            settings['solar_resid_flag']  = 'Ignore'
            settings['poly_n']            = 3
            settings['shift']             = -0.2
            settings['stretch']           = 0.05
            settings['ring_amt']          = 1.0
            settings['so2_amt']           = 1e+16
            settings['no2_amt']           = 1e+17
            settings['o3_amt']            = 1e+19
            settings['bro_amt']           = 1e+15
            settings['Fit shift']         = 'Fit'
            settings['Fit stretch']       = 'Fit'
            settings['Fit ring']          = 'Fit'
            settings['Fit so2']           = 'Fit'
            settings['Fit no2']           = 'Fit'
            settings['Fit o3']            = 'Fit'
            settings['Fit bro']           = 'Fit'
            settings['model_res']         = 0.01
            settings['model_pad']         = 3.0
            settings['sol_path']          = 'data_bases/gas data/sao2010.txt'
            settings['ring_path']         = 'data_bases/gas data/qdoas_ring.dat'
            settings['so2_path']          = 'data_bases/gas data/SO2_293K.dat'
            settings['no2_path']          = 'data_bases/gas data/No2_223l.dat'
            settings['o3_path']           = 'data_bases/gas data/O3_xsec.dat'
            settings['o3_temp']           = '233K'
            settings['bro_path']          = 'data_bases/gas data/BrO_Cross_298K.txt'
            settings['solar_resid_path']  = 'data_bases/gas data/solar_resid.txt'
            settings['Notebook page']     = 0
            settings['Spectra Filepaths'] = ''
            settings['Dark Filepaths']    = ''
 
        # Set notebook tab, and update on change
        def on_nb_change(event):
            settings['Notebook page'] = self.nb.index(self.nb.select())
            
        self.nb.select(int(settings['Notebook page']))
        self.nb.bind("<<NotebookTabChanged>>", func = on_nb_change)
        
        # Create loop counter to keep track of the analysis
        self.loop = 0
        
        # Create flag to ensure only one output file is created per program launch
        self.create_out_flag = True
        
        # Create flag to see if darks have been measured yet
        self.rt_dark_flag = False
        
        # Create flag to control whether or not to build the forward model
        self.build_model_flag = True
        self.common = {}
        
#========================================================================================
#=================================== Create plot canvas =================================
#========================================================================================        
            
        # Translate gas choice to parameter names and graph print
        gas_choice = {'so2' : r'SO$_2$',
                      'no2' : r'NO$_2$',
                      'o3'  : r'O$_3$' ,
                      'bro' : 'BrO'    ,
                      'ring': 'Ring'   }
        
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
        
        self.ax3.set_ylabel(gas_choice[settings['analysis_gas']] + ' Absorbance', 
                            fontsize = 10)
        self.ax3.set_xlabel('Wavelength (nm)', fontsize=10)
        
        self.ax4.set_ylabel(gas_choice[settings['analysis_gas']] + ' amt (ppm.m)',
                            fontsize = 10)
        if settings['x_plot'] == 'Number':
            self.ax4.set_xlabel('Spectrum number', fontsize=10)
        if settings['x_plot'] == 'Time':
            self.ax4.set_xlabel('Time (decimal hours)', fontsize=10)
        
        # Create lines to plot data series
        
        # Spectral data
        self.line0, = self.ax0.plot(0, 0, label = 'Spectrum')
        self.line1, = self.ax0.plot(0, 0, label = 'Fit')
        self.ax0.legend(loc = 0)
        
        # Full spectrum
        self.line2, = self.ax1.plot(0, 0)
        
        # Residual
        self.line3, = self.ax2.plot(0, 0)
        
        # SO2 transmittance data
        self.line4, = self.ax3.plot(0, 0, label = 'Meas Abs')
        self.line5, = self.ax3.plot(0, 0, label = 'Synth Abs')
        self.ax3.legend(loc = 0)
        
        # SO2 Time series and error bars
        self.line6, = self.ax4.plot(0, 0)
        
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
#================================== Post Analysis Frame =================================
#======================================================================================== 
#======================================================================================== 





#========================================================================================
#=================================== Create GUI frames ==================================
#========================================================================================

        # Create frames for diferent sections for post analysis     
        setup_frame = tk.LabelFrame(page1, text = 'Setup', font = LARG_FONT)
        setup_frame.grid(row=0, column=0, padx = 10, pady = 10, sticky="ew")
      
#========================================================================================
#================================= Create control inputs ================================
#======================================================================================== 
               
        # Find available flat spectra and form into list
        options = [settings['Spectrometer']]       

        for i, name in enumerate(glob.glob('data_bases/Spectrometer/flat_*')):
            options.append(name[29:-4])
            
        # Create function to turn on the fwd model flag if the spectrometer is changed
        def on_change(event):
            self.build_model_flag = True
        
        # Create entry to select spectrometer
        self.spec_name = tk.StringVar(setup_frame, value = options[0])
        make_input(frame = setup_frame, 
                   text = 'Spectrometer:', 
                   row = 0, column = 0, 
                   var = self.spec_name, 
                   input_type = 'OptionMenu',
                   options = options, 
                   sticky = 'W')
        
        # Create entry to select spectra type
        spec_options = [settings['Spectra Type'],
                        'iFit',
                        'Master.Scope',
                        'Jai Spec',
                        'Spectrasuite',
                        'GSJ',
                        'Ind']
        
        self.spec_type = tk.StringVar(setup_frame, value = spec_options[0])
        make_input(frame = setup_frame, 
                   text = 'Spectra Type:', 
                   row = 1, column = 0, 
                   var = self.spec_type, 
                   input_type = 'OptionMenu',
                   options = spec_options, 
                   sticky = 'W')
        
#========================================================================================
#================================= Create file dialouges ================================
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
        self.specfp_l = tk.Entry(setup_frame, font = NORM_FONT, width = 30,
                                 text = self.spec_ent)
        self.specfp_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W', 
                      columnspan = 2)
        specfp_b = ttk.Button(setup_frame, text="Select Spectra", 
                              command = lambda: spec_fp(self))
        specfp_b.grid(row = 2, column = 2, padx = 5, pady = 5, sticky = 'W')
        
        # File dialouge for darks
        if self.dark_fpaths == ['']:
            message = 'No spectra selected'
        else:
            message = str(len(self.dark_fpaths))+' spectra selected'
        self.dark_ent = tk.StringVar(value = message)
        self.darkfp_l = tk.Entry(setup_frame, font = NORM_FONT, width = 30,
                                 text = self.dark_ent)
        self.darkfp_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W', 
                      columnspan = 2)
        darkfp_b = ttk.Button(setup_frame, text = "Select Darks",
                              command = lambda: dark_fp(self))
        darkfp_b.grid(row = 3, column = 2, padx = 5, pady = 5, sticky = 'W')
        

        
        
     

#========================================================================================         
#========================================================================================
#=============================== Real Time Analysis Frame ===============================
#======================================================================================== 
#======================================================================================== 





#========================================================================================
#=================================== Create GUI frames ==================================
#========================================================================================

        # Create frames for diferent sections for post analysis     
        setup_frame2 = tk.LabelFrame(page2, text='Setup', font = LARG_FONT)
        setup_frame2.grid(row=0, column=0, padx = 10, pady = 10, sticky="ew")

#========================================================================================
#================================= Create control inputs ================================
#========================================================================================
        
        # Create label to display the spectrometer name
        self.c_spec = tk.StringVar(setup_frame2, value = 'Not Connected')
        make_input(frame = setup_frame2, 
                   text = 'Device:', 
                   row = 0, column = 0, 
                   var = self.c_spec, 
                   input_type = 'Label',
                   sticky = ['W', None])
        
        # Integration Time
        self.int_time = tk.DoubleVar(self, value = settings['int_time'])

        make_input(frame = setup_frame2, 
                   text = 'Integration\ntime (ms):', 
                   row = 1, column = 0, 
                   var = self.int_time, 
                   input_type = 'Spinbox',
                   width = 15,
                   vals = [50, 1000],
                   increment = 50)
        
        # Coadds
        self.coadds = tk.DoubleVar(self, value = settings['coadds'])
        make_input(frame = setup_frame2, 
                   text = 'Scans to\nAverage', 
                   row = 2, column = 0, 
                   var = self.coadds, 
                   input_type = 'Spinbox',
                   width = 15,
                   vals = [1, 100])
        
        # Number of darks to get
        self.no_darks = tk.DoubleVar(self, value = settings['no_darks'])
        make_input(frame = setup_frame2, 
                   text = 'Number\nof Darks', 
                   row = 3, column = 0, 
                   var = self.no_darks, 
                   input_type = 'Spinbox',
                   width = 15,
                   vals = [1, 100])
        
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
        
        # Create switch to toggle fitting on or off
        self.toggle_b = tk.Button(setup_frame2, text = 'FITTING OFF', width = 12, 
                                  height = 1, bg = 'red', font = LARG_FONT,
                                  command = lambda: fit_toggle(self, settings))
        self.toggle_b.grid(row=4, column=0, padx=5, pady=5, columnspan=3)
        
        
        
        
        
        
        
 
#========================================================================================       
#========================================================================================
#============================ Program start/stop and analysis ===========================
#========================================================================================
#========================================================================================
               
#========================================================================================
#================================= Create control buttons ===============================
#========================================================================================         
        
        # Frame to hold the buttons
        button_frame = ttk.Frame(quick_frame)
        button_frame.grid(row=0, column=0, padx = 10, pady = 10, columnspan = 5,
                          sticky="ew")        
        
        # Create button to start
        start_b = ttk.Button(button_frame, command = self.begin, text = 'Begin!')
        start_b.grid(row = 0, column = 0, padx = 25, pady = 5)
        
        # Create button to stop
        stop_b = ttk.Button(button_frame, command = lambda: stop(self, settings),
                            text = 'Stop')
        stop_b.grid(row = 0, column = 1, padx = 25, pady = 5)
        
        # Create button for advanced settings
        adv_set_b = ttk.Button(button_frame, text = 'Adv. Settings', 
                            command = lambda: adv_settings(self, settings, 'iFit'))
        adv_set_b.grid(row = 1, column = 0, padx = 25, pady = 5)

        # Create button to save settings
        save_b = ttk.Button(button_frame, text = 'Save Settings', command = self.save)
        save_b.grid(row = 1, column = 1, padx = 25, pady = 5)
               
#========================================================================================
#================================= Progress and Analysis ================================
#======================================================================================== 
        
        # Create progress bar
        self.progress = ttk.Progressbar(quick_frame, orient = tk.HORIZONTAL, length=300,
                                        mode = 'determinate')
        self.progress.grid(row = 1, column = 0, padx = 5, pady = 5, columnspan = 4)
        
        # Create status indicator
        self.status = tk.StringVar(quick_frame, value = 'Standby')
        self.status_e = tk.Label(quick_frame, textvariable = self.status)
        self.status_e.grid(row=1, column=4, padx=5, pady=5, sticky="EW")
        
        # Create ouput for last so2 amount
        self.last_amt = tk.StringVar(self, value = '-')
        make_input(frame = quick_frame, 
                   text = 'Last amt:', 
                   row = 2, column = 0, 
                   var = self.last_amt, 
                   input_type = 'Label',
                   sticky = 'W')
        
        # Create ouput for last so2 error
        self.last_err = tk.StringVar(self, value = '-')
        make_input(frame = quick_frame, 
                   text = '+/-', 
                   row = 2, column = 2, 
                   var = self.last_err, 
                   input_type = 'Label',
                   sticky = 'W')
        
        
        
        
        
        
        
        
        
        


        
#========================================================================================         
#========================================================================================
#===================================== GUI Operations ===================================
#======================================================================================== 
#======================================================================================== 

    # Report exceptions in a new window
    def report_callback_exception(self, *args):
        
        # Report error
        err = traceback.format_exception(*args)
        tkMessageBox.showerror('Exception', err)
        
        # Reset formation of the forward model
        self.build_model_flag = True

    # Close program on 'x' button
    def handler(self):
        
        # Turn on stopping flag
        self.stop_flag = True
        
        # Open save dialouge
        if tkMessageBox.askyesno('Exit', 'Would you like to\nsave the settings?'):
        
            self.save()
            self.quit()
            
        else:
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
            with open(self.notes_fname, 'a') as a:
                a.write(text)
        
        except AttributeError:
            pass
        
        # Force gui to update
        mygui.update(self)  

#========================================================================================
#======================================= Begin iFit =====================================
#========================================================================================
     
    # Function to begin analysis loop
    def begin(self):
        
        # get program mode from notebook frame
        if self.nb.index(self.nb.select()) == 0:
            mode = 'rt'
        else:
            mode = 'post'
        
        # Turn off stopping flag
        self.stop_flag = False
       
        # Create common dictionary
        if self.build_model_flag == True:
            common = {}
        else:
            common = self.common

        # Populate common with other data from the GUI
        common['poly_n']           = int(settings['poly_n']) + 1
        common['ils_width']        = float(settings['ils_width'])
        common['gauss_weight']     = float(settings['gauss_weight'])
        common['ldf']              = float(settings['ldf'])
        common['dark_flag']        = bool(settings['dark_flag'])
        common['flat_flag']        = bool(settings['flat_flag'])
        common['Fit shift']        = str(settings['Fit shift'])
        common['fit_weight']       = str(settings['fit_weight'])
        common['solar_resid_flag'] = str(settings['solar_resid_flag'])

#========================================================================================
#=============================== Build parameter dictionary =============================
#========================================================================================
        
        # Create parameter array
        common['params'] = OrderedDict()

        for i in range(common['poly_n']):
            common['params']['p'+str(i)] = [1.0, 'Fit']
            
        # Add other parameters
        common['params']['shift']     = [settings['shift'],     settings['Fit shift']   ]
        common['params']['stretch']   = [settings['stretch'],   settings['Fit stretch'] ]
        common['params']['ring_amt']  = [settings['ring_amt'],  settings['Fit ring']    ]
        common['params']['so2_amt']   = [settings['so2_amt'],   settings['Fit so2']     ]
        common['params']['no2_amt']   = [settings['no2_amt'],   settings['Fit no2']     ]
        common['params']['o3_amt']    = [settings['o3_amt'],    settings['Fit o3']      ]
        common['params']['bro_amt']   = [settings['bro_amt'],   settings['Fit bro']     ]
        common['params']['ils_width'] = [settings['ils_width'], settings['Fit ILS']     ]
        common['params']['ldf']       = [settings['ldf'],       settings['Fit LDF']     ]
        
        # Make sure first guesses are floats
        for key, val in common['params'].items():
            common['params'][key][0] = float(common['params'][key][0])
                                
        # Create empty last spec
        common['last_spec'] = []
        
        # Save initial guess parameters
        initial_params = common['params'].copy()
        
        # Create residual counter for forming solar residual
        resid_count = 0
        
#========================================================================================
#======================================= Run setup ======================================
#========================================================================================        

        if mode == 'rt':
            setup = rt_setup(self, settings, common, mygui)
            
        if mode == 'post':
            setup = post_setup(self, settings, common, mygui)
            
        if setup == False:
            return
        
        # Read in xsecs and flat spectrum
        if self.build_model_flag:
            common = build_fwd_data(common, settings, self)
        
            
#========================================================================================
#================================== Start Analysis Loop =================================
#========================================================================================

        # Open excel file
        with open(self.out_excel_fname, 'a') as writer:

            # Print output message to begin
            self.print_output('Loop Started\n' +\
                              'Spectrum number ' + str(self.loop)) 
            
            # Create empty arrays to hold the loop number and so2_amt values
            gas = {}
            spec_nos = []
            spec_times = []
            for g in ['so2', 'no2', 'o3', 'bro', 'ring']:
                gas[g + '_amts'] = []
                gas[g + '_errs'] = []
            
            # If forming solar residual create empty array and counter
            if common['solar_resid_flag'] == 'Generate':
                common['solar_resid'] = np.zeros(len(common['grid']))
                resid_count = 0
        
            # Update status bar message
            if mode == 'rt':
                self.status.set('Acquiring')
            else:
                self.status.set('Analysing')
            mygui.update(self)

            # Begin analysis loop
            while not self.stop_flag:
                
#========================================================================================
#=============================== Read and analyse spectra ===============================
#========================================================================================
                
#==================================== Post analysis =====================================
               
                if mode == 'post':
                    
                    try:
                        # Read the spectrum
                        fname = self.spec_fpaths[self.loop]
                        spec_data = read_spectrum(fname, self.spec_type.get())
                        x, y, read_date, read_time, spec_no, read_err = spec_data
                
                    except IndexError:
                        # Stop if reached the end of the spectra
                        self.print_output('Fitting complete')
                        break
                    
                    if read_err[0] == False:
                    
                        # If spectrum was read ok then fit!
                        fit_results = fit_spec(common, [x, y], common['grid'])
                        fit_dict, err_dict, y_data, fit, gas_T, fit_flag = fit_results
                        proceed_flag = True
                        
                    else:
                        # If a read error occured, skip the analysis step
                        self.print_output('Error reading spectrum:\n' + str(read_err[1]))
                        proceed_flag = False
                     
                    # Update progress bar
                    prog = (self.loop+1)/len(self.spec_fpaths) * 100
                    
                
#================================== Real time analysis ==================================
                
                if mode == 'rt':
                    
                    # Measure spectrum while fitting the last if fitting is on
                    [x, y], fit_results, spec_info = rt_analyse(self, settings, common, 
                                                                mygui)
                    spec_no, read_date, read_time, fname = spec_info
                    
                    # Unpack fit results
                    if len(fit_results) != 0:
                        fit_dict, err_dict, y_data, fit, gas_T, fit_flag = fit_results
                        proceed_flag = True
                    
                    # If not fitting, skip the analysis step
                    else:
                        proceed_flag = False
                    
                    # Set progress
                    prog = self.loop
                    
                # Update progress bar
                self.progress['value'] = prog
                    
#========================================================================================
#================================== Look at fit results =================================
#========================================================================================               
                
                if proceed_flag == True:
                
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
    
                    # Check fit quality and update first guess params if required
                    if fit_flag == False:
                        fit_msg = 'Failed'
                        if bool(settings['update_params']) == True:
                            common['params'] = initial_params.copy()
                            self.print_output('Fitting for spectrum '+str(spec_no)+\
                                              ' failed, resetting parameters')
                            
                    elif max_resid > float(settings['good_fit_bound'])/100:
                        fit_msg = 'Bad'
                        if bool(settings['update_params']) == True:
                            common['params'] = initial_params.copy()
                            self.print_output('Fitting for spectrum '+str(spec_no)+\
                                              ' bad, resetting parameters')
                            
                    else:
                        fit_msg = 'Good'
                        if bool(settings['update_params']) == True:
                            # Update first guesses with last fitted params
                            for key, val in fit_dict.items():
                                common['params'][key][0] = val
    
                    # Write results to excel file, starting with spectrum info
                    writer.write(str(fname)     + ',' + \
                                 str(spec_no)   + ',' + \
                                 str(read_date) + ',' + \
                                 str(read_time))          
        
                    # Print fit results and error for each parameter          
                    for key, val in common['params'].items():
    
                        if val[1] == 'Fit':
                            writer.write(','+str(fit_dict[key])+','+str(err_dict[key]))
                        if val[1] in ['Fix', 'Pre-calc', 'File']:
                            writer.write(','+str(common['params'][key][0])+',NaN')
                        if val[1] == 'N/A':
                            writer.write(',NaN,NaN')
                        
                    # Write fit quality and start new line
                    writer.write(',' + fit_msg + '\n')
                    
                    # Add values to array for plotting
                    spec_nos.append(spec_no)
                    try:
                        spec_times.append(hms_to_julian(read_time))
                    except AttributeError:
                        spec_times.append(hms_to_julian(read_time, 
                                                        str_format = '%H:%M:%S'))
                        
                    for parameter in ['so2', 'no2', 'o3', 'bro', 'ring']:
                        
                        # Make dictionary key
                        key1 = parameter + '_amt'
                        key2 = parameter + '_amts'
                        key3 = parameter + '_errs'
                        
                        # Choose conversion factor
                        if parameter == 'ring':
                            conv = 1
                        else:
                            conv = 2.463e15

                        if common['params'][key1][1] == 'Fit':
                            gas[key2].append(fit_dict[key1]/conv)
                            gas[key3].append(err_dict[key1]/conv)
                            
                        if common['params'][key1][1] == 'Fix':
                            gas[key2].append(common['params'][key1][0]/conv)
                            gas[key3].append(0)
                            
                        if common['params'][key1][1] == 'N/A':
                            gas[key2].append(0)
                            gas[key3].append(0)
    
                    # Update quick analysis with values
                    last_amt="{0:0.2f}".format(gas[settings['analysis_gas']+'_amts'][-1])
                    last_err="{0:0.2f}".format(gas[settings['analysis_gas']+'_errs'][-1])
                    self.last_amt.set(last_amt + ' ppm.m')
                    self.last_err.set(last_err + ' ppm.m')
                    
                    # Cut if too long to avoid slowing program
                    if bool(settings['scroll_flag']) == True:
                        lim = int(settings['scroll_spec_no'])
                        if len(spec_nos) > lim:
                            # Find the difference
                            diff = len(spec_nos) - lim
                            
                            # Extract that range
                            spec_nos = spec_nos[diff:]
                            spec_times = spec_times[diff:]
                            for m in gas:
                                gas[m] = gas[m][diff:]
                                
                    # Select whether to show so2 time series in time or number
                    if settings['x_plot'] == 'Number':
                        x_plot = spec_nos
                    if settings['x_plot'] == 'Time':
                        x_plot = spec_times
                                
#========================================================================================
#====================================== Update plot =====================================
#========================================================================================

                # Replot data
                if bool(settings['Show Graphs']) == True:            
                    
                    if proceed_flag == True:
                        
                        # Get selected transmittance data
                        meas_abs = gas_T['meas_abs_' + settings['analysis_gas']]
                        synth_abs = gas_T['synth_abs_' + settings['analysis_gas']]
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
                        t_lo = min(meas_abs) - abs((0.1*min(meas_abs)))
                        t_hi = max(meas_abs) + (0.1*max(meas_abs)) 
                        s_lo = min(synth_abs) - abs((0.1*min(synth_abs)))
                        s_hi = max(synth_abs) + (0.1*max(synth_abs))
                        t_lo, t_hi = min([t_lo, s_lo]), max([t_hi, s_hi])
                        
                        grid = common['grid']
                        # Build data array to pass to graphing function
                        #                 x data  y data     x limits     y limits
                        data = np.array(([grid,   y_data,    'auto'     , [y_lo,y_hi]],
                                         [grid,   fit,       'auto'     , [y_lo,y_hi]],
                                         [x   ,   y,         'auto'     , 'auto'     ],
                                         [grid,   resid,     'auto'     , 'auto'     ],
                                         [grid,   meas_abs,  'auto'     , [t_lo,t_hi]],
                                         [grid,   synth_abs, 'auto'     , [t_lo,t_hi]],
                                         [x_plot, gas_amts,  'auto'     , 'auto'     ]))
                    
                    elif mode == 'rt':
                        
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
                               
                # Add to the count cycle
                self.loop += 1
                
                # Keep common in memory
                self.common = common
                
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
#====================================== Save Settings ===================================
#========================================================================================
    
    # Function to save setting to the ifit_settings.txt file        
    def save(self):
        
        # Create or overright settings file
        with open('data_bases/ifit_settings.txt', 'w') as w:
            
            # Save each setting from the gui into settings
            settings['Spectrometer']      = str(self.spec_name.get())       
            settings['Spectra Type']      = str(self.spec_type.get())       
            settings['int_time']          = str(self.int_time.get())         
            settings['coadds']            = str(self.coadds.get())
            settings['no_darks']          = str(self.no_darks.get())
            settings['Spectra Filepaths'] = str(self.spec_fpaths)
            settings['Dark Filepaths']    = str(self.dark_fpaths)
            
            # Add all of the settings dictionary
            for s in settings:
               w.write(s + ';' + str(settings[s]) + ';' + str(type(settings[s])) + '\n')
                
        self.print_output('Settings saved')

    
# Run the App!
if __name__ == '__main__':    
    mygui().mainloop()
