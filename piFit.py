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
from queue import Queue
from threading import Thread
import datetime
from seabreeze.cseabreeze.wrapper import SeaBreezeError
import matplotlib.pyplot as plt
from collections import OrderedDict
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from ifit_lib.build_fwd_data import build_fwd_data
from ifit_lib.fit import fit_spec
from ifit_lib.acquire_spectrum import acquire_spectrum
from ifit_lib.file_control import make_directory, make_csv_file
from ifit_lib.gui_funcs import adv_settings, fit_toggle, stop, connect_spec, test_spec, \
                               update_int_time, read_darks, read_settings

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
        tk.Tk.wm_title(self, 'piFit-2-3')
        #tk.Tk.iconbitmap('data_bases/icon.ico')
        
        # Create frame to hold control variables
        control_frame = ttk.Frame(self, relief = 'groove')
        control_frame.grid(column=0, row=0, padx=10, pady=10, rowspan=2, sticky='NW')
        
        # Create frame to hold graphs
        graph_frame = ttk.Frame(self, relief = 'groove')
        graph_frame.grid( column=1,row=0, padx=10, pady=10, sticky="NW")
        graph_frame.columnconfigure(index=0, weight=1)
        graph_frame.rowconfigure(index = 0, weight = 1)
        
        # Create frame to hold text output
        text_frame = ttk.Frame(self, relief = 'groove')
        text_frame.grid(column=1, row=1, padx=10, pady=10, sticky="NW")
        
        # Frame for quick analysis
        quick_frame = tk.Frame(control_frame, relief = 'groove')
        quick_frame.grid(column=0, row=2, padx=10, pady=10, sticky="NW")
        
        mygui.columnconfigure(index=1, weight=1, self = self)
        mygui.rowconfigure(index = 5, weight = 1, self = self)
        
#========================================================================================
#====================================Create text output==================================
#========================================================================================      

        # Create button for advanced settings
        adv_set_b = ttk.Button(text_frame, text = 'Adv. Settings', width = 15, 
                            command=lambda: adv_settings(self, settings, 'piFit'))
        adv_set_b.grid(row = 0, column = 0, padx = 40, pady = 5)

        # Create button to save settings
        save_b = ttk.Button(text_frame, text = 'Save Settings', width = 15,
                            command = self.save)
        save_b.grid(row = 0, column = 1, padx = 40, pady = 5)  
                 
        # Build text box
        self.text_box = tkst.ScrolledText(text_frame, width = 45, height = 5)
        self.text_box.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W',
                           columnspan = 4)
        self.text_box.insert('1.0', 'Welcome to piFit! Written by Ben Esse\n\n')
        
#========================================================================================
#==============================Create quick analysis outputs=============================
#========================================================================================
        
        # Create progress bar
        self.progress = ttk.Progressbar(quick_frame, orient = tk.HORIZONTAL, length=200,
                                        mode = 'determinate')
        self.progress.grid(row = 0, column = 0, padx = 5, pady = 5, columnspan = 4)
        
        # Create status indicator
        self.status = tk.StringVar(quick_frame, value = 'Standby')
        self.status_e = tk.Label(quick_frame, textvariable = self.status)
        self.status_e.grid(row=0, column=5, padx=5, pady=5, sticky="EW")
        
        # Create ouput for last so2 amount
        self.last_so2_amt = tk.StringVar(self, value = '-')
        last_so2_amt_l = tk.Label(quick_frame, text = 'Last amt:', font = NORM_FONT)
        last_so2_amt_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        last_so2_amt_e = tk.Label(quick_frame, textvariable = self.last_so2_amt)
        last_so2_amt_e.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        # Create ouput for last so2 error
        self.last_so2_err = tk.StringVar(self, value = '-')
        last_so2_err_l = tk.Label(quick_frame, text = '+/-', 
                                  font = NORM_FONT)
        last_so2_err_l.grid(row = 1, column = 2, pady = 5, sticky = 'W')
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
            settings = read_settings('data_bases/pifit_settings.txt', settings)
   
        except FileNotFoundError:
            self.print_output('No settings file found, reverting to default')
            settings['wave_start']        = 305
            settings['wave_stop']         = 318
            settings['int_time']          = 100
            settings['coadds']            = 10
            settings['no_darks']          = 10
            settings['ils_width']         = 0.52
            settings['gauss_weight']      = 1.0
            settings['Fit ILS']           = 'Fix'
            settings['ldf']               = 0.0
            settings['Fit LDF']           = 'N/A'
            settings['dark_flag']         = True
            settings['flat_flag']         = True
            settings['update_params']     = True
            settings['good_fit_bound']    = 10
            settings['fit_weight']        = 'None'
            settings['Show Graphs']       = True
            settings['Show Error Bars']   = 0
            settings['analysis_gas']      = 'SO2'
            settings['scroll_flag']       = True
            settings['scroll_spec_no']    = 200
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
            settings['model_pad']         = 3
            settings['sol_path']          = 'data_bases/gas data/sao2010.txt'
            settings['ring_path']         = 'data_bases/gas data/qdoas_ring.dat'
            settings['so2_path']          = 'data_bases/gas data/SO2_293K.dat'
            settings['no2_path']          = 'data_bases/gas data/No2_223l.dat'
            settings['o3_path']           = 'data_bases/gas data/O3_xsec.dat'
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
        
        # Create graph label dictionary
        self.graph_labels = {'Spectrum':   ['Wavelength (nm)', 'Intensity (arb)'], 
                             'Fit':        ['Wavelength (nm)', 'Intensity (arb)'],
                             'Residual':   ['Wavelength (nm)', 'Residual'],
                             'Absorbance': ['Wavelength (nm)', 'Absorbance'],
                             'Gas amount': ['Spectrum Number', 'Amount (ppm.m)']}    
                    
        # Create figure to hold the graphs
        plt.rcParams.update({'font.size': 8} )
        self.fig = plt.figure(figsize = (4,2.5))
        
        # Create plot axes
        self.ax = self.fig.add_subplot(111)
        
        # Set axis labels
        self.ax.set_xlabel(self.graph_labels[settings['graph_view']][0])
        self.ax.set_ylabel(self.graph_labels[settings['graph_view']][1])
        
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
        setup_frame = tk.LabelFrame(control_frame, text='Spectrometer Setup', 
                                     font = LARG_FONT)
        setup_frame.grid(row=0, column=0, padx = 10, pady = 10, sticky="ew")
        
        button_frame = tk.LabelFrame(control_frame, text='Control', font = LARG_FONT)
        button_frame.grid(row=1, column=0, padx = 10, pady = 10, sticky="ew")

#========================================================================================
#==================================Create control inputs=================================
#========================================================================================
        
        # Create label to display the spectrometer name
        self.c_spec = tk.StringVar(setup_frame, value = 'None connected')
        c_spec_l = ttk.Label(setup_frame, text="Device: ", font = NORM_FONT)
        c_spec_l.grid(row=0, column=0, pady=5, padx=5, sticky='W')
        c_spec_e = ttk.Entry(setup_frame, textvariable = self.c_spec, width = 10)
        c_spec_e.grid(row = 0, column = 1, padx = 5, pady = 5)
        
        # Integration Time
        self.int_time = tk.IntVar(setup_frame, value = settings['int_time'])
        int_time_vals = [1,2,3,4,5,10,20,30,40,50,100,200,300,400,
                         500,600,700,800,900,1000]
        int_time_l = tk.Label(setup_frame, text = 'Int. Time:', font = NORM_FONT)
        int_time_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        int_time_sb = tk.Spinbox(setup_frame, textvariable=self.int_time,
                                   values = int_time_vals, width = 10)
        int_time_sb.grid(row = 1, column = 1, padx = 5, pady = 5)
        self.int_time.set(settings['int_time'])
        
        # Coadds
        self.coadds = tk.IntVar(setup_frame, value = settings['coadds'])
        coadd_vals=[1,2,3,4,5,10,20,30,40,50,100,200,300,400,500,1000]
        coadds_l = tk.Label(setup_frame, text = 'Coadds:', font = NORM_FONT)
        coadds_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        coadds_sb = tk.Spinbox(setup_frame, textvariable = self.coadds, 
                               values = coadd_vals, width = 10)
        coadds_sb.grid(row = 2, column = 1, padx = 5, pady = 5)
        self.coadds.set(settings['coadds'])
        
        # Number of darks to get
        self.no_darks = tk.IntVar(setup_frame, value = settings['no_darks'])
        no_dark_vals = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        no_darks_l = tk.Label(setup_frame, text = 'No. Darks:', font = NORM_FONT)
        no_darks_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
        no_darks_sb = tk.Spinbox(setup_frame, textvariable = self.no_darks,
                                 values = no_dark_vals, width = 10)
        no_darks_sb.grid(row = 3, column = 1, padx = 5, pady = 5)
        self.no_darks.set(settings['no_darks'])
        
        # Create button to connect to spectrometer
        connect_spec_b = ttk.Button(setup_frame, text = 'Connect', width = 10,
                                    command = lambda: connect_spec(self, settings))
        connect_spec_b.grid(row = 0, column = 2, pady = 5)
        
        # Create button to update integration time
        update_int_time_b = ttk.Button(setup_frame, text = 'Update', width = 10,
                                       command = lambda: update_int_time(self, settings))
        update_int_time_b.grid(row = 1, column = 2, pady = 5)
        
        # Create button to read a single spectrum
        test_spec_b = ttk.Button(setup_frame, text = 'Test Spec', width = 10,
                                 command = lambda: test_spec(self, settings, mygui, 
                                                             self.line0, self.ax))
        test_spec_b.grid(row = 2, column = 2, pady = 5)
        
        # Create button to read darks
        read_darks_b = ttk.Button(setup_frame, text = 'Read Darks', width = 10,
                                  command = lambda: read_darks(self, settings, mygui, 
                                                               self.line0, self.ax))
        read_darks_b.grid(row = 3, column = 2, pady = 5)
        
#========================================================================================
#==============================Create start and stop buttons=============================
#========================================================================================         
        
        # Create button to start
        start_aq_b = ttk.Button(button_frame, text = 'Begin!', command = self.begin,
                                width = 10)
        start_aq_b.grid(row = 0, column = 0, padx = 5, pady = 5, columnspan=2)
        
        # Create button to stop
        stop_aq_b = ttk.Button(button_frame, text = 'Stop', width = 10,
                               command = lambda: stop(self, settings))
        stop_aq_b.grid(row = 0, column = 1, padx = 5, pady = 5, columnspan=2)
        
        # Create switch to toggle fitting on or off
        self.toggle_button = tk.Button(button_frame, text = 'FITTING OFF', width = 12,
                                       height = 1, bg = 'red', font = LARG_FONT,
                                       command = lambda: fit_toggle(self, settings))
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
                
        self.graph_view = tk.StringVar(button_frame, value = settings['graph_view'])
        graph_view_l = tk.Label(button_frame, text = 'Graph View:', font = NORM_FONT)
        graph_view_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        graph_view_m = ttk.OptionMenu(button_frame, self.graph_view, *graph_options,
                                      command = update)
        graph_view_m.grid(row = 1, column = 1, padx = 5, pady = 5)
        





        
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
    def begin(self):
        
        # Turn off stopping flag
        settings['stop_flag'] = False
       
        # Create common dictionary
        common = {}
        
        # Populate common with other data from the GUI
        common['wave_start']       = float(settings['wave_start'])
        common['wave_stop']        = float(settings['wave_stop'])
        common['poly_n']           = int(settings['poly_n']) + 1
        common['ils_width']        = float(settings['ils_width'])
        common['gauss_weight']     = float(settings['gauss_weight'])
        common['ldf']              = float(settings['ldf'])
        common['dark_flag']        = bool(settings['dark_flag'])
        common['flat_flag']        = bool(settings['flat_flag'])
        common['solar_resid_flag'] = str(settings['solar_resid_flag'])
        common['Fit shift']       = str(settings['Fit shift'])
        common['fit_weight']       = str(settings['fit_weight'])

        # Turn of dark flag if in real time and no darks have been taken
        if settings['rt_dark_flag'] == False:
            common['dark_flag'] = 0
            self.print_output('WARNING! No dark spectra aquired!')

#========================================================================================
#================================Build parameter dictionary==============================
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
        common['last_spec'] = np.array([0])
        
        # Save initial guess parameters
        initial_params = common['params'].copy()

#========================================================================================
#=================================Read in xsecs and flat=================================
#========================================================================================

        # Update status
        self.status.set('Building Model')
        mygui.update(self)

        # Get spectrometer serial number to get flat and ILS
        common['spec_name'] = str(self.c_spec.get())
        
        # Load fitting data files
        common = build_fwd_data(common, settings, self)

#========================================================================================
#===================================Build dark spectrum==================================
#========================================================================================

        if common['dark_flag'] == True:
            
            # Assign dark spectrum
            common['dark'] = self.dark_spec

#========================================================================================
#========================Read test spectrum to get wavelength grid=======================
#========================================================================================
 
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
        common['fit_idx'] = np.where(np.logical_and(common['wave_start'] <= x, 
                                                    x <= common['wave_stop']))
        grid = x[common['fit_idx']]
        
        # Find stray light window
        common['stray_idx'] = np.where(np.logical_and(280 <= x, x <= 290))
        
        # If no stray light pixels available, turn off the flag
        if len(common['stray_idx'][0]) == 0:
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
            
        # Create new output folder
        if settings['create_out_flag'] == True:
           
            # Reset loop counter
            settings['loop'] = 0
            
            # Create filename for output file
            out_excel_fname = settings['folder'] + 'iFit_out.csv'
            
            # Open excel file and write header line
            make_csv_file(out_excel_fname, common)
            
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
                    t1 = Thread(target = acquire_spectrum, args=(self,
                                                                 settings['spec'],
                                                                 settings['int_time'],
                                                                 int(self.coadds.get()),
                                                                 True,
                                                                 True,
                                                                 result_queue))
                    
                    t2 = Thread(target = fit_spec, args = (common,
                                                           [x, common['last_spec']],
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
                    x, y, header, t = acquire_spectrum(self, settings['spec'], 
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
                if int(settings['Show Graphs']) == 1 and wait_flag == False:            
                    
                    if self.toggle_button.config('text')[-1]=='FITTING ON':
                        
                        # Get selected transmittance data
                        gas_abs  = gas_T[settings['analysis_gas'] + '_tran']
                        gas_spec = gas_T[settings['analysis_gas'] + '_spec']
                        gas_amts = gas[settings['analysis_gas'] + '_amts']
                        
                        graph_view = self.graph_view.get()
                        
                        if graph_view == 'Spectrum':
                            line0_x_data = grid
                            line0_y_data = y_data
                            line1_x_data = grid
                            line1_y_data = fit
                            
                        if graph_view == 'Fit':
                            line0_x_data = x
                            line0_y_data = y
                            line1_x_data = grid[0]
                            line1_y_data = y_data[0]
                            
                        if graph_view == 'Residual':
                            line0_x_data = grid
                            line0_y_data = resid
                            line1_x_data = grid[0]
                            line1_y_data = resid[0]
                            
                        if graph_view == 'Absorbance':
                            line0_x_data = grid
                            line0_y_data = gas_abs
                            line1_x_data = grid
                            line1_y_data = gas_spec
                            
                        if graph_view == 'Gas amount':
                            line0_x_data = spec_nos
                            line0_y_data = gas_amts
                            line1_x_data = spec_nos[0]
                            line1_y_data = gas_amts[0]
                            
                        # Update graph with correct data
                        self.line0.set_data(line0_x_data, line0_y_data)
                        self.line1.set_data(line1_x_data, line1_y_data)
                             
                    
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
#=======================================Save Settings====================================
#========================================================================================
    
    # Function to save setting to the ifit_settings.txt file        
    def save(self):
        
        # Create or overright settings file
        with open('data_bases/pifit_settings.txt', 'w') as w:
            
            # Add all of the settings dictionary
            for s in settings:
                w.write(s + ';' + str(settings[s]) + '\n')
                
        self.print_output('Settings saved')
    
# Run the App!
if __name__ == '__main__':    
    mygui().mainloop()
