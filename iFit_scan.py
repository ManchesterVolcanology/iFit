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
from ifit_lib.fit import fit_spec
from ifit_lib.read_binary_block import read_binary_block
from ifit_lib.update_graph import update_graph
from ifit_lib.gui_funcs import adv_settings, stop, read_settings

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
        self.text_box = tkst.ScrolledText(text_frame, width = 50, height = 8)
        self.text_box.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W',
                           columnspan = 4)
        self.text_box.insert('1.0', 'Welcome to iFit_scan! Written by Ben Esse\n\n')    
        
        # Create button to start
        start_b = ttk.Button(text_frame, text = 'Begin!', command = self.begin)
        start_b.grid(row = 0, column = 0, padx = 5, pady = 5)
        
        # Create button to stop
        stop_b = ttk.Button(text_frame, text = 'Stop', 
                            command = lambda: stop(self, settings))
        stop_b.grid(row = 1, column = 0, padx = 5, pady = 5)
        
        # Create button to save settings
        save_b = ttk.Button(text_frame, text = 'Save Settings', command = self.save)
        save_b.grid(row = 0, column = 2,  padx = 5, pady = 5)
        
        # Create button to open advanced settings
        adv_set_b = ttk.Button(text_frame, text = 'Adv. Settings', 
                               command = lambda: adv_settings(self, settings, 'iFit'))
        adv_set_b.grid(row = 1, column = 2,  padx = 5, pady = 5)
        
#========================================================================================
#===================================Set program settings=================================
#========================================================================================

        # Create settings dictionary
        global settings
        settings = {}
        
        # Read in settings file
        try:
            settings = read_settings('data_bases/ifit_scan_settings.txt', settings)
   
        except FileNotFoundError:
            self.print_output('No settings file found, reverting to origional')
            settings['wave_start']        = 305
            settings['wave_stop']         = 318
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
        self.ax4.set_xlabel('Scan Angle (degrees)', fontsize=10)
        
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
        self.canvas.draw()
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
#==============================Create quick analysis outputs=============================
#========================================================================================
               
        # Create progress bar
        self.progress = ttk.Progressbar(quick_frame, orient = tk.HORIZONTAL, length=350,
                                        mode = 'determinate')
        self.progress.grid(row = 0, column = 0, padx = 5, pady = 5, columnspan = 4)
        
        
        # Create loop counter for scans
        self.scan_count = tk.StringVar(self, value = '1 / 1')
        scan_count_l = tk.Label(quick_frame, text = 'Scan:', font = NORM_FONT)
        scan_count_l.grid(row = 0, column = 4, padx = 5, pady = 5, sticky = 'W')
        scan_count_e = tk.Label(quick_frame, textvariable = self.scan_count)
        scan_count_e.grid(row = 0, column = 5, padx = 5, pady = 5)
        
        # Create ouput for last so2 amount
        self.last_so2_amt = tk.DoubleVar(self, value = 0)
        last_so2_amt_l = tk.Label(quick_frame, text = 'Last amt:', 
                                  font = NORM_FONT)
        last_so2_amt_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        last_so2_amt_e = tk.Label(quick_frame, textvariable = self.last_so2_amt)
        last_so2_amt_e.grid(row = 1, column = 1, padx = 5, pady = 5)
        
        # Create ouput for last so2 error
        self.last_so2_err = tk.DoubleVar(self, value = 0)
        last_so2_err_l = tk.Label(quick_frame, text = 'Last error:', 
                                  font = NORM_FONT)
        last_so2_err_l.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
        last_so2_err_e = tk.Label(quick_frame, textvariable = self.last_so2_err)
        last_so2_err_e.grid(row = 1, column = 3, padx = 5, pady = 5)
        
        


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
            with open(settings['notes_fname'], 'a') as a:
                a.write(text)
        
        except KeyError:
            pass
        
        # Force gui to update
        mygui.update(self)
        
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
            
            # Add all of the settings dictionary
            for s in settings:
                w.write(s + ';' + str(settings[s]) + '\n')
            
            try:
                w.write('Scan Filepaths;' + str(self.scan_fpaths) + '\n')
            except AttributeError:
                w.write('Scan Filepaths; \n') 
                
        self.print_output('Settings saved')  
        
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
        common['wave_start']       = float(settings['wave_start'])
        common['wave_stop']        = float(settings['wave_stop'])
        common['poly_n']           = int(settings['poly_n']) + 1
        common['ils_width']        = float(settings['ils_width'])
        common['gauss_weight']     = float(settings['gauss_weight'])
        common['ldf']              = float(settings['ldf'])
        common['dark_flag']        = bool(settings['dark_flag'])
        common['flat_flag']        = bool(settings['flat_flag'])
        common['solar_resid_flag'] = str(settings['solar_resid_flag'])
        common['Fit shift']        = str(settings['Fit shift'])
        common['fit_weight']       = str(settings['fit_weight'])

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
#========================Read test spectrum to get wavelength grid=======================
#========================================================================================
        
        # Get filepaths to scan files
        scan_files = self.scan_fpaths
        
        # Read in first scan in array
        err,x,header,info_block,spec_block = read_binary_block(scan_files[0])
        
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
#=================================Read in xsecs and flat=================================
#========================================================================================

        # Get spectrometer serial number to get flat and ILS
        common['spec_name'] = 'I2J5046'#str(self.spec_name.get())
        
        # Load fitting data files
        common = build_fwd_data(common, settings, self)
 
#========================================================================================
#===================================Create ouput folder==================================
#========================================================================================
            
        # Reset loop counter
        settings['loop'] = 1
        
        # Get date
        meas_date = scan_files[0].split('/')[-1][:6]

        # Create filepath to directory to hold program outputs
        results_folder = 'Results/iFit_scan/' + str(meas_date) + '/'
        
        # Create folder if it doesn't exist
        if not os.path.exists(results_folder):
                os.makedirs(results_folder)                     

#========================================================================================
#===================================Start Analysis Loop==================================
#========================================================================================

        # Initiate loop counter display
        msg = '1 / ' + str(len(scan_files))
        self.scan_count.set(msg)

        # Output initiation message
        self.print_output('Loop Started\nSpectrum number ' + str(settings['loop']))

        for fpath in scan_files:
            
            # Create empty arrays to hold the loop number and so2_amt values
            gas = {}
            spec_nos = []
            gas['SO2_amts'] = []
            gas['SO2_errs'] = []
            gas['O3_amts']  = []
            gas['O3_errs']  = []
            gas['BrO_amts'] = []
            gas['BrO_errs'] = []
            gas['Ring_amts'] = []
            gas['Ring_errs'] = []
            
            # End loop if finished
            if settings['stop_flag'] == True:
                break

            # Read in scan block
            err, x, header, info_block, spec_block = read_binary_block(fpath)
            
            if err == 0:
            
                # First spectrum is the dark
                common['dark'] = spec_block[:,0]
                
                # Create output file
                fname = fpath.split('/')[-1]
                out_file = results_folder + fname[:-4] + '_output.csv'
                
                with open(out_file, 'w') as w:
                    w.write('Scan no,Time,Motor Position,Angle,so2 (ppm.m),so2 error,')
                            
                    for i in common['params'].keys():
                        w.write(i + ',' + i + '_e,')
                        
                    # Starrt new line 
                    w.write('\n')
                
                    # Analyse rest (106 spectra to a scan, first is dark, last is crap)
                    for n in range(1, len(spec_block.T)):
                        
                        # End loop if finished
                        if settings['stop_flag'] == True:
                            break
    
                        # Load spectrum
                        y = spec_block[:,n]

                        # Fit the spectrum
                        fit_data = fit_spec(common, [x, y], grid)
                        fit_dict, err_dict, y_data, fit, gas_T, fit_flag = fit_data    
                        
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
                        if common['params']['ring_amt'][1] == 'Fit':
                            gas['Ring_amts'].append(fit_dict['ring_amt'])
                            gas['Ring_errs'].append(err_dict['ring_amt'])
                            
                        if common['params']['so2_amt'][1] == 'Fix':
                            gas['SO2_amts'].append(common['params']['so2_amt'][0]/2.463e15)
                            gas['SO2_errs'].append(0)
                        if common['params']['o3_amt'][1] == 'Fix':
                            gas['O3_amts'].append(common['params']['o3_amt'][0]/2.463e15)
                            gas['O3_errs'].append(0)
                        if common['params']['bro_amt'][1] == 'Fix':
                            gas['BrO_amts'].append(common['params']['bro_amt'][0]/2.463e15)
                            gas['BrO_errs'].append(0)
                        if common['params']['ring_amt'][1] == 'Fix':
                            gas['Ring_amts'].append(common['params']['ring_amt'][0])
                            gas['Ring_errs'].append(0)
                            
                        if common['params']['so2_amt'][1] == 'N/A':
                            gas['SO2_amts'].append(0)
                            gas['SO2_errs'].append(0)
                        if common['params']['o3_amt'][1] == 'N/A':
                            gas['O3_amts'].append(0)
                            gas['O3_errs'].append(0)
                        if common['params']['bro_amt'][1] == 'N/A':
                            gas['BrO_amts'].append(0)
                            gas['BrO_errs'].append(0)
                        if common['params']['ring_amt'][1] == 'N/A':
                            gas['Ring_amts'].append(0)
                            gas['Ring_errs'].append(0)
    
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
                            if settings['update_params'] == True:
                                common['params'] = initial_params.copy()
                                self.print_output('Fitting for spectrum '+str(spec_no)+\
                                                  ' failed, resetting parameters')
                                
                        elif max_resid > float(settings['good_fit_bound'])/100:
                            fit_msg = 'Bad'
                            if settings['update_params'] == True:
                                common['params'] = initial_params.copy()
                                self.print_output('Fitting for spectrum '+str(spec_no)+\
                                                  ' bad, resetting parameters')
                                
                        else:
                            fit_msg = 'Good'
                            if settings['update_params'] == True:
                                # Update first guesses with last fitted params
                                for key, val in fit_dict.items():
                                    common['params'][key][0] = val
                                
                        # Print fit results and error for each parameter          
                        for key, val in common['params'].items():
    
                            if val[1] == 'Fit':
                                w.write(','+str(fit_dict[key])+','+str(err_dict[key]))
                            if val[1] == 'Fix':
                                w.write(','+str(common['params'][key][0])+',NaN')
                            if val[1] == 'N/A':
                                w.write(',NaN,NaN')
                            
                        # Write fit quality and start new line
                        w.write(',' + fit_msg + '\n')
                
#========================================================================================
#=======================================Update plot======================================
#========================================================================================
        
                        # Replot data
                        if settings['Show Graphs']:    
                            
                            
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
                        prog = ((n+1) / len(spec_block.T)) * 100
                        self.progress['value'] = prog
                        
                        # Force gui to update
                        mygui.update(self)
                        
            else:
                self.print_output('Error in file ' + str(settings['loop']))  
           
            # Update loop counter display
            msg = str(settings['loop']) + ' / ' + str(len(scan_files))
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

          
# Tkinter stuff 
if __name__ == '__main__':    
    mygui().mainloop()
