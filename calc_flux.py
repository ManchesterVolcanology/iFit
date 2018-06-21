# -*- coding: utf-8 -*-
"""
Created on Thu May 31 16:36:22 2018

@author: mqbpwbe2
"""

# Import required libraries
import numpy as np
import os
import traceback
import tkinter.messagebox as tkMessageBox
import tkinter.scrolledtext as tkst
import matplotlib
matplotlib.use('TkAgg')
from tkinter import ttk
import tkinter as tk
import datetime
from tkinter import filedialog as fd
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
from matplotlib.widgets import Slider
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from pandas import read_csv

from ifit_lib.find_nearest import extract_window
from ifit_lib.read_gps import read_txt_gps, gps_vector, haversine
from ifit_lib.center_of_grav import cog

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
        self.protocol("WM_DELETE_WINDOW",self.handler)
        
        # Button Style
        ttk.Style().configure('TButton', width = 20, height = 20, relief = 'flat')
        
        # Add a title and icon
        tk.Tk.wm_title(self, 'calc_flux 2.0')
        tk.Tk.iconbitmap(self, default = 'data_bases/icon.ico')

        # Create control Frame
        cont_frame = ttk.Frame(self, relief = 'groove')
        cont_frame.grid(row=0, column=0, padx=10, pady=10, rowspan=2, sticky="NW")
        
        # Create frame to hold graphs
        graph_frame = ttk.Frame(self, relief = 'groove')
        graph_frame.grid(row=1, column=1, padx=10, pady=0, rowspan=10, sticky="N")
        graph_frame.columnconfigure(index=0, weight=1)
        graph_frame.rowconfigure(index = 0, weight = 1)
        
        # Create frame to hold text output
        text_frame = ttk.Frame(self, relief = 'groove')
        text_frame.grid(row=2, column=0, padx=10, pady=10, sticky="NW")
        
        # Create global common dictionary to hold settings
        global common
        common = {}
        
        # Create loop counter to hold traverse number
        common['loop'] = 1
        
        # Create array to hold fluxes
        common['fluxes'] = []
        
        # Read in settings file
        choice = ['--select--']
        volc_data = {}
        try:
            with open('data_bases/volcano_data.txt', 'r') as r:
                
                # Read header line
                r.readline()
                
                # Read data
                data = r.readlines()
                
                # Unpack and save to dictionary
                for line in data:
                    
                    # Read each line
                    v_name, v_lon, v_lat, t_diff = line.strip().split(';')
                    
                    # Remove white space from parameters
                    v_name = v_name.strip()
                    v_lon  = v_lon.strip()
                    v_lat  = v_lat.strip()
                    t_diff = t_diff.strip()
                    
                    # Add to choise array and populate volcano dictionary
                    choice.append(v_name)
                    volc_data[v_name] = [v_lon, v_lat, t_diff]
                    
        except  FileNotFoundError:
            print('No input file found')
            
            # Set all values to zero
            volc_data['--select--'] = [0,0,0]
              
#========================================================================================
#====================================Create text output==================================
#========================================================================================     
        
        # Build text box
        self.text_box = tkst.ScrolledText(text_frame, width = 60, height = 10)
        self.text_box.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W',
                           columnspan = 4)
        self.text_box.insert('1.0', 'Welcome to calc_flux! Written by Ben Esse\n\n')
        
#========================================================================================
#=====================================Build Graph========================================
#========================================================================================
        
        # Create figure to hold the graphs
        plt.rcParams.update({'font.size': 8} )
        self.fig = plt.figure(figsize = (7,5))
        
        # Create the canvas to hold the graph in the GUI
        self.canvas = FigureCanvasTkAgg(self.fig, graph_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=0, column=0, columnspan = 3, padx=10, 
                                         pady = 10, sticky = 'NW')
        
        # Create plot axes
        self.ax0 = self.fig.add_subplot(111)
        self.fig.subplots_adjust(bottom=0.3)
        
        self.line0, = self.ax0.plot(0, 0, 'g', lw=1.5)
        self.vline0 = self.ax0.axvline(0)
        self.vline1 = self.ax0.axvline(1)
        self.fill   = self.ax0.fill_between([0], [0], [1])
        self.ax0.set_xlabel('Julian Time (Fraction of Day)', fontsize=10)
        self.ax0.set_ylabel('SO2 column amount (ppm.m)', fontsize=10)
        
        # Create sliders
        self.ax_lo = self.fig.add_axes([0.12, 0.15, 0.78, 0.03])
        self.slide_lo = Slider(self.ax_lo, 'Low Bound', 0, 1, valinit = 0.05)
        self.ax_hi = self.fig.add_axes([0.12, 0.1, 0.78, 0.03])
        self.slide_hi = Slider(self.ax_hi, 'High Bound', 0, 1, valinit = 0.95)
        
        def update(val):
            
            # Get slider values
            pos = self.slide_lo.val, self.slide_hi.val
            
            # Get graph limits
            lims = self.ax0.get_xlim()
            
            # Caclulate position from limits and slider value
            self.t0 = lims[0] + (lims[1] - lims[0]) * pos[0]
            self.t1 = lims[0] + (lims[1] - lims[0]) * pos[1]
            self.vline0.set_xdata(self.t0)
            self.vline1.set_xdata(self.t1)
            
        self.slide_lo.on_changed(update)
        self.slide_hi.on_changed(update)
        
        # Make it look nice
        #plt.tight_layout()
        
        # Add matplotlib toolbar above the plot canvas
        toolbar_frame = tk.Frame(self, bg = 'black')  
        toolbar_frame.grid(row=0,column=1, sticky = 'NW', padx = 6)                             
        toolbar = NavigationToolbar2TkAgg(self.canvas, toolbar_frame)
        toolbar.update()
        
#========================================================================================
#==================================Traverse controls=====================================
#========================================================================================

        # Create frame
        trav_frame = tk.LabelFrame(cont_frame, text='Traverse Settings', font=LARG_FONT)
        trav_frame.grid(row=0, column=0, padx=10, pady=6, columnspan = 2, sticky="NW")
        
        # Create inputs for iFit output file
        self.ifit_path = tk.StringVar(value = 'No file selected')
        ifit_path_l = tk.Label(trav_frame, text = 'iFit file:', font = NORM_FONT)
        ifit_path_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        ifit_path_e = tk.Entry(trav_frame, font = NORM_FONT, width = 20, 
                               text = self.ifit_path)
        ifit_path_e.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = 'W')
        ifit_path_b = ttk.Button(trav_frame, text="Select iFit output", 
                                 command = lambda: self.get_fp(self.ifit_path))
        ifit_path_b.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = 'W')
        
        # Create inputs for gps file
        self.gps_path = tk.StringVar(value = 'No file selected')
        gps_path_l = tk.Label(trav_frame, text = 'GPS file:', font = NORM_FONT)
        gps_path_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        gps_path_e = tk.Entry(trav_frame, font = NORM_FONT, width = 20, 
                              text = self.gps_path)
        gps_path_e.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = 'W')
        gps_path_b = ttk.Button(trav_frame, text="Select GPS file", 
                                command = lambda: self.get_fp(self.gps_path))
        gps_path_b.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'W')
        
        # Create inputs for wind speed
        self.wind_speed = tk.IntVar(value = 10.0)
        wind_speed_l = tk.Label(trav_frame, text = 'Wind Speed (ms-1):', 
                                font = NORM_FONT)
        wind_speed_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        wind_speed_e = tk.Entry(trav_frame, font = NORM_FONT, width = 20, 
                                text = self.wind_speed)
        wind_speed_e.grid(row = 2, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        # Create inputs for wind error
        self.wind_error = tk.DoubleVar(value = 10)
        wind_error_l = tk.Label(trav_frame, text = 'Wind Error (%):', 
                                font = NORM_FONT)
        wind_error_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
        wind_error_e = tk.Entry(trav_frame, font = NORM_FONT, width = 20, 
                                text = self.wind_error)
        wind_error_e.grid(row = 3, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        # Create button to read in files
        read_b = ttk.Button(trav_frame, text="Read Traverse Data", 
                            command = self.read_trav_data)
        read_b.grid(row = 2, column = 2, padx = 5, pady = 5, columnspan=3, sticky = 'W')
        
#========================================================================================
#===================================Volcano controls=====================================
#========================================================================================
        
        # Create frame
        volc_frame = tk.LabelFrame(cont_frame, text='Volcano Settings', font=LARG_FONT)
        volc_frame.grid(row=1, column=0, padx=10, pady=6, columnspan = 2, sticky="NW")
        
        # Create inputs for volcano lonitude
        self.volc_lon = tk.DoubleVar(value = 0.0)
        volc_lon_l = tk.Label(volc_frame, text = 'Volcano Lon:', font = NORM_FONT)
        volc_lon_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        volc_lon_e = tk.Entry(volc_frame, font = NORM_FONT, width = 25, 
                                 text = self.volc_lon)
        volc_lon_e.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        # Create inputs for volcano latitude
        self.volc_lat = tk.DoubleVar(value = 0.0)
        volc_lat_l = tk.Label(volc_frame, text = 'Volcano Lat:', font = NORM_FONT)
        volc_lat_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        volc_lat_e = tk.Entry(volc_frame, font = NORM_FONT, width = 25, 
                                 text = self.volc_lat)
        volc_lat_e.grid(row = 2, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        # Create input for time difference
        self.time_diff = tk.DoubleVar(value = 0.0)
        time_diff_l = tk.Label(volc_frame,text='Time Difference (hours):',font=NORM_FONT)
        time_diff_l.grid(row = 3, column = 0, padx = 5, pady = 5, sticky = 'W')
        time_diff_e = tk.Entry(volc_frame, font = NORM_FONT, width = 25, 
                                 text = self.time_diff)
        time_diff_e.grid(row = 3, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        # Update parameters on choice
        def volc_update(event):
            self.volc_lon.set(volc_data[self.volc_name.get()][0])
            self.volc_lat.set(volc_data[self.volc_name.get()][1])
            self.time_diff.set(volc_data[self.volc_name.get()][2])
        
        # Create choice of imported volcanos
        self.volc_name = tk.StringVar(volc_frame, value = choice[0])
        self.volc_name_l = tk.Label(volc_frame, text = 'Select Volcano:', 
                                    font = NORM_FONT)
        self.volc_name_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.volc_name_c = ttk.OptionMenu(volc_frame, self.volc_name, *choice,
                                          command = volc_update)
        self.volc_name_c.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        
        
#========================================================================================
#===================================Analysis controls====================================
#========================================================================================
        
        # Create button to analyse
        begin_b = ttk.Button(cont_frame, text="Analyse", command = self.analyse)
        begin_b.grid(row = 2, column = 0, padx = 5, pady = 5)
        
        # Create a button to exit
        exit_b = ttk.Button(cont_frame, text = 'Exit', command = self.close_window)
        exit_b.grid(row = 2, column = 1, padx = 5, pady = 5)
        
        
        
        
        
        
#========================================================================================         
#========================================================================================
#======================================GUI Operations====================================
#======================================================================================== 
#======================================================================================== 

    # Report exceptions in a new window
    def report_callback_exception(self, *args):
        err = traceback.format_exception(*args)
        tkMessageBox.showerror('Exception', err)
        

    def handler(self):

        self.quit()        
        
        
        
        
        
        
        
#========================================================================================         
#========================================================================================
#=====================================Button Functions===================================
#======================================================================================== 
#======================================================================================== 

#========================================================================================
#=====================================Filepath Buttons===================================
#========================================================================================
   
    # Function to select spectra to analyse
    def get_fp(self, label):
            
        # Open dialouge to get files
        fpath = fd.askopenfilename()

        # Check not empty and update
        if fpath != '':
            
            # Save output to input 
            label.set(str(fpath))  
            
#========================================================================================
#======================================Print Outputs=====================================
#========================================================================================
     
    # Function to print text to the output box          
    def text_output(self, text, add_line = True):
        
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
#======================================Close window======================================
#========================================================================================
    
    def close_window(self):
        
        if len(common['fluxes']) != 0:
            # Calculate the average flux and uncertainty
            av_flux = 0
            av_err  = 0
            n = 0 
            
            for i in common['fluxes']:
                av_flux += i[0]
                av_err  += (i[1]/i[0])**2
                n += 1
                
            av_flux = av_flux / n
            av_err  = (av_err / n)**0.5 * av_flux
            
            # Save the fluxes in a text file for ease of access
            header = 'Results from calc_flux.py\n' + \
                     'NOTE errors are from SO2 fitting and wind speed only\n' + \
                     'Flux (tonnes/day),     Error (+/- t/d)'
            
            np.savetxt(common['out_folder'] + 'flux_results.txt', common['fluxes'],
                       header = header)
            
            with open(common['out_folder'] + 'flux_results.txt', 'a') as r:
                r.write('Average flux = ' + str(int(av_flux)) + ' +/- ' + \
                        str(int(av_err)) + ' (t/day)')
        
        self.quit()
        
#========================================================================================
#==================================Read Traverse data====================================
#========================================================================================
         
    def read_trav_data(self):
        
        # Reset loop counter
        common['loop'] = 0
        
        self.text_output('Reading traverse data...', add_line = False)
           
        # Read data file
        try:
            data = read_csv(str(self.ifit_path.get()))

            # Load required data into arrays
            txt_time    = data['Time'] 
            txt_so2_amt = data['so2 (ppm.m)']
            txt_so2_err = data['so2 error']
            
        except KeyError:
            self.text_output('ERROR: Wrong iFit file format')
            return
        
        # Create empty time array
        common['time'] = np.ones(len(txt_time))
        
        for n, t in enumerate(txt_time):
        
            # Convert the time to julian time
            try:
                timestamp = datetime.datetime.strptime(t, '%H:%M:%S.%f').time()
                
            except ValueError:
                timestamp = datetime.datetime.strptime(t, '%H:%M:%S').time()
            
            # Split into hours, minutes and seconds
            h = timestamp.hour
            m = timestamp.minute
            s = timestamp.second + timestamp.microsecond/1e6
            
            # Convert to julian time
            common['time'][n] = (h * 3600.0 + m * 60.0 + s) / 86400.0
            
         # Convert others from string to number arrays
        common['so2_amt'] = np.ones(len(txt_so2_amt))
        for n, i in enumerate(txt_so2_amt):
            common['so2_amt'][n] = float(i)
            
        common['so2_err'] = np.ones(len(txt_so2_err))
        for n, i in enumerate(txt_so2_err):
            common['so2_err'][n] = float(i)
        
        self.text_output('Done!', add_line = False) 
            
        # Read GPS data from file
        self.text_output('Reading GPS data...', add_line = False)

        try:
            gps_out = read_txt_gps(str(self.gps_path.get()))
            common['gps_time'], common['lat'], common['lon'] = gps_out
            
        except:
            self.text_output('ERROR: Wrong GPS file format')
            return
        
        self.text_output('Done!')        
        
        # Create output folder
        common['out_folder'] = 'Results/calc_flux/' + str(data['Date'][0]) + '/'
        if not os.path.exists(common['out_folder']):
            os.makedirs(common['out_folder'])
            
        # Create output .csv file if it doesn't exist
        common['out_fname'] = common['out_folder'] + 'flux_out.csv'
        
        with open(common['out_fname'], 'w') as writer:
            writer.write('calc_flux.py output file \n')
        
        # Update graph
        self.line0.set_xdata(common['time'])
        self.line0.set_ydata(common['so2_amt'])
        
        self.vline0.set_xdata(common['time'][0])
        self.vline1.set_xdata(common['time'][-1])
        
        self.fill.remove()
        lo_bound = np.subtract(common['so2_amt'], common['so2_err'])
        hi_bound = np.add(common['so2_amt'], common['so2_err'])
        self.fill = self.ax0.fill_between(common['time'], lo_bound, hi_bound, 
                                          color = 'lightgreen')
        
        # Reset sliders
        self.slide_hi.reset()
        self.slide_lo.reset()
        
        # Rescale the axes
        self.ax0.relim()
        self.ax0.autoscale_view() 
        
        # Apply changes
        self.fig.canvas.draw()
        
#========================================================================================
#=======================================Analysis=========================================
#========================================================================================
    
    def analyse(self):
        
        try:
            # Read in origional data set
            time       = common['time']
            so2_amt    = common['so2_amt']
            so2_err    = common['so2_err']
            gps_time   = common['gps_time']
            lon        = common['lon']
            lat        = common['lat']
            volc_lon   = self.volc_lon.get()
            volc_lat   = self.volc_lat.get()
            wind_speed = float(self.wind_speed.get())
            wind_err   = float(self.wind_error.get())
            
        except KeyError:
            self.text_output('Please import data first')
            return
        
        # Correct for time difference
        gps_time = np.subtract(gps_time, int(self.time_diff.get()) / 24)
        
        # Check for change of day
        for n, t in enumerate(gps_time):
            if t > 1:
                gps_time[n] = t - 1
            
        # Get indices of selected data
        time_grid, idx0, idx1 = extract_window(time, self.t0, self.t1)
        idx1 += 1
        
        # Catch error if bounds are the wrong way around
        if idx0 > idx1:
            self.text_output('Bounds wrong way around, please try again')
            return
        
        # Cut to extract the selected traverse
        time    = time[idx0:idx1]
        so2_amt = so2_amt[idx0:idx1]
        so2_err = so2_err[idx0:idx1]

#========================================================================================
#============================Define peak of SO2 concentration============================
#========================================================================================

        peak_idx = cog(so2_amt)   

#========================================================================================
#========================Interpolate GPS data onto SO2 time grid=========================
#========================================================================================

        # Interpolate the GPS data onto the so2 time grid
        modlon_old = griddata(gps_time,lon,common['time'])
        modlat_old = griddata(gps_time,lat,common['time'])
        modlat = griddata(gps_time,lat,time)
        modlon = griddata(gps_time,lon,time)

#========================================================================================
#==============================Perform geometric correction==============================
#========================================================================================

        # Define wind vector as the vector from the volcano to the CoM
        wind_vector = haversine(volc_lon, volc_lat, modlon[peak_idx], modlat[peak_idx])

        # Find the distance-bearing vectors between each reading. Bearing is radians 
        #  anticlockwise from East from 0 to 2pi
        displacement, bearing, dir_corr = gps_vector(modlon, modlat, wind_vector[1])
        
        # Correction factor to account for the non-orthogonality of the plume
        orthog_correction = np.ones(len(bearing))
        for n, i in enumerate(bearing):
            correction = abs(np.cos((np.pi / 2) - (wind_vector[1] - i)))
            orthog_correction[n] = correction
            
        # Correct for if going the wrong way
        orthog_correction = np.multiply(orthog_correction, dir_corr)    

        # Apply the correction
        corr_dists = np.multiply(displacement, orthog_correction)    
        
        # Calculate cumulative distance
        cuml_dist = np.zeros(len(corr_dists) + 1)
        for n, i in enumerate(corr_dists):
            cuml_dist[n+1] = cuml_dist[n] + i

#========================================================================================
#=====================================Calculate flux=====================================
#========================================================================================

        # Convert so2 amounts from ppm.m to molecules.cm-2
        so2_amt_molec = np.multiply(so2_amt, 2.463e15)

        # Convert so2 amounts from molecules.cm-2 to molecules.m-2
        molec_m2 = so2_amt_molec[1:] * 1.0e4
    
        # Convert to molecules.m-1 and sum over the traverse
        molec_m = np.sum(np.multiply(molec_m2, corr_dists))

        # Multiply by windspeed
        molec_per_s = molec_m * wind_speed
    
        # Convert to moles.s-1
        moles_per_s = molec_per_s / 6.022e23
    
        # Convert to kg.s-1. Molar mass of so2 is 64.066g/mole
        kg_per_s = moles_per_s * 0.064066
    
        # Convert to tonnes per day, make sure it's positive
        try:
            flux = abs(int((kg_per_s * 86400.0) / (1000.0)))
            
        except ValueError:
            self.text_output('Error in flux calulation.\n' +\
                             'Please ensure volcano parameters are correct')
            return
        
#========================================================================================
#==================================Calculate flux error==================================
#========================================================================================
    
        # Calculate error in so2 amount
        total_err = np.sum(so2_err)
        total_so2 = np.sum(so2_amt)
        
        delta_A = total_err / total_so2
        
        # Combine with wind speed error
        delta_F = ( (delta_A)**2 + (wind_err/100)**2 )**0.5
        
        flux_err = int(flux * delta_F)
        
        self.text_output('Flux = '+str(flux)+' (+/- ' + str(flux_err)+') tonnes/day')

#========================================================================================
#=============================Plot selected SO2 and GPS data=============================
#========================================================================================
           
        # Create dictionary to pass to plotting function
        d = {}
        d['time']              = time
        d['so2_amt']           = so2_amt
        d['so2_err']           = so2_err
        d['peak_idx']          = peak_idx
        d['modlon']            = modlon
        d['modlat']            = modlat
        d['modlon_old']        = modlon_old
        d['modlat_old']        = modlat_old
        d['volc_lon']          = volc_lon
        d['volc_lat']          = volc_lat
        d['flux']              = flux
        d['flux_err']          = flux_err
        d['wind_speed']        = wind_speed
        d['bearing']           = bearing
        d['displacement']      = displacement
        d['orthog_correction'] = orthog_correction
        d['cuml_dist']         = cuml_dist
        

        make_graph(d)
        
        
#========================================================================================
#=======================================Make Graph=======================================
#========================================================================================
            
            
# Define function to open another window to show output graphs
def make_graph(d):
    popup = tk.Tk()
    tk.Tk.wm_title(popup, 'Results')
    
    def save(fig, d):
        
        # Save the figure
        fig.savefig(common['out_folder'] + 'traverse_' + str(common['loop']) + '.png')
        
        # Append flux and error to list
        common['fluxes'].append([d['flux'], d['flux_err']]) 

        # Make string of cog position
        cog_pos = str(d['modlon'][d['peak_idx']])+','+str(d['modlat'][d['peak_idx']])
        volc_loc = str(d['volc_lon']) + '/' + str(d['volc_lat'])

        with open(common['out_fname'], 'a') as writer:
            # Write header information to the file
            writer.write('Peak number,' + str(common['loop'])                 + '\n' \
                         'Flux (Tonnes/day),' + str(d['flux'])                + '\n' \
                         'Windspeed (m/s),' + str(d['wind_speed'])            + '\n' \
                         'Volcano Long/Lat,' + volc_loc                       + '\n' \
                         'CoM position,' + cog_pos                            + '\n' \
                         'Plume Start Time,' + str(d['time'][0])              + '\n' \
                         'Plume Centre Time,' + str(d['time'][d['peak_idx']]) + '\n' \
                         'Plume Stop Time,' + str(d['time'][-1])              + '\n' )
             
            # Write column headers            
            writer.write('Julian Time (GMT),SO2 Amount,Error,Longitude,Latitude,' + \
                         'Bearing,Distance,Correction Factor,Cumulative Distance\n')
            
            # Write information to the file
            for i in range(1,len(d['modlon'])):
                writer.write(str(d['time'][i])                + ',' + \
                             str(d['so2_amt'][i])             + ',' + \
                             str(d['so2_err'][i])             + ',' + \
                             str(d['modlon'][i])              + ',' + \
                             str(d['modlat'][i])              + ',' + \
                             str(d['bearing'][i-1])           + ',' + \
                             str(d['displacement'][i-1])      + ',' + \
                             str(d['orthog_correction'][i-1]) + ',' + \
                             str(d['cuml_dist'][i])           + '\n')

            writer.write('\n')
        
        # Add to the loop counter
        common['loop'] += 1
        
        # Close the window
        popup.destroy()
    
    
    
    # Create figure to hold the graphs
    fig = plt.figure(figsize = (10,6))
    
    # Create plot grid
    gs  = gridspec.GridSpec(2, 3, width_ratios = [1,0.95,0.05])
    ax1 = plt.subplot(gs[0,0])  
    ax2 = plt.subplot(gs[0,1])
    ax3 = plt.subplot(gs[1,0])
    ax4 = plt.subplot(gs[1,1])   

    # Plot the full SO2 amount time series
    ax1.plot(d['time'], d['so2_amt'],  'g', lw = 1.5)
    ax1.fill_between(d['time'], np.subtract(d['so2_amt'], d['so2_err']), 
                     np.add(d['so2_amt'], d['so2_err']), color = 'lightgreen')
    ax1.plot(d['time'][d['peak_idx']], d['so2_amt'][d['peak_idx']], 'ko')
    ax1.set_xlabel('Julian Time (Fraction of Day)', fontsize = 10)
    ax1.set_ylabel('SO2 column amount (ppm.m)', fontsize = 10)

    # Plot the selected traverse points and SO2 amounts wrt the volcano
    sc = ax2.scatter(d['modlon'], d['modlat'], c = d['so2_amt'], s=50,
                     lw=0.0,alpha = 0.5)
    cbax = plt.subplot(gs[0,2])
    cb = Colorbar(ax = cbax, mappable = sc)
    cb.set_label('SO2 column amount (ppm.m)', fontsize = 10)
    
    ax2.scatter(d['volc_lon'], d['volc_lat'], c='darkorange', s=180, label = 'Volcano')
    
    ax2.set_xlabel('Longitude (deg)', fontsize = 10)
    ax2.set_ylabel('Latitude (deg)', fontsize = 10)

    # Plot the full SO2 amount time series overlaid with the selected section
    ax3.plot(common['time'], common['so2_amt'], 'g-')
    ax3.fill_between(common['time'], np.subtract(common['so2_amt'], common['so2_err']),
                     np.add(common['so2_amt'], common['so2_err']), color = 'lightgreen')
    ax3.plot(d['time'], d['so2_amt'], 'r-')
    ax3.fill_between(d['time'], np.subtract(d['so2_amt'], d['so2_err']),
                     np.add(d['so2_amt'], d['so2_err']), 
                     color = 'coral')
    ax3.set_xlabel('Julian Time (Fraction of Day)', fontsize = 10)
    ax3.set_ylabel('SO2 column amount (ppm.m)', fontsize = 10)
    
    # Plot the full GPS track and overlay the selected window wrt the volcano location
    ax4.plot(d['modlon_old'], d['modlat_old'], 'g-', label = 'Full Track')
    ax4.plot(d['modlon'], d['modlat'], 'r', lw = 1.5, label = 'Traverse')
    ax4.scatter(d['modlon'][0], d['modlat'][0], c='g', s=100, label = 'Start')
    ax4.scatter(d['modlon'][-1], d['modlat'][-1], c='r', s=100, label = 'Stop')
    ax4.scatter(d['volc_lon'], d['volc_lat'], c='darkorange', s=200, label = 'Volcano')
    ax4.scatter(d['modlon'][d['peak_idx']], d['modlat'][d['peak_idx']], c='k', s=100,
                label = 'Plume CoM')
    ax4.set_xlabel('Longitude (deg)', fontsize = 10)
    ax4.set_ylabel('Latitude (deg)', fontsize = 10)
    ax4.legend(loc='center right', bbox_to_anchor=(1.5, 0.5), fontsize = 10)

    # Display the plots
    plt.tight_layout()
    
    # Create the canvas to hold the graph in the GUI
    canvas = FigureCanvasTkAgg(fig, popup)
    canvas.show()
    canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady = 10, columnspan = 2)
    
    # Create a button to save 
    save_b = ttk.Button(popup, text="Save", command = lambda: save(fig, d))
    save_b.grid(row = 1, column = 0, padx = 5, pady = 5)
    
    # Create a button to exit
    cancel_b = ttk.Button(popup, text = 'Cancel', command = popup.destroy)
    cancel_b.grid(row = 1, column = 1, padx = 5, pady = 5)

# Run the App!     
mygui().mainloop()