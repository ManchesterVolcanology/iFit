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
from tkinter import filedialog as fd
from math import degrees
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
from matplotlib.widgets import Slider
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from ifit_lib.read_csv import read_csv
import folium

from ifit_lib.find_nearest import extract_window
from ifit_lib.read_gps import read_gps, gps_vector, haversine
from ifit_lib.center_of_grav import cog
from ifit_lib.julian_time import hms_to_julian, julian_to_hms
from ifit_lib.build_gui import make_input

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

        # Close program on closure of window
        self.protocol("WM_DELETE_WINDOW",self.handler)

        # Button Style
        ttk.Style().configure('TButton', width = 20, height = 20, relief = 'flat')

        # Add a title and icon
        tk.Tk.wm_title(self, 'calc_flux 2.1')
        try:
            tk.Tk.iconbitmap(self, default = 'data_bases/calc_flux_icon.ico')
        except tk.TclError:
            pass

        # Create a left and right master frame
        left = ttk.Frame(self)
        left.grid(row=0, column=0)
        right = ttk.Frame(self)
        right.grid(row=0, column=1)

        # Create control Frame
        cont_frame = ttk.Frame(left)
        cont_frame.grid(row=0, column=0, padx=10, pady=10, rowspan=2)

        # Create frame to hold graphs
        graph_frame = ttk.Frame(right, relief = 'groove')
        graph_frame.grid(row=1, column=2, padx=10, pady=10, rowspan=2)
        graph_frame.columnconfigure(index=0, weight=1)
        graph_frame.rowconfigure(index = 0, weight = 1)

        # Create frame to hold text output
        text_frame = ttk.Frame(left, relief = 'groove')
        text_frame.grid(row=2, column=0, padx=10, pady=10, sticky = 'W')

        # Create flag to controll whether or not to save outputs on loading new data
        self.save_flag = False

        # Create global common dictionary to hold settings
        global common
        common = {}

        # Create loop counter to hold traverse number
        common['loop'] = 1

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
#=================================== Create text output =================================
#========================================================================================

        # Build text box
        self.text_box = tkst.ScrolledText(text_frame, width = 50, height = 10)
        self.text_box.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W',
                           columnspan = 4)
        self.text_box.insert('1.0', 'Welcome to calc_flux! Written by Ben Esse\n\n')

#========================================================================================
#==================================== Build Graph =======================================
#========================================================================================

        # Create figure to hold the graphs
        plt.rcParams.update({'font.size': 8} )
        self.fig = plt.figure(figsize = (7,5))

        # Create the canvas to hold the graph in the GUI
        self.canvas = FigureCanvasTkAgg(self.fig, graph_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=0, column=0, padx=10,
                                         pady = 10, sticky = 'NW')

        # Create plot axes
        self.ax0 = self.fig.add_subplot(111)
        self.fig.subplots_adjust(bottom=0.3)

        self.line0, = self.ax0.plot(0, 0, lw=1.5)
        self.vline0 = self.ax0.axvline(0, color = 'r')
        self.vline1 = self.ax0.axvline(1, color = 'r')
        self.fill   = self.ax0.fill_between([0], [0], [1])
        self.ax0.set_xlabel('Time (From Spectra)', fontsize=10)
        self.ax0.set_ylabel(r'SO$_2$ column amount (ppm.m)', fontsize=10)

        # Get origional limits
        common['lims'] = self.ax0.get_xlim()

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

        # Add buttons to zoom and reset
        button_frame = tk.Frame(right)
        button_frame.grid(row=0, column=2, padx=10, pady=10)
        zoom_b = ttk.Button(button_frame, text = 'Zoom', command = self.zoom)
        zoom_b.grid(row=0, column=0, sticky='E', padx=10)
        reset_b = ttk.Button(button_frame, text = 'Reset', command = self.reset)
        reset_b.grid(row=0, column=2, sticky='E', padx=10)

#========================================================================================
#================================= Traverse controls ====================================
#========================================================================================

        # Create frame
        trav_frame = tk.LabelFrame(cont_frame, text='Traverse Settings', font=LARG_FONT)
        trav_frame.grid(row=0, column=0, padx=10, pady=6, columnspan = 2, sticky="NW")

        # Create inputs for iFit output file
        self.ifit_path = tk.StringVar(value = 'No file selected')
        tk.Label(trav_frame, text = 'SO2 data:', font = NORM_FONT
                 ).grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        tk.Entry(trav_frame, font = NORM_FONT, width = 30, text = self.ifit_path
                 ).grid(row = 0, column = 1, padx = 5, pady = 5, sticky = 'W',
                        columnspan = 2)
        ttk.Button(trav_frame, text="Select SO2 File",
                   command = lambda: self.get_fp(self.ifit_path)
                   ).grid(row = 1, column=2, padx=5, pady=5, sticky='W')

        # Control whether or not to remove poor fits
        self.de_spike = tk.BooleanVar(trav_frame, value = True)
        make_input(frame = trav_frame,
                   text = 'Remove Bad\nSpectra?',
                   row = 1, column = 0,
                   var = self.de_spike,
                   input_type = 'Checkbutton')

        # Create inputs for gps file
        self.gps_path = tk.StringVar(value = 'No file selected')
        gps_path_l = tk.Label(trav_frame, text = 'GPS data:', font = NORM_FONT)
        gps_path_l.grid(row = 2, column = 0, padx = 5, pady = 5, sticky = 'W')
        gps_path_e = tk.Entry(trav_frame, font = NORM_FONT, width = 30,
                              text = self.gps_path)
        gps_path_e.grid(row = 2, column = 1, padx = 5, pady = 5, sticky = 'W',
                        columnspan = 2)
        gps_path_b = ttk.Button(trav_frame, text="Select GPS file",
                                command = lambda: self.get_fp(self.gps_path))
        gps_path_b.grid(row=3, column=2, padx=5, pady=5, sticky='W')

        # Control type of GPS data
        gps_types = ['text', 'text', 'NMEA']
        self.gps_dtype = tk.StringVar(trav_frame, value = 'text')
        make_input(frame = trav_frame,
                   text = 'GPS Data\nType',
                   row = 3, column = 0,
                   var = self.gps_dtype,
                   input_type = 'OptionMenu',
                   options = gps_types,
                   width = 6)

#========================================================================================
#================================== Volcano controls ====================================
#========================================================================================

        # Create frame
        volc_frame = tk.LabelFrame(cont_frame, text='Volcano Settings', font=LARG_FONT)
        volc_frame.grid(row=1, column=0, padx=10, pady=10, sticky="NW", columnspan=2)

        # Update parameters on choice
        def volc_update(event):
            self.volc_lon.set(volc_data[self.volc_name.get()][0])
            self.volc_lat.set(volc_data[self.volc_name.get()][1])
            self.time_diff.set(volc_data[self.volc_name.get()][2])

        # Create choice of imported volcanos
        self.volc_name = tk.StringVar(volc_frame, value = choice[0])
        make_input(frame = volc_frame,
                   text = 'Volcano:',
                   row = 0, column = 0,
                   var = self.volc_name,
                   input_type = 'OptionMenu',
                   options = choice,
                   width = 13,
                   command = volc_update,
                   columnspan = 3,
                   sticky = 'W')

        # Create inputs for volcano latitude
        self.volc_lat = tk.DoubleVar(value = 0.0)
        make_input(frame = volc_frame,
                   text = 'Volcano\nLatitude:',
                   row = 1, column = 0,
                   var = self.volc_lat,
                   input_type = 'Entry',
                   width = 15,
                   sticky = 'W')

        # Create inputs for volcano lonitude
        self.volc_lon = tk.DoubleVar(value = 0.0)
        make_input(frame = volc_frame,
                   text = 'Volcano\nLongitude:',
                   row = 2, column = 0,
                   var = self.volc_lon,
                   input_type = 'Entry',
                   width = 15,
                   sticky = 'W')

        # Create input for time difference
        self.time_diff = tk.DoubleVar(value = 0.0)
        make_input(frame = volc_frame,
                   text = 'Time\nDifference:',
                   row = 3, column = 0,
                   var = self.time_diff,
                   input_type = 'Entry',
                   width = 15,
                   sticky = 'W')

        # Create inputs for wind speed and units
        self.wind_speed = tk.DoubleVar(value = 0)
        make_input(frame = volc_frame,
                   text = 'Wind\nSpeed:',
                   row = 1, column = 2,
                   var = self.wind_speed,
                   input_type = 'Entry',
                   width = 10,
                   sticky = 'W')
        self.wind_unit = tk.StringVar(value = 'm/s')
        wind_speed_u = ttk.OptionMenu(volc_frame, self.wind_unit,
                                      *['m/s','m/s','knots'])
        wind_speed_u.config(width = 6)
        wind_speed_u.grid(row = 1, column = 4, padx = 5, pady = 5, sticky = 'EW')

        # Create inputs for wind error and units
        self.wind_error = tk.DoubleVar(value = 0)
        make_input(frame = volc_frame,
                   text = 'Wind\nError:',
                   row = 2, column = 2,
                   var = self.wind_error,
                   input_type = 'Entry',
                   width = 10,
                   sticky = 'W')
        self.error_unit = tk.StringVar(value = '%')
        wind_error_u = ttk.OptionMenu(volc_frame, self.error_unit,
                                      *['%','%','abs'])
        wind_error_u.config(width = 6)
        wind_error_u.grid(row = 2, column = 4, padx = 5, pady = 5, sticky = 'EW')



#========================================================================================
#================================== Analysis controls ===================================
#========================================================================================

        # Create button to read in files
        read_b = ttk.Button(cont_frame, text="   Load Traverse",
                            command = self.read_trav_data)
        read_b.grid(row = 2, column = 0, padx = 5, pady = 5)

        # Create button to analyse
        begin_b = ttk.Button(cont_frame, text = "Calculate Flux",
                             command = self.analyse)
        begin_b.grid(row = 2, column = 1, padx = 5, pady = 5)






#========================================================================================
#========================================================================================
#===================================== GUI Operations ===================================
#========================================================================================
#========================================================================================

    # Report exceptions in a new window
    def report_callback_exception(self, *args):
        err = traceback.format_exception(*args)
        tkMessageBox.showerror('Exception', err)

    def handler(self):

        # Close
        self.quit()







#========================================================================================
#========================================================================================
#==================================== Button Functions ==================================
#========================================================================================
#========================================================================================

#========================================================================================
#=================================== Zoom/reset Buttons =================================
#========================================================================================

    def zoom(self):

        # Get slider values
        pos = self.slide_lo.val, self.slide_hi.val

        # Get graph limits
        lims = self.ax0.get_xlim()
        lim_range = lims[1] - lims[0]

        # Update axis limits
        new_lims = [lims[0] + lim_range * pos[0],
                    lims[0] + lim_range * pos[1]]
        self.ax0.set_xlim(new_lims)
        self.rescale()

        # Reset sliders
        self.slide_hi.reset()
        self.slide_lo.reset()

        # Apply changes
        self.fig.canvas.draw()

    def reset(self):

        try:
            # Reset limits
            self.ax0.set_xlim(common['lims'])
            self.rescale()

            # Apply changes
            self.fig.canvas.draw()

        except KeyError:
            self.text_output('Please import data first')

    def rescale(self):

        try:
            # Get graph x limits
            lims = self.ax0.get_xlim()

            # Get current data
            idx = np.where(np.logical_and(common['time'] > lims[0],
                                          common['time'] < lims[1]))
            y = common['so2_amt'][idx]

            # Find 5% of the range of data
            y_pad = (max(y) - min(y)) * 0.05

            # Set limits
            y_lims = [min(y) - y_pad, max(y) + y_pad]

            # Set axis y limits
            self.ax0.set_ylim(y_lims)

            # Apply changes
            self.fig.canvas.draw()

        except KeyError:
            self.text_output('Please import data first')

#========================================================================================
#==================================== Filepath Buttons ==================================
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
#===================================== Print Outputs ====================================
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
#================================= Read Traverse data ===================================
#========================================================================================

    def read_trav_data(self):

        # Reset the graph
        self.reset()

        # Reset loop counter
        common['loop'] = 1

        # Reset flux amount array
        common['times']     = []
        common['locations'] = []
        common['azimuths']  = []
        common['w_speeds']  = []
        common['fluxes']    = []
        common['flux_errs'] = []

        self.text_output('Reading traverse data...', add_line = False)

        # Read data file
        try:
            data = read_csv(str(self.ifit_path.get()))

            # Load required data into arrays
            txt_time    = data['Time']
            txt_so2_amt = data['so2_amt']/2.463e15
            txt_so2_err = data['so2_amt_e']/2.463e15
            fit_flag    = data['Fit Quality']

        except KeyError:
            self.text_output('ERROR: Wrong iFit file format')
            return

        # Convert timestamps to decimal hours, first trying with us then without
        try:
            common['time'] = hms_to_julian(txt_time, str_format = '%H:%M:%S',
                                           out_format = 'decimal hours')
        except ValueError:
            common['time'] = hms_to_julian(txt_time, str_format = '%H:%M:%S.%f',
                                           out_format = 'decimal hours')

         # Convert others from string to number arrays
        common['so2_amt'] = np.asarray(txt_so2_amt, dtype = float)
        common['so2_err'] = np.asarray(txt_so2_err, dtype = float)

        self.text_output('Done!', add_line = False)

        # Remove bad datapoints
        if self.de_spike.get():

            # Get indicies of bad spectra
            spec_qual = np.zeros(len(txt_so2_amt))

            for n, q in enumerate(fit_flag):
                if q == 'Good':
                    spec_qual[n] = 1

            idx = np.where(spec_qual==1)

            # Remove bad data
            common['time'] = common['time'][idx]
            common['so2_amt'] = common['so2_amt'][idx]
            common['so2_err'] = common['so2_err'][idx]

        # Read GPS data from file
        self.text_output('Reading GPS data...', add_line = False)

        try:
            gps_out = read_gps(self.gps_path.get(), self.gps_dtype.get())
            common['gps_time'], common['lat'], common['lon'], common['alt'] = gps_out

        except:
            self.text_output('ERROR: Wrong GPS file format')
            return

        self.text_output('Done!')

        # Get date of the data
        common['analysis_date'] = str(data['Date'][0])

        # Create output folder
        common['out_folder'] = 'Results/calc_flux/' + common['analysis_date'] + '/'
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
                                          color = 'lightblue')

        # Reset sliders
        self.slide_hi.reset()
        self.slide_lo.reset()

        # Find 5% of the range of data
        x = common['time']
        y = common['so2_amt']
        x_pad = (max(x) - min(x))*0.05
        y_pad = (max(y) - min(y))*0.05

        # Set limits
        x_lims = [min(x) - x_pad, max(x) + x_pad]
        y_lims = [min(y) - y_pad, max(y) + y_pad]

        # Rescale the axes
        self.ax0.set_xlim(x_lims)
        self.ax0.set_ylim(y_lims)

        # Get origional limits
        common['lims'] = self.ax0.get_xlim()

        # Apply changes
        self.fig.canvas.draw()

        # Turn on save flag
        self.save_flag = True

#========================================================================================
#====================================== Analysis ========================================
#========================================================================================

    def analyse(self):

        try:
            # Read in origional data set
            time       = common['time']
            so2_amt    = common['so2_amt']
            so2_err    = common['so2_err']
            gps_time   = common['gps_time']
            lon        = np.array(common['lon'])
            lat        = np.array(common['lat'])
            volc_lon   = self.volc_lon.get()
            volc_lat   = self.volc_lat.get()
            wind_speed = self.wind_speed.get()
            wind_err   = self.wind_error.get()

        except KeyError:
            self.text_output('Please import data first')
            return

        # If wind error is absolute, turn into %
        if self.error_unit.get() == 'abs':
            wind_err = wind_err / wind_speed * 100

        # Get indices of selected data
        time_grid, idx0, idx1 = extract_window(time, self.t0, self.t1)
        idx1 += 1

        # Correct for time difference
        gps_time = np.add(gps_time, int(self.time_diff.get()))

        # Check for change of day
        for n, t in enumerate(time):
            if t > 24:
                time[n] = t - 24

        # Catch error if bounds are the wrong way around
        if idx0 > idx1:
            self.text_output('Bounds wrong way around, please try again')
            return

        # Cut to extract the selected traverse
        time    = time[idx0:idx1]
        so2_amt = so2_amt[idx0:idx1]
        so2_err = so2_err[idx0:idx1]

#========================================================================================
#=========================== Define peak of SO2 concentration ===========================
#========================================================================================

        # Get the centre of the plume
        peak_idx = cog(so2_amt)

#========================================================================================
#======================= Interpolate GPS data onto SO2 time grid ========================
#========================================================================================

        # Interpolate the GPS data onto the so2 time grid
        modlon_old = griddata(gps_time, lon, common['time'])
        modlat_old = griddata(gps_time, lat, common['time'])
        modlat     = griddata(gps_time, lat, time)
        modlon     = griddata(gps_time, lon, time)

#========================================================================================
#============================= Perform geometric correction =============================
#========================================================================================

        # Define wind vector as the vector from the volcano to the CoM
        wind_dist, wind_az = haversine(volc_lon, volc_lat, modlon[peak_idx], modlat[peak_idx])

        # Add plume peak location and azimuths to arrays
        common['locations'].append([modlat[peak_idx], modlon[peak_idx]])
        common['azimuths'].append(degrees(wind_az))

        # Find the distance-bearing vectors between each reading. Bearing is radians
        #  anticlockwise from East from 0 to 2pi
        displacement, bearing, dir_corr = gps_vector(modlon, modlat, wind_az)

        # Correction factor to account for the non-orthogonality of the plume
        orthog_correction = np.ones(len(bearing))
        for n, i in enumerate(bearing):
            correction = abs(np.cos((np.pi / 2) - (wind_az - i)))
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
#==================================== Calculate flux ====================================
#========================================================================================

        # If wind speed is in knotts, convert to m/s
        if self.wind_unit.get() == 'knots':
            wind_speed = wind_speed * 0.5144444444

        # Add wind speed to array
        common['w_speeds'].append([wind_speed, wind_speed * (wind_err / 100) ])

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
#================================= Calculate flux error =================================
#========================================================================================

        # Get error in column densities
        xsec_err = np.sum(np.power(so2_err, 2)) ** 0.5
        delta_A = xsec_err / np.sum(so2_amt)

        # Combine with wind speed error
        delta_F = ( (delta_A)**2 + (wind_err/100)**2 )**0.5

        flux_err = int(flux * delta_F)

        self.text_output('Peak ' + str(common['loop']) + \
                         ': Flux = '+str(flux)+' (+/- ' + str(flux_err)+') tonnes/day')

#========================================================================================
#============================ Plot selected SO2 and GPS data ============================
#========================================================================================

        # Remove NaN values
        modlon_old = modlon_old[~np.isnan(modlon_old)]
        modlat_old = modlat_old[~np.isnan(modlat_old)]

        # Create dictionary to pass to plotting function
        d = {}
        d['time']              = time
        d['time_diff']         = self.time_diff.get()
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
        d['wind_az']           = wind_az
        d['bearing']           = bearing
        d['displacement']      = displacement
        d['orthog_correction'] = orthog_correction
        d['cuml_dist']         = cuml_dist


        make_graph(d)


#========================================================================================
#====================================== Make Graph ======================================
#========================================================================================


# Define function to open another window to show output graphs
def make_graph(d):
    popup = tk.Tk()
    tk.Tk.wm_title(popup, 'Results')

    def save(fig, d):

        # Save the figure
        fig.savefig(common['out_folder'] + 'traverse_' + str(common['loop']) + '.png')

        # Append flux and error to list
        common['times'].append(julian_to_hms(d['time'][d['peak_idx']] - d['time_diff']))
        common['fluxes'].append(d['flux'])
        common['flux_errs'].append(d['flux_err'])

        # Make string of cog position
        cog_pos = str(d['modlon'][d['peak_idx']])+','+str(d['modlat'][d['peak_idx']])
        volc_loc = str(d['volc_lon']) + ',' + str(d['volc_lat'])

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
            writer.write('Julian Time,SO2 Amount,Error,Longitude,Latitude,Bearing,' + \
                         'Distance,Correction Factor,Cumulative Distance\n')

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


        # Produce map with gps track
        trav_map = folium.Map(location = [d['volc_lat'], d['volc_lon']],
                              zoom_start = 13,
                              tiles = 'Stamen Terrain')

        folium.Circle(radius = 100,
                      location = [d['volc_lat'], d['volc_lon']],
                      color = 'orange',
                      fill = True,
                      tooltip = 'Source').add_to(trav_map)

        folium.PolyLine(locations = np.column_stack((d['modlat_old'], d['modlon_old'])),
                        tooltip = 'GPS Track').add_to(trav_map)

        folium.PolyLine(locations = np.column_stack((d['modlat'], d['modlon'])),
                        tooltip = 'Traverse',
                        color = 'orange').add_to(trav_map)

        trav_map.save(common['out_folder'] + 'Map_' + str(common['loop']) + '.html')

        # Update text results file
        with open(common['out_folder'] + 'flux_results.txt', 'w') as w:

            # Write header lines
            w.write('Results from calc_flux.py for ' + common['analysis_date'] + '\n' + \
                    'NOTE errors are from SO2 fitting and wind speed only\n\n' + \
                    'Time (UTC)    Centre Location (lat/lon)   Plume Azimuth    ' + \
                    'Wind speed (m/s)    Flux (t/day)\n')

            # Write each time, flux and error
            for n, t in enumerate(common['times']):

                # Create central location string
                c_lat = "{:.4f}".format(common['locations'][n][0])
                c_lon = "{:.4f}".format(common['locations'][n][1])
                cent_loc = '{0: <28}'.format(c_lat + ' / ' + c_lon)

                # Create wind speed string
                w_spd = "{:.2f}".format(common['w_speeds'][n][0])
                w_err = "{:.1f}".format(common['w_speeds'][n][1])

                w.write(str(t)[:12] + '  ' +  \
                        cent_loc + \
                        '{0: <17}'.format("{:.1f}".format(common['azimuths'][n])) + \
                        '{0: <20}'.format(w_spd + ' (+/- ' + w_err + ')') + \
                        str(common['fluxes'][n]) + '(+/- ' + \
                        str(common['flux_errs'][n]) + ')\n')

            # Calculate the average flux and uncertainty and write
            av_flux = np.average(common['fluxes'])
            av_err  = np.average(np.power(common['flux_errs'], 2)**0.5)
            w.write('\nAverage flux = ' + str(int(av_flux)) + ' +/- ' + \
                    str(int(av_err)) + ' (t/day)')

        # Add to the loop counter
        common['loop'] += 1

        # Close the window
        popup.destroy()



    # Create figure to hold the graphs
    fig = plt.figure(figsize = (10,6))

    # Create plot grid
    gs  = gridspec.GridSpec(2, 2)
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[0,1])
    ax3 = plt.subplot(gs[1,0])
    ax4 = plt.subplot(gs[1,1])

    # Plot the full SO2 amount time series
    ax1.fill_between(d['time'], np.subtract(d['so2_amt'], d['so2_err']),
                     np.add(d['so2_amt'], d['so2_err']), color = '#fcc006')
    ax1.plot(d['time'], d['so2_amt'], 'C1', lw = 1.5)
    ax1.plot(d['time'][d['peak_idx']], d['so2_amt'][d['peak_idx']], 'ko')
    ax1.set_xlabel('Julian Time (Decimal hours)', fontsize = 10)
    ax1.set_ylabel(r'SO$_2$ SCD (ppm.m)', fontsize = 10)

    # Plot the selected traverse points and SO2 amounts wrt the volcano
    sc = ax2.scatter(d['modlon'], d['modlat'], c = d['so2_amt'], s = 50, lw = 0.0,
                     alpha = 0.5)

    cax = make_axes_locatable(ax2).append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(sc, cax = cax)
    cb.set_label(label = r'SO$_2$ SCD (ppm.m)', fontsize = 10)

    ax2.scatter(d['volc_lon'], d['volc_lat'], c='darkorange', s=180, label = 'Volcano')

    ax2.set_xlabel('Longitude (deg)', fontsize = 10)
    ax2.set_ylabel('Latitude (deg)', fontsize = 10)

    # Plot the full SO2 amount time series overlaid with the selected section
    ax3.plot(common['time'], common['so2_amt'])
    ax3.fill_between(common['time'], np.subtract(common['so2_amt'], common['so2_err']),
                     np.add(common['so2_amt'], common['so2_err']), color = 'lightblue')
    ax3.plot(d['time'], d['so2_amt'])
    ax3.fill_between(d['time'], np.subtract(d['so2_amt'], d['so2_err']),
                     np.add(d['so2_amt'], d['so2_err']),
                     color = '#fcc006')
    ax3.set_xlabel('Julian Time (Decimal hours)', fontsize = 10)
    ax3.set_ylabel(r'SO$_2$ SCD (ppm.m)', fontsize = 10)

    # Plot the full GPS track and overlay the selected window wrt the volcano location
    ax4.plot(d['modlon_old'], d['modlat_old'], label = 'Full Track')
    ax4.plot(d['modlon'], d['modlat'], lw = 2.0, label = 'Traverse')
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
    canvas.get_tk_widget().grid(row=1, column=0, padx=10, pady = 10, columnspan = 2)

    # Add matplotlib toolbar above the plot canvas
    toolbar_frame = tk.Frame(popup, bg = 'black')
    toolbar_frame.grid(row=0,column=0, sticky = 'NW', padx = 6, columnspan=2)
    toolbar = NavigationToolbar2TkAgg(canvas, toolbar_frame)
    toolbar.update()

    # Create a button to save
    save_b = ttk.Button(popup, text="Save", command = lambda: save(fig, d))
    save_b.grid(row = 2, column = 0, padx = 5, pady = 5)

    # Create a button to exit
    cancel_b = ttk.Button(popup, text = 'Cancel', command = popup.destroy)
    cancel_b.grid(row = 2, column = 1, padx = 5, pady = 5)

# Run the App!
mygui().mainloop()