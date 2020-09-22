# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 17:11:08 2020

@author: mqbpwbe2
"""

import logging
import traceback
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import tkinter as tk
from tkinter import ttk
import tkinter.messagebox as tkMessageBox
import tkinter.scrolledtext as ScrolledText
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from ifitgui.build_gui import make_input
from ifitgui.gui_functions import select_files
from ifit.haversine import haversine

# Define some fonts to use in the program
NORM_FONT = ('TkFixedFont', 10)
LARGE_FONT = ('TkFixedFont', 12, 'bold')

# =============================================================================
# -----------------------------------------------------------------------------
# Set up logging handler
# -----------------------------------------------------------------------------
# =============================================================================

class TextHandler(logging.Handler):

    def __init__(self, text, gui):

        # Run the regular Handler __init__
        logging.Handler.__init__(self)

        # Store a reference to the Text it will log to
        self.text = text
        self.gui = gui

    def emit(self, record):
        msg = self.format(record)
        def append():
            # Add message
            self.text.configure(state='normal')
            self.text.insert(tk.END, msg + '\n')
            self.text.configure(state='disabled')

            # Autoscroll to the bottom
            self.text.yview(tk.END)

        # This is necessary because we can't modify the Text from other threads
        self.text.after(0, append)

        # Force the gui to update
        self.gui.update()

# =============================================================================
# -----------------------------------------------------------------------------
# Main GUI Program
# -----------------------------------------------------------------------------
# =============================================================================

class mygui(tk.Tk):

    def __init__(self, parent, *args, **kwargs):
        
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.root = parent

        # Cause exceptions to report in a new window
        tk.Tk.report_callback_exception = self.report_callback_exception

        # Close program on closure of window
        self.root.protocol("WM_DELETE_WINDOW", self.handler)
        
        # Read in pre-saved volcano settings
        self.volc_choice = ['--select--']
        self.volc_data = {}
        try:
            
            with open('bin/volcano_data.txt', 'r') as r:

                # Read header line
                r.readline()

                # Read data
                data = r.readlines()

                # Unpack and save to dictionary
                for line in data:

                    # Read each line
                    v_name, v_lon, v_lat, t_diff = line.strip().split('\t')

                    # Remove white space from parameters
                    v_name = v_name.strip()
                    v_lon  = v_lon.strip()
                    v_lat  = v_lat.strip()
                    t_diff = t_diff.strip()

                    # Add to choise array and populate volcano dictionary
                    self.volc_choice.append(v_name)
                    self.volc_data[v_name] = [v_lon, v_lat, t_diff]

        except  FileNotFoundError:
            print('No input file found')

            # Set all values to zero
            self.volc_data['--select--'] = [0,0,0]

        # Build the GUI
        self.build_gui()
        
# =============================================================================
#   Build Gui
# =============================================================================

    def build_gui(self):
        
        # Set default button Style
        ttk.Style().configure('TButton', width=10, height=20, relief="flat")
        
# =============================================================================
#         Build Containers
# =============================================================================
        
        # Create frame for setup
        setup_frame = tk.LabelFrame(self.root, text='Setup', font=LARGE_FONT)
        setup_frame.grid(row=0, column=0, padx=10, pady=10, sticky="W",
                           columnspan=2)
        
        # Create the frame for the volcano settings
        volc_frame = tk.LabelFrame(self.root, text='Volcano', font=LARGE_FONT)
        volc_frame.grid(row=1, column=0, padx=10, pady=10, sticky="W")
        
        # Create frame for control
        control_frame = tk.LabelFrame(self.root, text='Control', 
                                      font=LARGE_FONT)
        control_frame.grid(row=2, column=0, padx=10, pady=10, sticky="W")
        
        # Create the graph frame
        graph_frame = ttk.Frame(self.root)
        graph_frame.grid(row=1, column=1, padx=10, pady=10, rowspan=2)

# =============================================================================
#       Add Program Controls
# =============================================================================

        # Create input for the iFit output file
        so2_ftype = [("iFit output", ".csv")]
        # self.so2_path = tk.StringVar(value = 'No file selected')
        self.so2_path = tk.StringVar(value = 'output.csv')
        make_input(frame=setup_frame,
                   text='SO2 File',
                   var=self.so2_path, 
                   input_type='Entry',
                   row=0, column=0,
                   width=80)
        ttk.Button(setup_frame, text = "Browse", width = 8,
                   command = lambda: select_files(single_file=True,
                                                  holder=self.so2_path,
                                                  filetypes=so2_ftype)
                   ).grid(row=0, column=2, padx=5, pady=5, sticky='W')

        # Create input for the GPS file
        gps_ftype = [("GPS file", ".txt")]
        # self.gps_path = tk.StringVar(value = 'No file selected')
        self.gps_path = tk.StringVar(value = 'gps-data.txt')
        make_input(frame=setup_frame,
                   text='GPS File',
                   var=self.gps_path, 
                   input_type='Entry',
                   row=1, column=0,
                   width=80)
        ttk.Button(setup_frame, text = "Browse", width = 8,
                   command = lambda: select_files(single_file=True,
                                                  holder=self.gps_path,
                                                  filetypes=gps_ftype)
                   ).grid(row=1, column=2, padx=5, pady=5, sticky='W')

        # Create input for the output file
        out_ftype = [("text file", ".txt")]
        self.out_path = tk.StringVar(value = 'No file selected')
        make_input(frame=setup_frame,
                   text='Output File',
                   var=self.out_path, 
                   input_type='Entry',
                   row=2, column=0,
                   width=80)
        ttk.Button(setup_frame, text = "Browse", width = 8,
                   command = lambda: select_files(single_file=True,
                                                  holder=self.out_path,
                                                  save_flag=True,
                                                  filetypes=out_ftype)
                   ).grid(row=2, column=2, padx=5, pady=5, sticky='W')

        # Control whether poor fits are removed
        self.de_spike = tk.BooleanVar(setup_frame, value=True)
        make_input(frame=setup_frame,
                   text='Remove Bad\nSpectra?',
                   var=self.de_spike,
                   input_type='Checkbutton',
                   row=0, column=3)
        
        # Make a button to read in the data
        ttk.Button(setup_frame, text="Import Traverse", width=20, 
                   command=self.import_data
                   ).grid(row=1, column=3, padx=5, pady=5, sticky='W', 
                          columnspan=2)
                          
        # Set the flux units to calculate
        flux_unit_choice = ['kg/s', 'kg/s', 't/day']
        self.flux_units = tk.StringVar(setup_frame, value=flux_unit_choice[0])
        make_input(frame = setup_frame,
                   text = 'Flux Units:',
                   row = 2, column = 3,
                   var = self.flux_units,
                   input_type = 'OptionMenu',
                   options = flux_unit_choice,
                   width = 5,
                   sticky = 'W') 

# =============================================================================
#       Add Volcano and Wind Controls
# =============================================================================

        # Update parameters on choice
        def volc_update(event):
            self.volc_lon.set(self.volc_data[self.volc_name.get()][0])
            self.volc_lat.set(self.volc_data[self.volc_name.get()][1])
            self.time_diff.set(self.volc_data[self.volc_name.get()][2])

        # Create choice of imported volcanos
        self.volc_name = tk.StringVar(volc_frame, value=self.volc_choice[0])
        make_input(frame = volc_frame,
                   text = 'Volcano:',
                   row = 0, column = 0,
                   var = self.volc_name,
                   input_type = 'OptionMenu',
                   options = self.volc_choice,
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
        wind_speed_u.grid(row=1, column=4, padx=5, pady=5, sticky='EW')

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
        wind_error_u.grid(row=2, column=4, padx=5, pady=5, sticky='EW')

# =============================================================================
#         Set up text output and control frame
# =============================================================================

        # Make a button to read in the data
        ttk.Button(control_frame, text="Calculate Flux", width=20, 
                   command=self.calc_flux
                   ).grid(row=0, column=0, padx=5, pady=5, sticky='W')
                          
        # Add text widget to display logging info
        st = ScrolledText.ScrolledText(control_frame, state='disabled',
                                       width = 60, height = 12)
        st.configure(font='TkFixedFont')
        st.grid(row=1, column=0, padx=10, pady=10, sticky="NW", columnspan=5)

        # Create textLogger
        text_handler = TextHandler(st, self)

        # Logging configuration
        logging.basicConfig(filename='calc_flux.log',
                            level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S',
                            filemode='w')

        # Add the handler to logger
        logger = logging.getLogger()
        logger.addHandler(text_handler)

# =============================================================================
#       Add the graph
# =============================================================================
        
        # Create figure to hold the graphs
        plt.rcParams.update({'font.size': 8} )
        self.fig = plt.figure(figsize = (6.4, 4.0))

        # Create the canvas to hold the graph in the GUI
        self.canvas = FigureCanvasTkAgg(self.fig, graph_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=1, column=0, padx=10,
                                         pady=10, sticky='NW')

        # Create plot axes
        self.ax = self.fig.add_subplot(111)
        self.fig.subplots_adjust(bottom=0.3)

        self.line0, = self.ax.plot([], [], 'C0-')
        self.line1, = self.ax.plot([], [], 'C1-', lw=2)
        self.vline0 = self.ax.axvline(0, color = 'r')
        self.vline1 = self.ax.axvline(1, color = 'r')
        self.plume_cent, = self.ax.plot([], [], 'C2o', ms=10)
        self.ax.set_xlabel('Time (decimal hours)', fontsize=10)
        self.ax.set_ylabel('SO$_2$ SCD (molec cm$^{-2}$)', fontsize=10)

        # Create sliders
        self.ax_lo = self.fig.add_axes([0.12, 0.15, 0.78, 0.03])
        self.slide_lo = Slider(self.ax_lo, 'Low Bound', 0, 1, valinit = 0.05)
        self.ax_hi = self.fig.add_axes([0.12, 0.1, 0.78, 0.03])
        self.slide_hi = Slider(self.ax_hi, 'High Bound', 0, 1, valinit = 0.95)
        
        self.slide_lo.on_changed(self.slider_update)
        self.slide_hi.on_changed(self.slider_update)

        # Add matplotlib toolbar above the plot canvas
        toolbar_frame = tk.Frame(graph_frame, bg = 'black')
        toolbar_frame.grid(row=0, column=0, sticky='W', padx=5, pady=5)
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        toolbar.update()
        
# =============================================================================
#   Slider Update
# =============================================================================
        
    def slider_update(self, val):

        # Get slider values
        pos = self.slide_lo.val, self.slide_hi.val

        # Get graph limits
        lims = self.ax.get_xlim()

        # Caclulate position from limits and slider value
        self.t0 = lims[0] + (lims[1] - lims[0]) * pos[0]
        self.t1 = lims[0] + (lims[1] - lims[0]) * pos[1]
        self.vline0.set_xdata(self.t0)
        self.vline1.set_xdata(self.t1)
        
# =============================================================================
#   Import Data
# =============================================================================
        
    def import_data(self):

        logging.info('Reading in SO2 data...')
        
        # Read in the SO2 results
        so2_df = pd.read_csv(self.so2_path.get(), parse_dates=['Time'])
        
        self.so2_time = np.array([t.hour + t.minute/60 + t.second/3600 
                                  for t in so2_df['Time']])
        self.so2_scd = so2_df['SO2'].to_numpy()
        self.so2_err = so2_df['SO2_err'].to_numpy()
        
        if self.de_spike:
            mask = so2_df['fit_quality'] != 1
            self.so2_scd = np.ma.masked_where(mask, self.so2_scd)
            self.so2_err = np.ma.masked_where(mask, self.so2_err)
            
        logging.info('Done')
        
        logging.info('Importing GPS data...')
        
        gps_df = pd.read_table(self.gps_path.get(), sep='\t', 
                               parse_dates=['time'])
        
        self.gps_time = np.array([t.hour + t.minute/60 + t.second/3600 
                                  for t in gps_df['time']])
        self.lat = gps_df['latitude'].to_numpy()
        self.lon = gps_df['longitude'].to_numpy()
        self.alt = gps_df['altitude (m)'].to_numpy()
        
        logging.info('Done')
        
        # Update graph
        self.line0.set_xdata(self.so2_time)
        self.line0.set_ydata(self.so2_scd)

        self.vline0.set_xdata(self.so2_time[0])
        self.vline1.set_xdata(self.so2_time[-1])
        
        self.ax.relim()
        self.ax.autoscale_view()
        
        self.canvas.draw()
        self.update()
        
        logging.info('Traverse data imported')
        
        # Create a loop counter
        self.loop = 1
 
# =============================================================================
#   Calculate flux       
# =============================================================================

    def calc_flux(self):

        # If wind error is absolute, turn into fractional error
        if self.error_unit.get() == 'abs':
            wind_error = self.wind_error.get() / self.wind_speed.get()
        else:
            wind_error = self.wind_error.get() / 100
            
        # Check that the sliders are the right way around:
        if self.t1 < self.t0:
            logging.warn('Bounds the wrong way around, reversing')
            self.t0, self.t1 = self.t1, self.t0
            
        # Pull the selected data from the so2 data
        so2_idx = np.where(np.logical_and(self.so2_time >= self.t0,
                                          self.so2_time <= self.t1))
        
        so2_time = self.so2_time[so2_idx]
        so2_scd = self.so2_scd[so2_idx]
        so2_err = self.so2_err[so2_idx]

        # Find the centre of mass of the plume
        cum_so2_scd = np.cumsum(so2_scd)
        peak_idx = np.abs(cum_so2_scd - cum_so2_scd[-1]).argmin()
        
        # Correct for the time diffeence between the SO2 and GPS data
        gps_time = self.gps_time + self.time_diff.get()

        # Extract the relevant gps data and interpolate onto the so2 grid
        lat = griddata(gps_time, self.lat, so2_time)
        lon = griddata(gps_time, self.lon, so2_time)

        # Calculate the plume bearing
        volc_loc = [self.volc_lat.get(), self.volc_lon.get()]
        peak_loc = [lat[peak_idx], lon[peak_idx]]
        plume_dist, plume_bearing = haversine(volc_loc, peak_loc)
        
        # Calculate the distance and bearing of each measurement vector
        vect = [haversine([lat[i-1],lon[i-1]], [lat[i],lon[i]])
                for i in range(len(lat) - 1)]
        
        # Unpack the distance and bearing from the vectors
        trav_dist, trav_bearing = np.asarray(vect).T
        
        # Correct the distance for the angle between the travel and plume 
        # vectors
        corr_dist = np.multiply(trav_dist, np.sin(plume_bearing-trav_bearing))
        
        # If wind speed is in knotts, convert to m/s
        if self.wind_unit.get() == 'knots':
            wind_speed = self.wind_speed.get() * 0.5144444444
        else:
            wind_speed = self.wind_speed.get()
            
        # Convert so2 amounts from molecules.cm-2 to molecules.m-2
        so2_molec_per_m2 = so2_scd * 1.0e4
        
        # Multiply by the distance moved and sum
        so2_molec_per_m = np.sum(np.multiply(so2_molec_per_m2[1:], corr_dist))
        
        # Multiply by the wind speed
        so2_molec_per_s = so2_molec_per_m * wind_speed
        
        # Convert to moles
        so2_moles_per_s = so2_molec_per_s / 6.022e23
        
        # Convert to kg/s. Molar mass of SO2 is 64.066 g/mole
        so2_kg_per_s = so2_moles_per_s * 0.064066
        
        # Convert to t/day if required
        if self.flux_units.get() == 't/day':
            flux = abs(so2_kg_per_s * 60*60*24 / 1000.0)
        else:
            flux = abs(so2_kg_per_s)
        
        # Calculate the Flux Error
        tot_so2_err = np.sum(np.power(so2_err, 2)) ** 0.5
        frac_so2_err = tot_so2_err / np.sum(so2_scd)
        
        # Combine with the wind speed error
        frac_err = ( (frac_so2_err)**2 + (wind_error)**2 )**0.5
        flux_err = flux * frac_err
        
        msg=f'Peak {self.loop} flux: {flux:.02f} +/- {flux_err:.02f}'+\
            f' {self.flux_units.get()}'
        logging.info(msg)
        
        self.loop += 1

# =============================================================================
#   Exiting and Errors
# =============================================================================

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
        text = 'Would you like to quit?'
        message = tkMessageBox.askquestion('Exit',
                                           message = text,
                                           type = 'yesno')

        if message == 'yes':
            self.quit()

        if message == 'no':
            pass
        

if __name__ == '__main__':

    root = tk.Tk()
    mygui(root)
    root.mainloop()