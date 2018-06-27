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
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from pandas import read_csv

from ifit_lib.find_nearest import extract_window
from ifit_lib.read_gps import read_txt_gps
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
        tk.Tk.wm_title(self, 'calc_ratio 2.0')
        tk.Tk.iconbitmap(self, default = 'data_bases/icon.ico')

        # Create control Frame
        cont_frame = ttk.Frame(self, relief = 'groove')
        cont_frame.grid(row=0, column=0, padx=10, pady=10, rowspan=2, sticky="NW")
        
        # Create frame to hold graphs
        graph_frame = ttk.Frame(self, relief = 'groove')
        graph_frame.grid(row=0, column=1, padx=10, pady=10, rowspan=10, sticky="N")
        graph_frame.columnconfigure(index=0, weight=1)
        graph_frame.rowconfigure(index = 0, weight = 1)
        
        # Create frame to hold text output
        text_frame = ttk.Frame(self, relief = 'groove')
        text_frame.grid(row=2, column=0, padx=10, pady=10, sticky="NW")
        
        # Make graph widget scale        
        mygui.columnconfigure(index=1, weight=1, self = self)
        mygui.rowconfigure(index = 5, weight = 1, self = self)
        
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
        self.text_box = tkst.ScrolledText(text_frame, width = 50, height = 10)
        self.text_box.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W',
                           columnspan = 4)
        self.text_box.insert('1.0', 'Welcome to calc_ratio! Written by Ben Esse\n\n')
        
#========================================================================================
#=====================================Build Graph========================================
#========================================================================================
        
        # Create figure to hold the graphs
        plt.rcParams.update({'font.size': 8} )
        self.fig = plt.figure(figsize = (7,5))
        
        # Create the canvas to hold the graph in the GUI
        self.canvas = FigureCanvasTkAgg(self.fig, graph_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=1, column=0, columnspan = 3, padx=10, 
                                         pady=10, sticky = 'NW')
        
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
        self.t0 = 0
        self.t1 = -1
            
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
        
        # Add matplotlib toolbar above the plot canvas
        toolbar_frame = tk.Frame(graph_frame, bg = 'black')  
        toolbar_frame.grid(row=0,column=0, sticky = 'NW', padx = 6, pady = 10)                             
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
         
        # Create frame
        meas_frame = tk.LabelFrame(cont_frame, text='Analysis Controls', font=LARG_FONT)
        meas_frame.grid(row=2, column=0, padx=10, pady=6, columnspan = 2, sticky="NW")
        
        # Create spin boxes to select gases to calculate ratios
        options = ['so2','bro','o3','no2']
        
        self.gas0 = tk.StringVar(meas_frame, value = options)
        self.gas0_l = tk.Label(meas_frame, text = 'x-axis:', font = NORM_FONT)
        self.gas0_l.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.gas0_c = ttk.OptionMenu(meas_frame, self.gas0, *options)
        self.gas0_c.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        self.gas1 = tk.StringVar(meas_frame, value = options)
        self.gas1_l = tk.Label(meas_frame, text = 'y-axis:', font = NORM_FONT)
        self.gas1_l.grid(row = 1, column = 0, padx = 5, pady = 5, sticky = 'W')
        self.gas1_c = ttk.OptionMenu(meas_frame, self.gas1, *options)
        self.gas1_c.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = 'W')
        
        # Create button to calculate ratio
        begin_b = ttk.Button(meas_frame, text="Calculate Ratio", command=self.calc_ratio)
        begin_b.grid(row = 2, column = 0, padx = 5, pady = 5, rowspan=2)
        
        # Create a button to exit
        exit_b = ttk.Button(meas_frame, text = 'Exit', command = self.quit)
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
#==================================Read Traverse data====================================
#========================================================================================
         
    def read_trav_data(self):
        
        # Reset loop counter
        common['loop'] = 0
        
        self.text_output('Reading traverse data...', add_line = False)
           
        # Read data file
        try:
            common['data'] = read_csv(str(self.ifit_path.get()))

            # Load required data into arrays
            txt_time    = common['data']['Time'] 
            
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
        date_string = str(common['data']['Date'][0]).replace('/', '-')
        common['out_folder']='Results/calc_ratio/' + date_string + '/'
        if not os.path.exists(common['out_folder']):
            os.makedirs(common['out_folder'])
            
        # Create output .csv file if it doesn't exist
        common['out_fname'] = common['out_folder'] + 'ratio_out.csv'
        
        with open(common['out_fname'], 'w') as writer:
            writer.write('calc_ratio.py output file \n')
        
        # Update graph
        self.line0.set_xdata(common['time'])
        self.line0.set_ydata(common['data']['so2_amt'])
        
        self.vline0.set_xdata(common['time'][0])
        self.vline1.set_xdata(common['time'][-1])
        
        self.fill.remove()
        lo_bound = np.subtract(common['data']['so2_amt'], common['data']['so2_amt_e'])
        hi_bound = np.add(common['data']['so2_amt'], common['data']['so2_amt_e'])
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
#====================================Calculate Ratio=====================================
#========================================================================================
    
    def calc_ratio(self):
        
        # Read in origional data set
        try:
            data       = common['data']
            time       = common['time']
            gps_time   = common['gps_time']
            lon        = common['lon']
            lat        = common['lat']
            volc_lon   = self.volc_lon.get()
            volc_lat   = self.volc_lat.get()
            
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
        try:
            time_grid, idx0, idx1 = extract_window(time, self.t0, self.t1)
            idx1 += 1
            
        except AttributeError:
            self.text_output('Please set sliders')
            return
        
        # Catch error if bounds are the wrong way around
        if idx0 > idx1:
            self.text_output('Bounds wrong way around, please try again')
            return
        
        # Cut selected gases
        common['gas0_label'] = self.gas0.get() + '_amt'
        common['gas1_label'] = self.gas1.get() + '_amt'
        common['x']    = data[common['gas0_label']][idx0:idx1]
        common['y']    = data[common['gas1_label']][idx0:idx1]
        common['xerr'] = data[common['gas0_label'] + '_e'][idx0:idx1]
        common['yerr'] = data[common['gas1_label'] + '_e'][idx0:idx1] 
        
#========================================================================================
#=====================================Calculate ratio====================================
#========================================================================================
        
        # Perform liner regression
        m, c, R, p, stderr = linregress(common['x'], common['y'])
        common['m'], common['c'], common['R^2'] = m, c, R**2
        
        # Create line to plot
        lin_fit = np.polyval([m,c], common['x'])
        
#========================================================================================
#======================Interpolate GPS data onto spectra time grid=======================
#========================================================================================

        # Interpolate the GPS data onto the so2 time grid
        common['modlat'] = griddata(gps_time,lat,time)
        common['modlon'] = griddata(gps_time,lon,time)
        
        # Get plume position and time
        peak_idx = cog(data['so2_amt'])
        common['peak_idx'] = peak_idx
        common['volc_loc']  = str(volc_lat) + ' N, ' + str(volc_lon) + ' E'
        common['plume_loc'] = str(common['modlat'][peak_idx]) + ' N, ' + \
                              str(common['modlon'][peak_idx]) + ' E'
        common['plume_time'] = common['time'][idx0:idx1]
        common['p_start_t'] = str(common['time'][0])
        common['p_centr_t'] = str(common['time'][common['peak_idx']])
        common['p_stop_t']  = str(common['time'][-1]) 
                           
#========================================================================================
#=============================Plot selected SO2 and GPS data=============================
#========================================================================================
           
        # Make labels
        labels = {}
        labels['so2_amt'] = r'SO$_2$ (molecules$\cdot$cm$^-$$^2$)'
        labels['bro_amt'] = r'BrO (molecules$\cdot$cm$^-$$^2$)'
        labels['o3_amt']  = r'O$_3$ (molecules$\cdot$cm$^-$$^2$)'
        
        # Create plot grid
        fig, ax = plt.subplots()
        plt.rcParams.update({'font.size': 12})
        gs  = gridspec.GridSpec(1,1)
        ax0 = plt.subplot(gs[0,0]) 
        
        ax0.errorbar(common['x'], common['y'], fmt = 'o', xerr = common['xerr'],
                     yerr = common['yerr'], label = '_nolegend_', zorder = 1, capsize=2, 
                     markersize = 3, alpha = 0.5)
        
        graph_label = 'y = '+"{0:.2g}".format(m)+'x + '+"{0:.2g}".format(c)+'\n'+\
                      'R$^2$ = ' + "{0:.2g}".format(R**2)
                      
        ax0.plot(common['x'], lin_fit, 'r', label = graph_label, zorder = 2)
        
        try:
            ax0.set_xlabel(labels[common['gas0_label']])
            ax0.set_ylabel(labels[common['gas1_label']])
        except KeyError:
            ax0.set_xlabel(common['gas0_label'])
            ax0.set_ylabel(common['gas1_label'])
        
        '''
        ax0.plot(common['x'], np.divide(common['y'], common['x']))
        ax0.set_ylabel('BrO/SO2')
        ax0.set_xlabel('SO2 column amount ( molec/cm2)')
        '''
        plt.legend(loc=0)
        plt.tight_layout()   

        # Make graph in a new window
        make_graph(fig)
            

        
#========================================================================================
#=======================================Make Graph=======================================
#========================================================================================
            
            
# Define function to open another window to show output graphs
def make_graph(fig):
    popup = tk.Tk()
    tk.Tk.wm_title(popup, 'Results')
    
    def save(fig):
        
        # Save the figure
        fig.savefig(common['out_folder'] + 'ratio_' + str(common['loop']) + '.png')
        
        # Save the revelant data
        with open(common['out_fname'], 'a') as w:
            
            # Write header info
            w.write('Peak number,' + str(common['loop'])       + '\n' \
                    'Ratio,' + str(common['m'])                + '\n' \
                    'Intercept,' + str(common['c'])            + '\n' \
                    'R^2,' + str(common['R^2'])                + '\n' \
                    'Volcano Long/Lat,' + common['volc_loc']   + '\n' \
                    'Plume position,' + common['plume_loc']    + '\n' \
                    'Plume Start Time,' + common['p_start_t']  + '\n' \
                    'Plume Centre Time,' + common['p_centr_t'] + '\n' \
                    'Plume Stop Time,' + common['p_stop_t']    + '\n')
        
            # Write column headers            
            w.write('Julian Time (GMT),Latitude,Longitude,' + common['gas0_label'] + ','\
                    + common['gas0_label'] + '_e,' + common['gas1_label'] + ','\
                    + common['gas1_label'] + '_e\n')
            
            # Change to numpy arrays to avoid index errors
            x = np.asarray(common['x'])
            xerr = np.asarray(common['xerr'])
            y = np.asarray(common['y'])
            yerr = np.asarray(common['yerr'])
            plume_time = np.asarray(common['plume_time'])
            
            # Write information to the file
            for i in range(1, len(x)):
                w.write(str(plume_time[i])           + ',' + \
                        str(common['modlat'][i])     + ',' + \
                        str(common['modlon'][i])     + ',' + \
                        str(x[i])                    + ',' + \
                        str(xerr[i])                 + ',' + \
                        str(y[i])                    + ',' + \
                        str(yerr[i])                 + '\n')
        
                w.write('\n')
                
        # Add to the loop counter
        common['loop'] += 1
        
        # Close the window
        popup.destroy()
    
    # Create the canvas to hold the graph in the GUI
    canvas = FigureCanvasTkAgg(fig, popup)
    canvas.show()
    canvas.get_tk_widget().grid(row=1, column=0, padx=10, pady = 10, columnspan = 2)
      
    # Add matplotlib toolbar above the plot canvas
    toolbar_frame = tk.Frame(popup, bg = 'black')  
    toolbar_frame.grid(row=0,column=0, sticky = 'NW', padx = 6, pady = 10, columnspan=2)                             
    toolbar = NavigationToolbar2TkAgg(canvas, toolbar_frame)
    toolbar.update()
    
    # Create a button to save 
    save_b = ttk.Button(popup, text="Save", command = lambda: save(fig))
    save_b.grid(row = 2, column = 0, padx = 5, pady = 5)
    
    # Create a button to exit
    cancel_b = ttk.Button(popup, text = 'Cancel', command = popup.destroy)
    cancel_b.grid(row = 2, column = 1, padx = 5, pady = 5)

# Run the App!     
mygui().mainloop()