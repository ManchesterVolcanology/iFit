# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:36:03 2019

@author: mqbpwbe2
"""
import time
import os
import logging
import numpy as np
import pandas as pd
import datetime as dt
from tkinter import filedialog as fd
import tkinter.messagebox as tkMessageBox

from ifit.parameters import Parameters
from ifit.spectral_analysis import Analyser
from ifit.load_spectra import read_spectrum, average_spectra, read_scan

#==============================================================================
#=============================== analysis_loop ================================
#==============================================================================

def analysis_loop(gui):
    '''Function to handle the analysis loop'''

    # Set status
    gui.status.set('Loading')
    logging.info('\nStart\n')

    # Set the stopping flag to False
    gui.stop_flag = False

    # Pull the parameters from the parameter table
    params = Parameters()
    logging.info('Generating model parameters')
    
    # Add parameters from GUI table
    for row in gui.gastable._params:
        if row != []:
            name  = row[0].get().strip()
            value = float(row[1].get())
            vary  = bool(row[2].get())
            xpath = row[3].get().strip()
            params.add(name  = name,
                       value = value,
                       vary  = vary,
                       xpath = xpath)

    # Add polynomial parameters
    for i, p in enumerate(gui.polytable._params):
        if p != []:
            params.add(f'bg_poly{i}', value=p[1].get(), vary=p[2].get())

    for i, p in enumerate(gui.offsettable._params):
        if p != []:
            params.add(f'offset{i}', value=p[1].get(), vary=p[2].get())

    for i, p in enumerate(gui.shifttable._params):
        if p != []:
            params.add(f'shift{i}', value=p[1].get(), vary=p[2].get())

    # Check if ILS is in the fit
    if gui.widgets['ils_mode'].get() == 'Manual':
        params.add('fwem', value=gui.widgets['fwem'].get(),
                    vary=gui.widgets['fwem_fit'].get())
        params.add('k', value=gui.widgets['k'].get(),
                   vary=gui.widgets['k_fit'].get())
        params.add('a_w', value=gui.widgets['a_w'].get(),
                   vary=gui.widgets['a_w_fit'].get())
        params.add('a_k', value=gui.widgets['a_k'].get(),
                   vary=gui.widgets['a_k_fit'].get())

    # Generate the analyser
    analyser = Analyser(params=params,
                        fit_window    = [gui.widgets['fit_lo'].get(), 
                                         gui.widgets['fit_hi'].get()],
                        frs_path      = gui.widgets['frs_path'].get(),
                        model_padding = gui.widgets['model_padding'].get(),
                        model_spacing = gui.widgets['model_spacing'].get(),
                        flat_flag     = gui.widgets['flat_flag'].get(),
                        flat_path     = gui.widgets['flat_path'].get(),
                        stray_flag    = gui.widgets['stray_flag'].get(),
                        stray_window  = [gui.widgets['stray_lo'].get(),
                                         gui.widgets['stray_hi'].get()],
                        dark_flag     = gui.widgets['dark_flag'].get(),
                        ils_type      = gui.widgets['ils_mode'].get(),
                        ils_path      = gui.widgets['ils_path'].get())

    # Report fitting parameters
    logging.info(params.pretty_print(cols=['name', 'value', 'vary']))

    # Record the start time
    gui.start_time = time.time()

    # Begin the analysis
    spec_list = ['iFit','Master.Scope','Spectrasuite','Basic']
    if gui.widgets['spec_type'].get() in spec_list:
        spectra_loop(gui, analyser)

    elif gui.widgets['spec_type'].get() in ['FLAME', 'OpenSO2']:
        scan_loop(gui, analyser)

#==============================================================================
#================================ spectra_loop ================================
#==============================================================================

def spectra_loop(gui, analyser):

    # Make a list of column names
    cols = ['File', 'Number', 'Time']
    for par in analyser.params:
        cols += [par, f'{par}_err']
    cols += ['fit_quality', 'int_lo', 'int_hi', 'int_av']

    # Pull the spectra type and file paths
    spec_fnames = gui.spec_fnames
    dark_fnames = gui.dark_fnames
    spec_type = gui.widgets['spec_type'].get()

    # Make a dataframe to hold the fit results
    df = pd.DataFrame(index = np.arange(len(spec_fnames)), columns = cols)

    # Read in the dark spectra
    logging.info('Reading dark spectra')
    x, analyser.dark_spec = average_spectra(dark_fnames, spec_type)

    # Create a loop counter
    gui.loop = 0

    logging.info('Beginning analysis loop')
    gui.status.set('Analysing')

    # Begin the analysis
    while not gui.stop_flag:

        # Get the filename
        fname = spec_fnames[gui.loop]

        # Read in the spectrum
        logging.debug(f'Reading in spectrum {fname}')
        x, y, info, read_err = read_spectrum(fname, spec_type)
                    
        # Pull processing settings from the GUI
        update_flag = gui.widgets['update_flag'].get()
        resid_limit = gui.widgets['resid_limit'].get()
        resid_type  = gui.widgets['resid_type'].get()
        int_limit   = [gui.widgets['lo_int_limit'].get(),
                       gui.widgets['hi_int_limit'].get()]
        graph_p     = [gui.widgets['graph_param'].get()]
        interp_meth = gui.widgets['interp_method'].get()

        # Fit the spectrum
        logging.debug(f'Fitting spectrum {fname}')
        fit_result = analyser.fit_spectrum(spectrum=[x,y],
                                           update_params = update_flag,
                                           resid_limit   = resid_limit,
                                           resid_type    = resid_type,
                                           int_limit     = int_limit,
                                           calc_od       = graph_p,
                                           interp_method = interp_meth)

        # Add the the results dataframe
        row = [fname, info['spec_no'], info['time']]
        for par in fit_result.params.values():
            row += [par.fit_val, par.fit_err]
        row += [fit_result.nerr, fit_result.int_lo, fit_result.int_hi,
                fit_result.int_av]

        df.loc[gui.loop] = row

        # Update numerical outputs
        try:
            key = gui.widgets['graph_param'].get()
            gui.last_amt.set(f'{df[key][gui.loop]:.03g}')
            gui.last_err.set(f'{df[key+"_err"][gui.loop]:.03g}')

        except KeyError:
            gui.last_amt.set('-')
            gui.last_err.set('-')

        # Plot the graphs
        if gui.widgets['graph_flag'].get():

            # Pick the parameter to plot
            try:
                param = gui.widgets['graph_param'].get()
                plot_x = df['Number']
                plot_y = df[gui.widgets['graph_param'].get()]
                meas_od = fit_result.meas_od[param]
                synth_od = fit_result.synth_od[param]

                # Trim if required
                if gui.widgets['scroll_flag'].get():
                    if gui.loop > gui.widgets['scroll_amt'].get():
                        diff = gui.loop - gui.widgets['scroll_amt'].get()
                        plot_x = plot_x[diff:]
                        plot_y = plot_y[diff:]

            except KeyError:
                plot_x = []
                plot_y = []
                meas_od = []
                synth_od = []

            # Organise data to plot
            #        x_data, y_data
            data = [[fit_result.grid, fit_result.spec ],
                    [fit_result.grid, fit_result.fit  ],
                    [x,               y               ],
                    [fit_result.grid, fit_result.resid],
                    [fit_result.grid, meas_od         ],
                    [fit_result.grid, synth_od        ],
                    [plot_x,          plot_y          ]
                    ]

            gui.figure.update_plots(data)

            gui.canvas.draw()
            gui.update()

        # Update the progress bar
        gui.progress['value'] = (gui.loop+1)/len(gui.spec_fnames) * 100

        # Make GUI show updates
        gui.update()

        # Check if analysis is finished
        if gui.loop == len(spec_fnames)-1:
            gui.stop_flag = True
            gui.status.set('Standby')
            logging.info('Analysis finished!')

        # Update the loop counter
        gui.loop += 1

    # Report the time taken to analyse
    gui.end_time = time.time()
    logging.info(f'Analysis time: {gui.end_time-gui.start_time:.02f}')

    try:
        # Save the results
        df.to_csv(gui.widgets['save_path'].get())

    except PermissionError:

        # Open save dialouge
        text = 'Cannot save output: Permission Denied\nSelect new save location?'
        message = tkMessageBox.askquestion('Save Error',
                                           message = text,
                                           type = 'yesno')

        if message == 'yes':
            select_save(holder=gui.save_path)
            df.to_csv(gui.widgets['save_path'].get())

        if message == 'no':
            pass

#==============================================================================
#================================= scan_loop ==================================
#==============================================================================

def scan_loop(gui, analyser):

    # Make a list of column names
    cols = ['Number', 'Time', 'Motor_Pos']
    for par in analyser.params:
        cols += [par, f'{par}_err']
    cols += ['fit_quality', 'int_lo', 'int_hi', 'int_av']

    # Get the wavelength grid of the spectrometer
    x = np.loadtxt(gui.widgets['wl_calib'].get())

    # Create a loop counter
    gui.loop = 0
    spec_fnames = gui.spec_fnames
    spec_type = gui.widgets['spec_type'].get()

    logging.info('Beginning analysis loop')
    gui.status.set('Analysing')

    # Begin the analysis
    while not gui.stop_flag:

        # Get the filename
        fname = spec_fnames[gui.loop]

        # Read in the scan file
        logging.debug(f'Reading in scan {fname}')
        err, info_block, spec_block = read_scan(fname, spec_type)

        if not err:

            # Make a dataframe to hold the fit results
            df = pd.DataFrame(index=np.arange(len(spec_fnames)), columns=cols)

            # Get the dark
            analyser.dark_spec = spec_block[0]

            # Cycle through the scan block
            for n, y in enumerate(spec_block[1:]):

                # Extract spectrum info
                if spec_type == 'FLAME':
                    n_aq, h, m, s, motor_pos = info_block[:,n+1]
                if spec_type == 'OpenSO2':
                    n_aq, h, m, s, motor_pos, int_t, coadds = info_block[n+1]
                    
                # Pull processing settings from the GUI
                update_flag = gui.widgets['update_flag'].get()
                resid_limit = gui.widgets['resid_limit'].get()
                resid_type  = gui.widgets['resid_type'].get()
                int_limit   = [gui.widgets['lo_int_limit'].get(),
                               gui.widgets['hi_int_limit'].get()]
                graph_p     = [gui.widgets['graph_param'].get()]
                interp_meth = gui.widgets['interp_method'].get()

                # Fit the spectrum
                logging.debug(f'Fitting spectrum {fname}')
                fit_result = analyser.fit_spectrum(spectrum=[x,y],
                                                   update_params = update_flag,
                                                   resid_limit   = resid_limit,
                                                   resid_type    = resid_type,
                                                   int_limit     = int_limit,
                                                   calc_od       = graph_p,
                                                   pre_process   = True,
                                                   interp_method = interp_meth)

                # Add to the results dataframe
                time = dt.time(int(h), int(m), int(s))
                row = [n_aq, time, motor_pos]
                for par in fit_result.params.values():
                    row += [par.fit_val, par.fit_err]
                row += [fit_result.nerr, fit_result.int_lo, fit_result.int_hi,
                        fit_result.int_av]

                df.loc[n] = row

                # Update numerical outputs
                try:
                    key = gui.widgets['graph_param'].get()
                    gui.last_amt.set(f'{df[key][n]:.03g}')
                    gui.last_err.set(f'{df[key+"_err"][n]:.03g}')

                except KeyError:
                    gui.last_amt.set('-')
                    gui.last_err.set('-')

                # Plot the graphs
                if gui.widgets['graph_flag'].get():

                    # Pick the parameter to plot
                    try:
                        param = gui.widgets['graph_param'].get()
                        plot_x = df['Number']
                        plot_y = df[gui.widgets['graph_param'].get()]
                        meas_od = fit_result.meas_od[param]
                        synth_od = fit_result.synth_od[param]

                        # Trim if required
                        if gui.widgets['scroll_flag'].get():
                            if n > gui.widgets['scroll_amt'].get():
                                diff = n - gui.widgets['scroll_amt'].get()
                                plot_x = plot_x[diff:]
                                plot_y = plot_y[diff:]

                    except KeyError:
                        plot_x = []
                        plot_y = []
                        meas_od = np.full(fit_result.grid.shape, np.nan)
                        synth_od = np.full(fit_result.grid.shape, np.nan)


                    # Organise data to plot
                    #        x_data,          y_data
                    data = [[fit_result.grid, fit_result.spec ],
                            [fit_result.grid, fit_result.fit  ],
                            [x,               y               ],
                            [fit_result.grid, fit_result.resid],
                            [fit_result.grid, meas_od         ],
                            [fit_result.grid, synth_od        ],
                            [plot_x,          plot_y          ]
                            ]

                    gui.figure.update_plots(data)

                    gui.canvas.draw()
                    gui.update()

                # Update the progress bar
                gui.progress['value'] = (n+1)/len(spec_block) * 100
                gui.status.set(f'{gui.loop} / {len(gui.spec_fnames)}')

                # Make GUI show updates
                gui.update()

                # Check if stopped
                if gui.stop_flag:
                    break

            # Save the output
            try:

                # Get just the filename
                file_name = fname.split('/')[-1].split('.')[0]

                # Save the results
                df.to_csv(f'{gui.widgets["save_path"].get()}{file_name}.csv')

            except PermissionError:

                logging.warn(f'Unable to save file {fname}. File already open.')

        else:
            logging.warn(f'Error reading file {fname}')

        # Check if analysis is finished
        if gui.loop == len(gui.spec_fnames)-1:
            gui.stop_flag = True
            gui.status.set('Standby')
            logging.info('Analysis finished!')

        # Update the loop counter
        gui.loop += 1

#==============================================================================
#================================ select_files ================================
#==============================================================================

def select_files(single_file=False, holder=None, entry=None):

    '''
    Function to select files

    **Parameters**

    single_file : bool, optional, default = False
        Controls whether to select a single file or multiple

    holder : tk variable, optional, default = None
        The tk variable to asign the file name(s). Ignored if None.

    entry : tk.Entry
        The entry to which to print the number of spectra selected

    **Returns**

    spec_list : list or str
        The list of file paths, or a single file path
    '''

    # Get the cwd
    cwd = os.getcwd().replace('\\', '/')

    # Open dialouge to get a single file
    if single_file:
        fpath = fd.askopenfilename()

        if fpath != '' and holder != None:

            # Check if in the same cwd. If so trim the file path
            if cwd in fpath:
                fpath = fpath[len(cwd)+1:]

            holder.set(fpath)

        return fpath

    # Open dialouge to get multiple files
    else:
        fpaths = fd.askopenfilenames()

        if fpaths != '':

            if holder != None:
                # Clear the holder list
                holder.clear()

                for fname in fpaths:
                    holder.append(str(fname))

            # Save output to input
            if entry != None:
                entry.set(f'{len(holder)} files selected')

        return fpaths

#==============================================================================
#================================ select_files ================================
#==============================================================================

def select_save(holder=None):

    '''
    Function to select a folder

    **Parameters**

    holder : tk variable, optional, default = None
        The tk variable to asign the file name(s). Ignored if None.

    **Returns**

    spec_list : list or str
        The list of file paths, or a single file path
    '''

    # Get the cwd
    cwd = os.getcwd().replace('\\', '/')

    # Open a file dialouge to get a folder
    fpath = fd.asksaveasfilename(defaultextension='.csv')

    if fpath != None:

        # Check if in the same cwd. If so trim the file path
        if cwd in fpath:
            fpath = fpath[len(cwd)+1:]

        if holder != None:
            holder.set(fpath)

    return fpath

#==============================================================================
#=============================== import_params ================================
#==============================================================================

def import_params(fpath, table = None):
    '''Reads in the param file'''

    # Create a Parameters object
    params = Parameters()

    # Read in the file
    with open(fpath, 'r') as r:

        lines = r.readlines()

        # Skipping the header read each line
        for line in lines[1:]:
            cols = line.strip().split(',')
            name = cols[0].strip()
            value = float(cols[1])
            if cols[2].strip() == 'True':
                vary = True
            else:
                vary = False
            xpath = cols[3].strip()

            params.add(name, value, vary, xpath)

    if table != None:

        # Clear the table
        table.set_params(params)

    return params


#==============================================================================
#==================================== stop ====================================
#==============================================================================

def stop(self):
    '''Stops the analysis loop'''
    self.stop_flag = True
    logging.info('Analysis stoped')
