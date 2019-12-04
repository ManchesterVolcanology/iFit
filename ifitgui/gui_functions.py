# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:36:03 2019

@author: mqbpwbe2
"""

import os
import logging
import numpy as np
import pandas as pd
import datetime as dt
from tkinter import filedialog as fd
import tkinter.messagebox as tkMessageBox

from ifit.parameters import Parameters
from ifit.model_setup import model_setup
from ifit.spectral_analysis import pre_process, fit_spectrum
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

    # Pull the settings from the GUI
    logging.info('Reading model settings')
    settings = {'w_lo':          gui.widgets['w_lo'].get(),
                'w_hi':          gui.widgets['w_hi'].get(),
                'model_spacing': gui.widgets['model_spacing'].get(),
                'model_padding': gui.widgets['model_padding'].get(),
                'dark_flag':     gui.widgets['dark_flag'].get(),
                'flat_flag':     gui.widgets['flat_flag'].get(),
                'stray_flag':    gui.widgets['stray_flag'].get(),
                'ils_type':      gui.widgets['ils_mode'].get(),
                'frs_path':      gui.widgets['frs_path'].get(),
                'flat_path':     gui.widgets['flat_path'].get(),
                'ils_path':      gui.widgets['ils_path'].get()}

    # Pull the parameters from the parameter table
    params = Parameters()
    logging.info('Generating model parameters')

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
    if settings['ils_type'] == 'Manual':
        params.add('fwem', value=gui.widgets['fwem'].get(),
                    vary=gui.widgets['fwem_fit'].get())
        params.add('k', value=gui.widgets['k'].get(), 
                   vary=gui.widgets['k_fit'].get())
        params.add('a_w', value=gui.widgets['a_w'].get(), 
                   vary=gui.widgets['a_w_fit'].get())
        params.add('a_k', value=gui.widgets['a_k'].get(), 
                   vary=gui.widgets['a_k_fit'].get())

    settings['gas_data'] = {}

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
            settings['gas_data'][name] = [xpath]

    # Build the forward model
    logging.info('Setting up forward model')
    common = model_setup(settings)

    # Report fitting parameters
    logging.info(params.pretty_print(cols=['name', 'value', 'vary']))

    # Add the parameters to the common and make a copy of the initial params
    common['params'] = params
    common['x0'] = params.valueslist()

    # Begin the analysis
    spec_list = ['iFit','Master.Scope','Spectrasuite','Basic']
    if gui.widgets['spec_type'].get() in spec_list:
        spectra_loop(gui, common, settings)

    elif gui.widgets['spec_type'].get() in ['FLAME', 'OpenSO2']:
        scan_loop(gui, common, settings)

#==============================================================================
#================================ spectra_loop ================================
#==============================================================================

def spectra_loop(gui, common, settings):

    # Make a list of column names
    cols = ['File', 'Number', 'Time']
    for par in common['params']:
        cols += [par, f'{par}_err']
    cols += ['fit_quality']
    
    # Pull the spectra type and file paths
    spec_fnames = gui.spec_fnames
    dark_fnames = gui.dark_fnames
    spec_type = gui.widgets['spec_type'].get()
    
    # Make a dataframe to hold the fit results
    df = pd.DataFrame(index = np.arange(len(spec_fnames)), columns = cols)

    # Read in the dark spectra
    logging.info('Reading dark spectra')
    x, common['dark'] = average_spectra(dark_fnames, spec_type)

    # Find the fit and stray windows
    logging.info('Calculating the fit window')
    common['fit_idx'] = np.where(np.logical_and(x > settings['w_lo'],
                                                x < settings['w_hi']))

    if common['stray_flag']:
        logging.info('Calculating the stray light window')
        common['stray_idx'] = np.where(np.logical_and(x > 280,  x < 290))

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

        # Pre-process the spectrum before the fit
        logging.debug(f'Pre-processing spectrum {fname}')
        spectrum = pre_process([x,y], common)
        grid, spec = spectrum

        # Fit the spectrum
        logging.debug(f'Fitting spectrum {fname}')
        fit_result = fit_spectrum(spectrum,
                                  common,
                                  update_params=gui.widgets['update_flag'].get(),
                                  resid_limit=gui.widgets['resid_limit'].get(),
                                  calc_od = [gui.widgets['graph_param'].get()]
                                  )

        # Add the the results dataframe
        row = [fname, info['spec_no'], info['time']]
        for par in fit_result.params.values():
            row += [par.fit_val, par.fit_err]
        row += [fit_result.nerr]

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
                plot_x = df['Number']
                plot_y = df[gui.widgets['graph_param'].get()]

                # Trim if required
                if gui.widgets['scroll_flag'].get():
                    if gui.loop > gui.widgets['scroll_amt'].get():
                        diff = gui.loop - gui.widgets['scroll_amt'].get()
                        plot_x = plot_x[diff:]
                        plot_y = plot_y[diff:]

            except KeyError:
                plot_x = []
                plot_y = []

            # Organise data to plot
            #        x_data, y_data
            data = [[grid,   spec],
                    [grid,   fit_result.fit],
                    [x,      y],
                    [grid,   fit_result.resid],
                    [grid,   fit_result.meas_od['SO2']],
                    [grid,   fit_result.synth_od['SO2']],
                    [plot_x, plot_y]
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


    try:
        # Save the results
        df.to_csv(gui.widgets['save_path'].get())

    except PermissionError:

        # Open save dialouge
        text = 'Cannot save output: Permission Denied'
        message = tkMessageBox.askquestion('Select new save location?',
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

def scan_loop(gui, common, settings):

    # Make a list of column names
    cols = ['Number', 'Time', 'Motor_Pos']
    for par in common['params']:
        cols += [par, f'{par}_err']
    cols += ['fit_quality']

    # Get the wavelength grid of the spectrometer
    x = np.loadtxt(gui.widgets['wl_calib'].get())

    # Find the fit and stray windows
    logging.info('Calculating the fit window')
    common['fit_idx'] = np.where(np.logical_and(x > settings['w_lo'],
                                                x < settings['w_hi']))

    if common['stray_flag']:
        logging.info('Calculating the stray light window')
        common['stray_idx'] = np.where(np.logical_and(x>280, x<290))

        if len(common['stray_idx'][0]) == 0:
            logging.warn('No stray light window found, ' + \
                         'disabling stray light correction')
            common['stray_flag'] = False

    # Create a loop counter
    gui.loop = 0
    spec_fnames = gui.spec_fnames
    spec_type = gui.spec_type

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
            common['dark'] = spec_block[0]

            # Cycle through the scan block
            for n, y in enumerate(spec_block[1:]):

                # Extract spectrum info
                if spec_type == 'FLAME':
                    n_aq, h, m, s, motor_pos = info_block[:,n+1]
                if spec_type == 'OpenSO2':
                    n_aq, h, m, s, motor_pos, int_t, coadds = info_block[n+1]

                # Pre-process the spectrum before the fit
                logging.debug(f'Pre-processing spectrum {fname}')
                spectrum = pre_process([x,y], common)
                grid, spec = spectrum

                # Fit the spectrum
                logging.debug(f'Fitting spectrum {fname}')
                fit_result = fit_spectrum(spectrum,
                                          common,
                                          gui.widgets['update_flag'].get(),
                                          gui.widgets['resid_limit'].get(),
                                          [gui.widgets['graph_param'].get()]
                                          )

                # Add to the results dataframe
                time = dt.time(int(h), int(m), int(s))
                row = [n_aq, time, motor_pos]
                for par in fit_result.params.values():
                    row += [par.fit_val, par.fit_err]
                row += [fit_result.nerr]

                df.loc[n] = row

                # Update numerical outputs
                try:
                    key = gui.graph_param.get()
                    gui.last_amt.set(f'{df[key][n]:.03g}')
                    gui.last_err.set(f'{df[key+"_err"][n]:.03g}')

                except KeyError:
                    gui.last_amt.set('-')
                    gui.last_err.set('-')

                # Plot the graphs
                if gui.widgets['graph_flag'].get():

                    # Pick the parameter to plot
                    try:
                        plot_x = df['Number']
                        plot_y = df[gui.widgets['graph_param'].get()]

                        # Trim if required
                        if gui.widgets['scroll_flag'].get():
                            if n > gui.widgets['scroll_amt'].get():
                                diff = n - gui.widgets['scroll_amt'].get()
                                plot_x = plot_x[diff:]
                                plot_y = plot_y[diff:]

                    except KeyError:
                        plot_x = []
                        plot_y = []

                    # Organise data to plot
                    #        x_data, y_data
                    data = [[grid,   spec],
                            [grid,   fit_result.fit],
                            [x,      y],
                            [grid,   fit_result.resid],
                            [grid,   fit_result.meas_od['SO2']],
                            [grid,   fit_result.synth_od['SO2']],
                            [plot_x, plot_y]
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