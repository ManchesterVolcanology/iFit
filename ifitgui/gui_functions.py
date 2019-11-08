# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:36:03 2019

@author: mqbpwbe2
"""

import os
import logging
import numpy as np
import pandas as pd
from tkinter import filedialog as fd

from ifit.parameters import Parameters
from ifit.model_setup import model_setup
from ifit.spectral_analysis import FitSpectrum, pre_process
from ifit.load_spectra import read_spectrum, average_spectra

#==============================================================================
#=============================== analysis_loop ================================
#==============================================================================

def analysis_loop(gui):
    '''Function to handle the analysis loop'''

    # Set status
    gui.status.set('Loading')

    # Set the stopping flag to False
    gui.stop_flag = False

    # Pull the settings from the GUI
    logging.info('Reading model settings')
    settings = {'w_lo':          gui.w_lo.get(),
                'w_hi':          gui.w_hi.get(),
                'model_spacing': gui.model_spacing.get(),
                'model_padding': gui.model_padding.get(),
                'dark_flag':     gui.dark_flag.get(),
                'flat_flag':     True,
                'stray_flag':    gui.stray_flag.get(),
                'ils_type':      'File',
                'frs_path':      gui.frs_path.get(),
                'flat_path':     gui.flat_path.get(),
                'ils_path':      gui.ils_path.get()}

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
    msg = params.pretty_print(cols = ['name', 'value', 'vary'])
    logging.info(msg)

    common['params'] = params

    # Make a list of column names
    cols = ['File', 'Number', 'Time']
    for par in params:
        cols += [par, f'{par}_err']
    cols += ['fit_quality']

    n_spec = np.arange(len(gui.spec_fnames))

    # Make a dataframe to hold the fit results
    df = pd.DataFrame(index = n_spec, columns = cols)

    # Read in the dark spectra
    logging.info('Reading dark spectra')
    x, common['dark'], read_err = average_spectra(gui.dark_fnames)

    # Find the fit and stray windows
    logging.info('Calculating the fit window')
    common['fit_idx'] = np.where(np.logical_and(x > settings['w_lo'],
                                                x < settings['w_hi']))

    if common['stray_flag']:
        logging.info('Calculating the stray light window')
        common['stray_idx'] = np.where(np.logical_and(x > 280,  x < 290))

    # Check the parameter to graph is in the parameter list
    plot_param_flag = True
    if gui.graph_param.get() not in params:
        logging.warn(f'Graphing parameter {gui.graph_param.get()} not valid')
        plot_param_flag = False

    gui.loop = 0

    logging.info('Beginning analysis loop')
    gui.status.set('Analysing')

    # Begin the analysis
    while not gui.stop_flag:

        # Get the filename
        fname = gui.spec_fnames[gui.loop]

        # Read in the spectrum
        logging.debug(f'Reading in spectrum {fname}')
        x, y, info, read_err = read_spectrum(fname)

        # Pre-process the spectrum before the fit
        logging.debug(f'Pre-processing spectrum {fname}')
        spectrum = pre_process([x,y], common)
        grid, spec = spectrum

        # Fit the spectrum
        logging.debug(f'Fitting spectrum {fname}')
        fit_result = FitSpectrum(spectrum, common,
                                 update_params = gui.update_flag.get())

        # Add the the results dataframe
        row = [fname, info['spec_no'], info['time']]
        for par in fit_result.params.values():
            row += [par.fit_val, par.fit_err]
        row += [fit_result.nerr]

        df.loc[gui.loop] = row

        # Update numerical outputs
        gui.last_amt.set(f'{df[gui.graph_param.get()][gui.loop]:.03g}')
        gui.last_err.set(f'{df[gui.graph_param.get() + "_err"][gui.loop]:.03g}')

        # Plot the graphs
        if gui.graph_flag.get():

            # Pick the parameter to plot
            if plot_param_flag:
                plot_x = df['Number']
                plot_y = df[gui.graph_param.get()]

            # Trim if required
            if gui.scroll_flag.get():
                if gui.loop > gui.graph_data_n.get():
                    diff = gui.loop - gui.graph_data_n.get()
                    plot_x = plot_x[diff:]
                    plot_y = plot_y[diff:]

            # Organise data to plot
            #        x_data, y_data
            data = [[grid,   spec],
                    [grid,   fit_result.fit],
                    [x,      y],
                    [grid,   fit_result.resid],
                    [[],     []],
                    [[],     []],
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
        if gui.loop == len(gui.spec_fnames)-1:
            gui.stop_flag = True
            gui.status.set('Standby')
            logging.info('Analysis finished!')

        # Update the loop counter
        gui.loop += 1



    # Save the results
    df.to_csv('fit_results.csv')



#==============================================================================
#================================= get_fpaths =================================
#==============================================================================

# Function to select spectra to analyse
def get_fpaths(spec_list, entry = None):

    '''
    Function to get plume spectra file paths

    **Parameters**

    spec_list : list
        The list to hold the spectra file names

    entry : tk.Entry
        The entry to which to print the number of spectra selected

    **Returns**

    spec_list : list
        The list of the spectra to analyse
    '''

    # Open dialouge to get files
    fpaths = fd.askopenfilenames()

    if fpaths != '':
        spec_list = []
        for fname in fpaths:
            spec_list.append(str(fname))

        # Save output to input
        if entry != None:
            entry.set(f'{len(spec_list)} spectra selected')

    logging.info(f'{len(spec_list)} spectra selected')

    return spec_list

#==============================================================================
#================================ select_files ================================
#==============================================================================

def select_files(single_file=False, holder=None, entry=None):

    '''
    Function to get plume spectra file paths

    **Parameters**

    file_list : list
        The list to hold the spectra file names

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
                entry.set(f'{len(holder)} spectra selected')

        return fpaths

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