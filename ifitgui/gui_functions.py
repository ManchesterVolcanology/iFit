import time
import os
import logging
import numpy as np
import pandas as pd
import datetime as dt
# from queue import Queue
# from threading import Thread
from multiprocessing import Process, Queue
from tkinter import filedialog as fd
import tkinter.messagebox as tkMessageBox

from ifit.parameters import Parameters
from ifit.spectral_analysis import Analyser
from ifit.load_spectra import read_spectrum, average_spectra, read_scan
from .spectrometers import VSpectrometer


# =============================================================================
# Plot fit results
# =============================================================================

def plot_fit_results(gui, spectrum, fit_result, tseries):

    # Pick the parameter to plot
    try:
        param = gui.widgets['graph_param'].get()
        plot_x = tseries[0]
        plot_y = tseries[1]
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
    data = [[fit_result.grid, fit_result.spec],
            [fit_result.grid, fit_result.fit],
            spectrum,
            [fit_result.grid, fit_result.resid],
            [fit_result.grid, meas_od],
            [fit_result.grid, synth_od],
            [plot_x,          plot_y]
            ]

    gui.figure.update_plots(data)

    gui.canvas.draw()
    gui.update()


# =============================================================================
# Generate analyser
# =============================================================================

def generate_analyser(gui):

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
            name = row[0].get().strip()
            value = float(row[1].get())
            vary = bool(row[2].get())
            xpath = row[3].get().strip()
            params.add(name=name,
                       value=value,
                       vary=vary,
                       xpath=xpath)

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
                        fit_window=[gui.widgets['fit_lo'].get(),
                                    gui.widgets['fit_hi'].get()],
                        frs_path=gui.widgets['frs_path'].get(),
                        model_padding=gui.widgets['model_padding'].get(),
                        model_spacing=gui.widgets['model_spacing'].get(),
                        flat_flag=gui.widgets['flat_flag'].get(),
                        flat_path=gui.widgets['flat_path'].get(),
                        stray_flag=gui.widgets['stray_flag'].get(),
                        stray_window=[gui.widgets['stray_lo'].get(),
                                      gui.widgets['stray_hi'].get()],
                        dark_flag=gui.widgets['dark_flag'].get(),
                        ils_type=gui.widgets['ils_mode'].get(),
                        ils_path=gui.widgets['ils_path'].get(),
                        despike_flag=gui.widgets['despike_flag'].get(),
                        spike_limit=gui.widgets['spike_limit'].get())

    # Report fitting parameters
    logging.info(params.pretty_print(cols=['name', 'value', 'vary']))

    gui.analyser = analyser


# =============================================================================
# analysis_loop
# =============================================================================

def analysis_loop(gui, rt_flag):
    """Function to handle the analysis loop for post analysis"""

    # Determine whether the analyser needs to be generated
    if rt_flag and not gui.fitting_flag:
        gui.analyser = None

    else:
        generate_analyser(gui)

    # Record the start time
    gui.start_time = time.time()

    # Launch the analyser loop corresponding to the analysis type

    # Real-time analysis
    if rt_flag:
        rt_analysis_loop(gui)

    # Scanner analysis
    elif gui.widgets['spec_type'].get() in ['FLAME', 'OpenSO2']:
        scan_loop(gui)

    # Normal spectra analysis
    else:
        spectra_loop(gui)


# =============================================================================
# spectra_loop
# =============================================================================

def spectra_loop(gui):

    # Make a list of column names
    cols = ['File', 'Number', 'Time']
    for par in gui.analyser.params:
        cols += [par, f'{par}_err']
    cols += ['fit_quality', 'int_lo', 'int_hi', 'int_av']

    # Pull the spectra type and file paths
    spec_fnames = gui.spec_fnames
    dark_fnames = gui.dark_fnames
    spec_type = gui.widgets['spec_type'].get()

    # Make a dataframe to hold the fit results
    df = pd.DataFrame(index=np.arange(len(spec_fnames)), columns=cols)

    # Read in the dark spectra
    logging.info('Reading dark spectra')
    x, gui.analyser.dark_spec = average_spectra(dark_fnames, spec_type)

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
        resid_type = gui.widgets['resid_type'].get()
        int_limit = [gui.widgets['lo_int_limit'].get(),
                     gui.widgets['hi_int_limit'].get()]
        graph_p = [gui.widgets['graph_param'].get()]
        interp_meth = gui.widgets['interp_method'].get()

        # Fit the spectrum
        logging.debug(f'Fitting spectrum {fname}')
        fit_result = gui.analyser.fit_spectrum(spectrum=[x, y],
                                               update_params=update_flag,
                                               resid_limit=resid_limit,
                                               resid_type=resid_type,
                                               int_limit=int_limit,
                                               calc_od=graph_p,
                                               interp_method=interp_meth)

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
            tseries = [df['Number'], df[gui.widgets['graph_param'].get()]]
            plot_fit_results(gui, [x, y], fit_result, tseries)

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
        text = 'Cannot save output: Permission Denied\n' + \
               'Select new save location?'
        message = tkMessageBox.askquestion('Save Error',
                                           message=text,
                                           type='yesno')

        if message == 'yes':
            file_io(single_file=True, holder=gui.save_path, save_flag=True,
                    filetypes=(('Comma Separated', '.csv')))
            df.to_csv(gui.widgets['save_path'].get())

        if message == 'no':
            pass


# =============================================================================
# scan_loop
# =============================================================================

def scan_loop(gui):

    # Make a list of column names
    cols = ['Number', 'Time', 'Motor_Pos']
    for par in gui.analyser.params:
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
            gui.analyser.dark_spec = spec_block[0]

            # Cycle through the scan block
            for n, y in enumerate(spec_block[1:]):

                # Extract spectrum info
                if spec_type == 'FLAME':
                    n_aq, h, m, s, motor_pos = info_block[:, n+1]
                if spec_type == 'OpenSO2':
                    n_aq, h, m, s, motor_pos, int_t, coadds = info_block[n+1]

                # Pull processing settings from the GUI
                update_flag = gui.widgets['update_flag'].get()
                resid_limit = gui.widgets['resid_limit'].get()
                resid_type = gui.widgets['resid_type'].get()
                int_limit = [gui.widgets['lo_int_limit'].get(),
                             gui.widgets['hi_int_limit'].get()]
                graph_p = [gui.widgets['graph_param'].get()]
                interp_meth = gui.widgets['interp_method'].get()

                # Fit the spectrum
                logging.debug(f'Fitting spectrum {fname}')
                fit_res = gui.analyser.fit_spectrum(spectrum=[x, y],
                                                    update_params=update_flag,
                                                    resid_limit=resid_limit,
                                                    resid_type=resid_type,
                                                    int_limit=int_limit,
                                                    calc_od=graph_p,
                                                    pre_process=True,
                                                    interp_method=interp_meth)

                # Add to the results dataframe
                time = dt.time(int(h), int(m), int(s))
                row = [n_aq, time, motor_pos]
                for par in fit_res.params.values():
                    row += [par.fit_val, par.fit_err]
                row += [fit_res.nerr, fit_res.int_lo, fit_res.int_hi,
                        fit_res.int_av]

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
                    tseries = [df['Number'],
                               df[gui.widgets['graph_param'].get()]]
                    plot_fit_results(gui, [x, y], fit_res, tseries)

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

                logging.warn(f'Permission Error: Unable to save file {fname}.')

        else:
            logging.warn(f'Error reading file {fname}')

        # Check if analysis is finished
        if gui.loop == len(gui.spec_fnames)-1:
            gui.stop_flag = True
            gui.status.set('Standby')
            logging.info('Analysis finished!')

        # Update the loop counter
        gui.loop += 1


# =============================================================================
# Real Time Analysis
# =============================================================================

def rt_analysis_loop(gui, buffer=1000):
    """Handles real time spectra acquisition and analysis"""

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    gui.spec.fpath = 'Masaya_Traverse/spectrum_00366.txt'
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Log the start of acquisition
    logging.info('Beginning spectra acquisition')
    gui.status.set('Acquiring')

    # Create holder for the previous spectrum
    gui.last_spectrum = None

    # Create a loop counter and spectrum number
    gui.loop = 0
    spec_no = 0

    # Create the output spectra folder if it is missing
    savepath = gui.widgets["rt_save_path"].get() + 'spectra/'
    if not os.path.isdir(savepath):
        os.makedirs(savepath)

    # Begin the acquisition loop
    while not gui.stop_flag:

        # Generate the spectrum filename
        spec_fname = f'{savepath}spectrum_{spec_no:05d}.txt'

        # Check the file doesn't already exist
        while os.path.isfile(spec_fname):
            spec_no += 1
            spec_fname = f'{savepath}spectrum_{spec_no:05d}.txt'

        # Just acquire spectra
        if not gui.fitting_flag:

            # Take a spectrum
            x, y = gui.spec.get_spectrum(spec_fname)

            # Add to the plot
            if gui.widgets['graph_flag'].get():

                data = [[[], []], [[], []], [x, y], [[], []], [[], []],
                        [[], []], [[], []]]
                gui.figure.update_plots(data)

                gui.canvas.draw()
                gui.update()

        # Analyse in real-time
        else:

            if gui.last_spectrum is None:

                # Take a spectrum
                x, y = gui.spec.get_spectrum(spec_fname)
                gui.last_spectrum = [x, y]

            else:

                # Set the dark spectrum for the analyser
                gui.analyser.dark_spec = gui.dark_spectrum

                # Create a queue to hold results
                results_queue = Queue()

                # Pull processing settings from the GUI
                update_flag = gui.widgets['update_flag'].get()
                resid_limit = gui.widgets['resid_limit'].get()
                resid_type = gui.widgets['resid_type'].get()
                int_limit = [gui.widgets['lo_int_limit'].get(),
                             gui.widgets['hi_int_limit'].get()]
                graph_p = [gui.widgets['graph_param'].get()]
                interp_meth = gui.widgets['interp_method'].get()

                fit_args = (gui.last_spectrum, update_flag, resid_limit,
                            resid_type, int_limit, graph_p, True, interp_meth,
                            results_queue)

                # Set one process to acquire the spectrum
                # t1 = Process(target=gui.spec.get_spectrum,
                #              args=(spec_fname, results_queue))

                # Set another process to analyse the previous spectrum
                t2 = Process(target=gui.analyser.fit_spectrum,
                             args=fit_args)

                # Launch the processes
                # t1.start()
                t2.start()

                # Join once finished
                # t1.join()
                t2.join()

                # Get the results from the processes
                output = {}
                while not results_queue.empty():
                    result = results_queue.get()
                    output[result[0]] = result[1]

                spectrum = gui.last_spectrum
                fit_result = output['fit_result']

                # Update the numerical outputs
                key = gui.widgets['graph_param'].get()
                amt = fit_result.params[key].fit_val
                err = fit_result.params[key].fit_err
                gui.last_amt.set(f'{amt:.03g}')
                gui.last_err.set(f'{err:.03g}')

                # Plot the graphs
                if gui.widgets['graph_flag'].get():
                    tseries = [[], []]
                    plot_fit_results(gui, spectrum, fit_result, tseries)

                gui.last_spectrum = output['spectrum']

                print('mark')

        gui.update()


# =============================================================================
# Connect to spectrometer
# =============================================================================

def connect_spectrometer(gui):
    """Connects or dissconnects to the spectrometer"""

    if not gui.connect_flag.get():

        # Connect to the spectrometer
        gui.spec = VSpectrometer(integration_time=gui.widgets["IntTime"].get(),
                                 coadds=gui.widgets["Coadds"].get())

        # Update the GUI
        gui.spec_id.set(gui.spec.serial_number)
        gui.con_b.config(text='Disconnect')

        # Create a holder for the dark spectra
        gui.dark_spectrum = np.zeros(gui.spec.pixels)

        gui.connect_flag.set(True)

    else:
        # Disconnect the spectrometer
        gui.spec.close()

        # Update the GUI
        gui.spec_id.set('Not connected')
        gui.con_b.config(text='Connect')

        gui.connect_flag.set(False)


# =============================================================================
# Test Spectrum
# =============================================================================

def test_spectrum(gui):
    """Measure a test spectrum and display"""

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    gui.spec.fpath = 'Masaya_Traverse/spectrum_00000.txt'
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Read the spectrum
    # try:
    x, y = gui.spec.get_spectrum()

    # Plot the spectrum
    data = [[[], []], [[], []], [x, y], [[], []], [[], []], [[], []],
            [[], []]]
    gui.figure.update_plots(data)

    gui.canvas.draw()
    gui.update()

    # except AttributeError:
    #     logging.warning('No spectrometer connected!')


# =============================================================================
# Read Darks
# =============================================================================

def read_darks(gui):
    """Read in the dark spectra"""

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    gui.spec.fpath = 'Masaya_Traverse/dark.txt'
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Update the status
    gui.status.set('Acquiring')

    # Get the number of darks
    ndarks = gui.widgets["ndarks"].get()
    logging.info(f'Reading {ndarks} dark spectra')

    # Generate the save location for the dark spectra
    dark_path = gui.widgets["rt_save_path"].get() + 'dark/'

    # Generate the output folder, generating a new one if it already exists
    if not os.path.isdir(dark_path):
        os.makedirs(dark_path)

    else:
        n = 0
        while os.path.isdir(dark_path):
            n += 1
            dark_path = gui.widgets["rt_save_path"].get() + f'dark({n})/'
        os.makedirs(dark_path)

    try:
        # Create a holder for the averaged dark spectrum
        npixels = gui.spec.pixels
        dark_arr = np.zeros([ndarks, npixels])

        # Read the dark spectra
        for i in range(ndarks):
            fname = f'{dark_path}spectrum_{i:05d}.txt'
            x, y = gui.spec.get_spectrum(fname=fname)
            dark_arr[i] = y

            # Update the progress bar
            gui.progress['value'] = (i+1)/ndarks * 100

            # Make GUI show updates
            gui.update()

        # Average the spectra
        gui.dark_spectrum = np.average(dark_arr, axis=0)

        # Plot the spectrum
        data = [[[], []], [[], []], [x, gui.dark_spectrum], [[], []], [[], []],
                [[], []], [[], []]]
        gui.figure.update_plots(data)

        gui.canvas.draw()
        gui.update()

    except AttributeError:
        logging.warning('No spectrometer connected!')

    gui.status.set('Standby')


# =============================================================================
# Toggle fitting
# =============================================================================

def toggle_fitting(gui):
    """Toggle fitting on or off"""

    # If fitting is ON --> turn it OFF
    if gui.fitting_flag:
        gui.fitting_flag = False
        gui.fit_b.config(text='Fitting OFF', bg='red')
        logging.info('Fitting turned off')

    # If fitting is OFF --> turn it ON
    else:
        gui.fitting_flag = True
        gui.fit_b.config(text='Fitting ON', bg='green')
        logging.info('Fitting turned on')


# =============================================================================
# file_io
# =============================================================================

def file_io(folder=False, single_file=False, holder=None, entry=None,
            save_flag=False, filetypes=None):
    """Function to handle file dialouges for input/output selection.

    Parameters
    ----------
    folder : bool, optional
        If True then the dialouge selects a folder, otherwise a file/files are
        selected. Overrides single_file. Default is False
    single_file : bool, optional
        Controls whether multiple file selection is supported or not. Default
        is False
    holder : tkinter Variable or None, optional
        The tk variable to hold the returned file path. Ignored if None.
        Default is None
    entry : tk Object or None, optional
        A tk Label or Entry to hold the number of selected files if selecting
        multiple. Ignored if None or if single_file is True. Default is None
    save_flag : bool, optional
        If True then a "save as" dialouge is used. Default is False
    filetypes : list of tuples or None, optional
        The default file extensions to offer the user. Each tuple is a pair of
        strings with a description and the extension. Default is None

    Returns
    -------
    str
        The file path(s) selected by the user. If the file dialouge is
        cancelled then an empty string is returned.
    """

    if filetypes is None:
        filetypes = [("all files", "*.*")]
    else:
        filetypes.append(("all files", "*.*"))

    # Get the cwd
    cwd = os.getcwd().replace("\\", "/")

    # Folder ==================================================================

    # Open dialouge for a folder
    if folder:
        fpath = fd.askdirectory(initialdir=cwd, title="Select folder")

        if fpath != '' and holder is not None:
            fpath += '/'
            # Check if in the same cwd. If so trim the file path
            if cwd in fpath:
                fpath = fpath[len(cwd)+1:]

            holder.set(fpath)

        return fpath

    # Single File =============================================================

    # Open dialouge to get a single file
    if single_file:
        if save_flag:
            fname = fd.asksaveasfilename(initialdir=cwd,
                                         title="Select file",
                                         filetypes=filetypes)
        else:
            fname = fd.askopenfilename(filetypes=filetypes)

        if fname != '' and holder is not None:

            # Check if in the same cwd. If so trim the file path
            if cwd in fname:
                fname = fname[len(cwd)+1:]

            holder.set(fname)

        return fname

    # Multiple Files ==========================================================

    # Open dialouge to get multiple files
    else:
        fnames = fd.askopenfilenames(filetypes=filetypes)

        if fnames != '':

            if holder is not None:
                # Clear the holder list
                holder.clear()

                for fname in fnames:
                    holder.append(str(fname))

            # Save output to input
            if entry is not None:
                entry.set(f'{len(holder)} files selected')

        return fnames


# =============================================================================
# stop
# =============================================================================

def stop(self):
    """Stops the analysis loop"""
    self.stop_flag = True
    logging.info('Analysis stoped')
