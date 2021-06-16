import os
import sys
import time
import logging
import traceback
import numpy as np
import pandas as pd
from functools import partial
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt, QObject, QRunnable, pyqtSignal, pyqtSlot
from PyQt5.QtWidgets import (QComboBox, QTextEdit, QLineEdit, QDoubleSpinBox,
                             QSpinBox, QCheckBox, QFileDialog, QPushButton,
                             QTableWidgetItem, QMenu, QTableWidget,
                             QPlainTextEdit, QHeaderView)

from ifit.parameters import Parameters
from ifit.spectral_analysis import Analyser
from ifit.load_spectra import read_spectrum, average_spectra
from ifit.spectrometers import Spectrometer


logger = logging.getLogger(__name__)


class QTextEditLogger(logging.Handler, QObject):
    """Record logs to the GUI."""

    appendPlainText = pyqtSignal(str)

    def __init__(self, parent):
        super().__init__()
        QObject.__init__(self)
        self.widget = QPlainTextEdit(parent)
        self.widget.setReadOnly(True)
        self.widget.setFont(QFont('Courier', 10))
        self.appendPlainText.connect(self.widget.appendPlainText)

    def emit(self, record):
        """Emit the log."""
        msg = self.format(record)
        self.appendPlainText.emit(msg)


# Create a worker signals object to handle worker signals
class WorkerSignals(QObject):
    """Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data
    progress
        `int` indicating % progress
    plotter
        `list` data to be plotted on the analysis graphs
    error
        `tuple` (exctype, value, traceback.format_exc() )
    status
        'str' status message for the GUI
    spectrum
        `tuple` spectrum to be displayed on the scope graph
    """

    finished = pyqtSignal()
    progress = pyqtSignal(int)
    plotter = pyqtSignal(list)
    error = pyqtSignal(tuple)
    status = pyqtSignal(str)
    spectrum = pyqtSignal(tuple)
    data = pyqtSignal(np.ndarray)


# Create a worker to handle QThreads
class Worker(QRunnable):
    """Worker thread.

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up

    Parameters
    ----------
    fn : function
        The function to run on the worker thread
    mode : str
        Flags which signals to use. Must be one of 'analyse' or 'acquire', for
        a spectral analysis or acquisition thread respectively

    Attributes
    ----------
    args : list
        Arguments to  pass to the function
    kwargs : dict
        Keyword arguments to pass to the function
    signals : WorkerSignals object
        The worker signals
    is_paused : bool
        State to show if the worker has been paused
    is_killed : bool
        State to show if the worker has been killed
    spec_fname : str
        File path to the last measured spectrum
    """

    def __init__(self, fn, mode, *args, **kwargs):
        super(Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.mode = mode
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Set killed and paused flags
        self.is_paused = False
        self.is_killed = False

        # Create a holder for the spectrum filepath
        self.spec_fname = None

        # Add the callbacks to our kwargs
        if 'analyse' in mode:
            self.kwargs['progress_callback'] = self.signals.progress
            self.kwargs['status_callback'] = self.signals.status
            self.kwargs['plot_callback'] = self.signals.plotter

        if 'acquire' in mode:
            self.kwargs['spectrum_callback'] = self.signals.spectrum
            self.kwargs['progress_callback'] = self.signals.progress
            self.kwargs['status_callback'] = self.signals.status

    @pyqtSlot()
    def run(self):
        """Initialise the runner function with passed args and kwargs."""
        # Retrieve args/kwargs here; and fire processing using them
        try:
            self.fn(self, self.mode, *self.args, **self.kwargs)
        except Exception:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))

        # Done
        self.signals.finished.emit()

    def pause(self):
        """Pause the analysis/acquisition."""
        if self.is_paused:
            self.is_paused = False
        else:
            self.is_paused = True

    def resume(self):
        """Resume the analysis/acquisition."""
        self.is_paused = False

    def kill(self):
        """Terminate the analysis/acquisition."""
        if self.is_paused:
            self.is_paused = False
        self.is_killed = True

    def set_spectrum(self, spec_fname):
        """Set the spectrum file to be analysed in real time analysis."""
        self.spec_fname = spec_fname


# =============================================================================
# Generate analyser
# =============================================================================

def generate_analyser(widgetData):
    """Generate the iFit analyser.

    Parameters
    ----------
    widgetData : dict
        Contains the program settings from the GUI

    Returns
    -------
    analyser : Analyser
        Returns the constructed iFit analyser
    """
    # Pull the parameters from the parameter table
    params = Parameters()
    logger.info('Generating model parameters')

    # Build the parameters from GUI tables
    for line in widgetData['gas_params']:
        params.add(name=line[0], value=line[1], vary=line[2], xpath=line[4],
                   plume_gas=line[3])
    for i, line in enumerate(widgetData['bgpoly_params']):
        params.add(name=f'bg_poly{i}', value=line[0], vary=line[1])
    for i, line in enumerate(widgetData['offset_params']):
        params.add(name=f'offset{i}', value=line[0], vary=line[1])
    for i, line in enumerate(widgetData['shift_params']):
        params.add(name=f'shift{i}', value=line[0], vary=line[1])

    # Check if ILS is in the fit
    if widgetData['ils_mode'] == 'Manual':
        params.add('fwem', value=float(widgetData['fwem']),
                   vary=widgetData['fwem_fit'])
        params.add('k', value=float(widgetData['k']),
                   vary=widgetData['k_fit'])
        params.add('a_w', value=float(widgetData['a_w']),
                   vary=widgetData['a_w_fit'])
        params.add('a_k', value=float(widgetData['a_k']),
                   vary=widgetData['a_k_fit'])

    # Add the light dilution factor
    params.add('LDF', value=float(widgetData['ldf']),
               vary=widgetData['ldf_fit'])

    # Get the bad pixels
    if widgetData['bad_pixels'] != '':
        bad_pixels = [int(i) for i in widgetData['bad_pixels'].split(',')]
    else:
        bad_pixels = []

    # Generate the analyser
    analyser = Analyser(params=params,
                        fit_window=[widgetData['fit_lo'],
                                    widgetData['fit_hi']],
                        frs_path=widgetData['frs_path'],
                        model_padding=widgetData['model_padding'],
                        model_spacing=widgetData['model_spacing'],
                        flat_flag=widgetData['flat_flag'],
                        flat_path=widgetData['flat_path'],
                        stray_flag=widgetData['stray_flag'],
                        stray_window=[widgetData['stray_lo'],
                                      widgetData['stray_hi']],
                        dark_flag=widgetData['dark_flag'],
                        ils_type=widgetData['ils_mode'],
                        ils_path=widgetData['ils_path'],
                        despike_flag=widgetData['despike_flag'],
                        spike_limit=widgetData['spike_limit'],
                        bad_pixels=bad_pixels)

    # Report fitting parameters
    logger.info(params.pretty_print(cols=['name', 'value', 'vary']))

    return analyser


# =============================================================================
# Analysis setup
# =============================================================================

def analysis_setup(worker, analysis_mode, widgetData, buffer_cols):
    """Set up the mode dependant analysis settings."""
    # By default there is no dark spectrum
    dark_spec = None

    if analysis_mode == 'post_analyse':

        # Set acquisition parameters
        spec_fnames = widgetData['spec_fnames'].split('\n')
        dark_fnames = widgetData['dark_fnames'].split('\n')
        spec_type = widgetData['spec_type']
        wl_calib_file = widgetData['wl_calib']
        buffer = Buffer(len(spec_fnames), buffer_cols)
        save_path = widgetData['save_path']
        write_mode = 'w'
        worker.signals.status.emit('Analysing')

        # Read the dark spectrum
        if widgetData['dark_flag']:

            if len(dark_fnames) == 0:
                logger.warning('No dark spectra selected, disabling dark '
                               + 'correction')

            else:
                x, dark_spec = average_spectra(dark_fnames, spec_type,
                                               wl_calib_file)

    elif analysis_mode == 'rt_analyse':

        # Set acquisition parameters
        spec_fnames = None
        dark_fnames = ['bin/.dark']
        spec_type = 'iFit'
        wl_calib_file = ''
        buffer = Buffer(2000, buffer_cols)
        save_path = widgetData['rt_save_path'] + '/iFit_rt_output.csv'
        worker.signals.status.emit('Acquiring')

        # Read the dark spectrum
        if widgetData['dark_flag']:

            try:
                dark_spec = np.loadtxt('bin/.dark')
            except OSError:
                logger.warning('No dark spectrum found, disabling dark '
                               + 'correction')

        # Set the output write format
        if os.path.isfile(save_path):
            write_mode = 'a'
        else:
            write_mode = 'w'

    return (spec_fnames, dark_spec, spec_type, wl_calib_file, buffer,
            save_path, write_mode)


# =============================================================================
# Spectra analysis loop
# =============================================================================

def analysis_loop(worker, analysis_mode, widgetData, progress_callback,
                  plot_callback, status_callback):
    """Control loop for spectral analysis.

    Parameters
    ----------
    worker : Worker
        The worker running the thread
    analysis_mode : str
        Controls how the analysis is handled based on whether performing real
        time or post-analysis. Must be either 'post_analyse' or rt_analyse
    widgetData : dict
        Contains the program settings from the GUI
    progrss_callback : Signal
        Reports progress to the progress bar
    plot_callback : Signal
        Reports the data to display on the analysis plots
    status_callback : Signal
        Reports the status of the program

    Returns
    -------
    None
    """
    # Create a loop counter
    loop = 0

    # Update the status
    status_callback.emit('Loading')

    # Generate the anlyser
    analyser = generate_analyser(widgetData)

    # Pull processing settings
    update_flag = widgetData['update_flag']
    resid_limit = widgetData['resid_limit']
    resid_type = widgetData['resid_type']
    int_limit = [widgetData['lo_int_limit'],
                 widgetData['hi_int_limit']]
    graph_p = [r[0] for r in widgetData['gas_params']]
    interp_meth = widgetData['interp_method']

    # Make a list of the fit parameters to form the columns in the buffer
    buffer_cols = ['Number']
    for par in analyser.params:
        buffer_cols += [par, f'{par}_err']
    buffer_cols += ['fit_quality']

    # Set mode dependent settings
    [spec_fnames, dark_spec, spec_type, wl_calib_file, buffer, save_path,
     write_mode] = analysis_setup(worker,
                                  analysis_mode,
                                  widgetData,
                                  buffer_cols)

    # Set the dark spectrum
    if analyser.dark_flag and dark_spec is not None:
        analyser.dark_spec = dark_spec
    else:
        analyser.dark_flag = False

    # Set status
    logger.info('Beginning analysis loop')

    # Open the output file
    with open(save_path, write_mode) as w:

        # If writing a new file, add the columns
        if write_mode == 'w':

            # Make a list of column names
            cols = ['File', 'Number', 'Time']
            for par in analyser.params:
                cols += [par, f'{par}_err']
            cols += ['fit_quality', 'int_lo', 'int_hi', 'int_av']

            # Write the header
            w.write(cols[0])
            for c in cols[1:]:
                w.write(f',{c}')
            w.write('\n')

        # Begin the analysis loop
        while not worker.is_killed:

            # Check if analysis is paused
            while worker.is_paused:
                time.sleep(0.01)

            # Get the spectrum filename
            if analysis_mode == 'post_analyse':
                fname = spec_fnames[loop]
            elif analysis_mode == 'rt_analyse':
                fname = worker.spec_fname
                while fname is None:
                    fname = worker.spec_fname
                worker.spec_fname = None

            # Read in the spectrum
            x, y, info, read_err = read_spectrum(fname, spec_type,
                                                 wl_calib_file)

            # Fit the spectrum
            fit_result = analyser.fit_spectrum(spectrum=[x, y],
                                               update_params=update_flag,
                                               resid_limit=resid_limit,
                                               resid_type=resid_type,
                                               int_limit=int_limit,
                                               calc_od=graph_p,
                                               interp_method=interp_meth)

            # Write spectrum info to file
            w.write(f'{fname},{info["spec_no"]},{info["time"]}')

            # Add results to the buffer dataframe and write results to file
            row = [info["spec_no"]]
            for par in fit_result.params.values():
                row += [par.fit_val, par.fit_err]
                w.write(f',{par.fit_val},{par.fit_err}')
            row += [fit_result.nerr]
            buffer.update(row)

            # Write fit quality info and start a new line
            w.write(f',{fit_result.nerr},{fit_result.int_lo},'
                    + f'{fit_result.int_hi},{fit_result.int_av}')
            w.write('\n')

            # Emit the progress
            if analysis_mode == 'post_analyse':
                progress_callback.emit(((loop+1)/len(spec_fnames))*100)

            # Emit graph data
            plot_callback.emit([fit_result, [x, y], buffer.df,
                                info["spec_no"]])

            # Check if analysis is finished
            if analysis_mode == 'post_analyse' and loop == len(spec_fnames)-1:
                logger.info('Analysis finished!')
                worker.kill()

            # Update loop counter
            loop += 1


# =============================================================================
# Acquisition loop
# =============================================================================

def acquire_spectra(worker, acquisition_mode, widgetData, spectrometer,
                    spectrum_callback, progress_callback, status_callback):
    """Loop to handle spectra acquisition.

    Parameters
    ----------
    worker : Worker
        The worker running the thread
    aquisition_mode : str
        How spectra acquisition is run:
            - acquire_single: measure a single test spectrum
            - acquire_darks: measure dark spectra
            - acquire_cont: continuous acquisition
    widgetData : dict
        Contains the program settings from the GUI
    spectrum_callback : Signal
        Reports the measured spectrum to plot on the scope plot
    progrss_callback : Signal
        Reports progress to the progress bar
    status_callback : Signal
        Reports the status of the program

    Returns
    -------
    None
    """
    # Update the status
    status_callback.emit('Acquiring')

    # Read single spectrum
    if acquisition_mode == 'acquire_single':

        # Read the spectrum
        spectrum, info = spectrometer.get_spectrum()

        spectrum_callback.emit((spectrum, info, True))

    # Read dark spectra
    if acquisition_mode == 'acquire_darks':

        ndarks = widgetData['ndarks']
        logger.info(f'Reading {ndarks} dark spectra')

        # Set the save location for the dark spectra
        save_path = widgetData['rt_save_path'] + '/dark/'

        # Check the output folder exists and generate it if not
        if not os.path.isdir(save_path):
            os.makedirs(save_path)

        # If it exists, create a new folder to avoid overwriting the old one
        else:
            n = 0
            while os.path.isdir(save_path):
                n += 1
                if n == 1:
                    save_path = save_path[:-1] + f'({n})/'
                else:
                    save_path = save_path.split('(')[0] + f'({n})/'
            os.makedirs(save_path)

        # Create an empty array to hold the dark spectra
        dark_arr = np.zeros([ndarks, spectrometer.pixels])

        for i in range(ndarks):
            fname = f'{save_path}spectrum{i:05d}.txt'
            spectrum, info = spectrometer.get_spectrum(fname=fname)
            dark_arr[i] = spectrum[1]

            # Display the spectrum
            spectrum_callback.emit((spectrum, info, True))

            # Update the progress bar
            progress_callback.emit(((i+1)/ndarks)*100)

            # Get if the worker is killed
            if worker.is_killed:
                logger.info('Dark spectra acquisition interupted')

        dark_spec = np.average(dark_arr, axis=0)
        np.savetxt('bin/.dark', dark_spec)

    if acquisition_mode == 'acquire_cont':

        # Set the save location for the dark spectra
        save_path = widgetData['rt_save_path'] + '/spectra/'

        # Check the output folder exists and generate it if not
        if not os.path.isdir(save_path):
            os.makedirs(save_path)

        logger.info('Beginning spectra acquisition')

        # Create a spectrum number parameter
        nspec = 0

        while not worker.is_killed:

            # Check if the loop is paused
            while worker.is_paused:
                time.sleep(0.01)

            # Generate the spectrum filename
            spec_fname = f'{save_path}spectrum_{nspec:05d}.txt'

            # Check the file doesn't already exist
            while os.path.isfile(spec_fname):
                nspec += 1
                spec_fname = f'{save_path}spectrum_{nspec:05d}.txt'

            # Measure the spectrum
            spectrum, info = spectrometer.get_spectrum(fname=spec_fname)

            # Display the spectrum
            spectrum_callback.emit((spectrum, info, True))


# =============================================================================
# Connect to spectrometer
# =============================================================================

def connect_spectrometer(gui):
    """Connect or dissconnect the spectrometer."""
    if not gui.connected_flag:

        # Connect to the spectrometer
        spec = Spectrometer(integration_time=gui.widgets.get("int_time"),
                            coadds=gui.widgets.get("coadds"),
                            correct_dark_counts=gui.widgets.get("nonlin_flag"),
                            correct_nonlinearity=gui.widgets.get("eldark_flag")
                            )

        # Check if connection was successful
        if spec.serial_number is not None:

            # Add the spectrometer to the parent GUI
            gui.spectrometer = spec

            # Update the GUI
            gui.spec_id.setText(gui.spectrometer.serial_number)
            gui.connect_btn.setText('Disconnect')

            # Create a holder for the dark spectra
            gui.dark_spectrum = np.zeros(gui.spectrometer.pixels)

            # Update GUI features
            gui.connected_flag = True
            gui.acquire_test_btn.setEnabled(True)
            gui.update_inttime_btn.setEnabled(True)
            gui.rt_flag_btn.setEnabled(True)
            gui.update_coadds_btn.setEnabled(True)
            gui.acquire_darks_btn.setEnabled(True)
            gui.rt_start_btn.setEnabled(True)
            for k in ["nonlin_flag", "eldark_flag"]:
                gui.widgets[k].setEnabled(False)

    else:
        # Disconnect the spectrometer
        gui.spectrometer.close()

        # Update the GUI
        gui.spec_id.setText('Not connected')
        gui.connect_btn.setText('Connect')

        # Update GUI features
        gui.connected_flag = False
        gui.acquire_test_btn.setEnabled(False)
        gui.update_inttime_btn.setEnabled(False)
        gui.rt_flag_btn.setEnabled(False)
        gui.update_coadds_btn.setEnabled(False)
        gui.acquire_darks_btn.setEnabled(False)
        gui.rt_start_btn.setEnabled(False)
        for k in ["nonlin_flag", "eldark_flag"]:
            gui.widgets[k].setEnabled(True)


# =============================================================================
# Data Buffer
# =============================================================================

class Buffer:
    """A buffer to hold the analysis results for plotting.

    Parameters
    ----------
    size : int
        Size of the buffer (number of rows)
    columns : list
        Column names to use

    Attributes
    ----------
    df : Pandas DataFrame
        Holds the data
    fill : int
        Current fill level of the buffer
    """

    def __init__(self, size, columns):

        self.df = pd.DataFrame(index=np.arange(size), columns=columns)
        self.fill = 0
        self.size = size

    def update(self, data):
        """Add a row to the buffer.

        Parameters
        ----------
        data : list
            The data to add to the buffer. Must have same length and order as
            columns. If full then the oldest row is removed.

        Returns
        -------
        None
        """
        if self.fill < self.size:
            self.df.iloc[self.fill] = data
            self.fill += 1

        else:
            for key in self.df:
                self.df[key] = np.roll(self.df[key], -1)
            self.df.iloc[self.fill-1] = data


# =============================================================================
# Widgets Object
# =============================================================================

class Widgets(dict):
    """Object to allow easy config/info transfer with PyQT Widgets."""

    def __init__(self):
        super().__init__()

    def get(self, key):
        """Get the value of a widget."""
        if type(self[key]) == QTextEdit:
            return self[key].toPlainText()
        elif type(self[key]) == QLineEdit:
            return self[key].text()
        elif type(self[key]) == QComboBox:
            return str(self[key].currentText())
        elif type(self[key]) == QCheckBox:
            return self[key].isChecked()
        elif type(self[key]) in [SpinBox, DSpinBox, QSpinBox, QDoubleSpinBox]:
            return self[key].value()

    def set(self, key, value):
        """Set the value of a widget."""
        if type(self[key]) in [QTextEdit, QLineEdit]:
            self[key].setText(str(value))
        if type(self[key]) == QComboBox:
            index = self[key].findText(value, Qt.MatchFixedString)
            if index >= 0:
                self[key].setCurrentIndex(index)
        if type(self[key]) == QCheckBox:
            self[key].setChecked(value)
        if type(self[key]) in [SpinBox, DSpinBox, QSpinBox, QDoubleSpinBox]:
            self[key].setValue(value)


# =============================================================================
# Spinbox classes
# =============================================================================

# Create a Spinbox object for ease
class DSpinBox(QDoubleSpinBox):
    """Object for generating custom float spinboxes."""

    def __init__(self, value, range=None, step=1.0):
        super().__init__()
        if range is not None:
            self.setRange(*range)
        self.setValue(value)
        self.setSingleStep(step)


class SpinBox(QSpinBox):
    """Object for generating custom integer spinboxes."""

    def __init__(self, value, range):
        super().__init__()
        self.setRange(*range)
        self.setValue(value)


# =============================================================================
# Parameter Table
# =============================================================================

class Table(QTableWidget):
    """Object to build parameter tables."""

    def __init__(self, parent, type, width, height, pname=None):
        super().__init__(parent)

        self._width = width
        self._height = height
        self._type = type
        self._pname = pname

        if self._type == 'param':
            self._param_table()

        if self._type == 'poly':
            self._poly_table()

        header = self.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeToContents)
        # header.setSectionResizeMode(0, QHeaderView.Stretch)

    def _param_table(self):
        """Create a parameter table."""
        self.setFixedWidth(self._width)
        self.setFixedHeight(self._height)
        self.setColumnCount(6)
        self.setRowCount(0)
        self.setHorizontalHeaderLabels(['Name', 'Value', 'Vary?', 'Plume Gas?',
                                        'Xsec Path', ''])

    def _poly_table(self):
        """Create a polynomial table."""
        self.setFixedWidth(self._width)
        self.setFixedHeight(self._height)
        self.setColumnCount(2)
        self.setRowCount(0)
        self.setHorizontalHeaderLabels(['Value', 'Vary?'])

    def add_row(self):
        """Add row to the bottom of the table."""
        n = self.rowCount()
        self.setRowCount(n+1)

        if self._type == 'param':
            cb0 = QCheckBox()
            cb0.setChecked(True)
            self.setCellWidget(n, 2, cb0)
            cb1 = QCheckBox()
            cb1.setChecked(False)
            self.setCellWidget(n, 3, cb1)
            b = QPushButton('Browse')
            self.setItem(n, 1, QTableWidgetItem('0.0'))
            self.setCellWidget(n, 5, b)
            b.clicked.connect(partial(self.set_xsec, n))

        if self._type == 'poly':
            cb = QCheckBox()
            cb.setChecked(True)
            self.setItem(n, 0, QTableWidgetItem('0.0'))
            self.setCellWidget(n, 1, cb)

    def rem_row(self):
        """Remove the last row from the table."""
        rows = [i.row() for i in self.selectedIndexes()]
        for row in sorted(rows, reverse=True):
            self.removeRow(row)

    def contextMenuEvent(self, event):
        """Set up right click to add/remove rows."""
        menu = QMenu(self)
        addAction = menu.addAction('Add')
        remAction = menu.addAction('Remove')
        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == addAction:
            self.add_row()
        if action == remAction:
            self.rem_row()

    def set_xsec(self, n):
        """Command to set the cross-section in the table."""
        cwd = os.getcwd() + '/'
        cwd = cwd.replace("\\", "/")
        fname, _ = QFileDialog.getOpenFileName()
        if cwd in fname:
            fname = fname[len(cwd):]
        self.setItem(n, 4, QTableWidgetItem(fname))

    def setData(self, data):
        """Populate the table using saved config."""
        for i in range(len(data)):
            self.add_row()
            line = data[i]
            if self._type == 'param':
                self.setItem(i, 0, QTableWidgetItem(line[0]))
                self.setItem(i, 1, QTableWidgetItem(str(line[1])))
                self.cellWidget(i, 2).setChecked(line[2])
                self.cellWidget(i, 3).setChecked(line[3])
                self.setItem(i, 4, QTableWidgetItem(line[4]))

            elif self._type == 'poly':
                self.setItem(i, 0, QTableWidgetItem(str(line[0])))
                self.cellWidget(i, 1).setChecked(line[1])

    def getData(self):
        """Extract the information from the table."""
        # Get number of rows
        nrows = self.rowCount()
        data = []

        try:
            # Read the data from a param table
            if self._type == 'param' and nrows > 0:
                for i in range(nrows):
                    row = [self.item(i, 0).text(),
                           float(self.item(i, 1).text()),
                           self.cellWidget(i, 2).isChecked(),
                           self.cellWidget(i, 3).isChecked(),
                           self.item(i, 4).text()]
                    data.append(row)

            # Read the data from a poly table
            elif self._type == 'poly' and nrows > 0:
                for i in range(nrows):

                    row = [float(self.item(i, 0).text()),
                           self.cellWidget(i, 1).isChecked()]
                    data.append(row)
        except AttributeError:
            pass

        return data
