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
                             QPlainTextEdit)

from .parameters import Parameters
from .spectral_analysis import Analyser
from .load_spectra import read_spectrum, average_spectra
from .spectrometers import VSpectrometer


class LogSignals(QObject):
    """Defines signals for the ThreadLogger"""
    signal = pyqtSignal(str)


class ThreadLogger(logging.Handler):

    def __init__(self):
        super().__init__()
        self.log = LogSignals()

    def emit(self, record):
        msg = self.format(record)
        self.log.signal.emit(msg)


class QTextEditLogger(logging.Handler, QObject):
    appendPlainText = pyqtSignal(str)

    def __init__(self, parent):
        super().__init__()
        QObject.__init__(self)
        self.widget = QPlainTextEdit(parent)
        self.widget.setReadOnly(True)
        self.widget.setFont(QFont('Courier', 10))
        self.appendPlainText.connect(self.widget.appendPlainText)

    def emit(self, record):
        msg = self.format(record)
        self.appendPlainText.emit(msg)


# Create a worker signals object to handle worker signals
class WorkerSignals(QObject):
    """Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data
    error
        `tuple` (exctype, value, traceback.format_exc() )
    result
        `object` data returned from processing, anything
    progress
        `int` indicating % progress
    """
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(object)
    progress = pyqtSignal(int)
    plotter = pyqtSignal(list)
    status = pyqtSignal(str)


# Create a worker to handle QThreads
class Worker(QRunnable):
    """Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up

    :param callback: The function callback to run on this worker thread.
                     Supplied args and kwargs will be passed through to the
                     runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    """

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Set killed and paused flags
        self.is_paused = False
        self.is_killed = False

        # Add the callback to our kwargs
        self.kwargs['progress_callback'] = self.signals.progress
        self.kwargs['status_callback'] = self.signals.status
        self.kwargs['plotter'] = self.signals.plotter

    @pyqtSlot()
    def run(self):
        """Initialise the runner function with passed args, kwargs"""

        # Retrieve args/kwargs here; and fire processing using them
        try:
            result = self.fn(self, *self.args, **self.kwargs)
        except Exception:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))

        # Return the result of the processing
        else:
            self.signals.result.emit(result)

        # Done
        finally:
            self.signals.finished.emit()

    def pause(self):
        if self.is_paused:
            self.is_paused = False
        else:
            self.is_paused = True

    def resume(self):
        self.is_paused = False

    def kill(self):
        if self.is_paused:
            self.is_paused = False
        self.is_killed = True


# =============================================================================
# Generate analyser
# =============================================================================

def generate_analyser(widgetData):

    # Pull the parameters from the parameter table
    params = Parameters()
    logging.info('Generating model parameters')

    # Build the parameters from GUI tables
    for line in widgetData['gas_params']:
        params.add(name=line[0], value=line[1], vary=line[2], xpath=line[3])
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
        params.add('k', value=widgetData['k'],
                   vary=widgetData['k_fit'])
        params.add('a_w', value=widgetData['a_w'],
                   vary=widgetData['a_w_fit'])
        params.add('a_k', value=widgetData['a_k'],
                   vary=widgetData['a_k_fit'])

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
                        spike_limit=widgetData['spike_limit'])

    # Report fitting parameters
    logging.info(params.pretty_print(cols=['name', 'value', 'vary']))

    return analyser


# =============================================================================
# Spectra analysis loop
# =============================================================================

def analysis_loop(worker, widgetData, progress_callback, plotter,
                  status_callback):

    # Create a loop counter
    loop = 0

    # Generate the anlyser
    status_callback.emit('Loading')
    analyser = generate_analyser(widgetData)

    # Make a list of column names
    cols = ['File', 'Number', 'Time']
    for par in analyser.params:
        cols += [par, f'{par}_err']
    cols += ['fit_quality', 'int_lo', 'int_hi', 'int_av']

    # Pull the spectra type and file paths
    spec_fnames = widgetData['spec_fnames'].split('\n')
    dark_fnames = widgetData['dark_fnames'].split('\n')
    spec_type = widgetData['spec_type']

    # Pull processing settings
    update_flag = widgetData['update_flag']
    resid_limit = widgetData['resid_limit']
    resid_type = widgetData['resid_type']
    int_limit = [widgetData['lo_int_limit'],
                 widgetData['hi_int_limit']]
    graph_p = [r[0] for r in widgetData['gas_params']]
    interp_meth = widgetData['interp_method']

    # Make a dataframe to hold the fit results
    df = pd.DataFrame(index=np.arange(len(spec_fnames)), columns=cols)

    # Make the results file
    save_path = widgetData['save_path']

    # Read in the dark spectra
    if widgetData['dark_flag']:
        logging.info('Reading dark spectra')
        x, analyser.dark_spec = average_spectra(dark_fnames, spec_type)

    # Set status
    logging.info('Beginning analysis loop')
    worker.signals.status.emit('Analysing')

    # Open the output file
    with open(save_path, 'w') as w:
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
            fname = spec_fnames[loop]

            # Read in the spectrum
            x, y, info, read_err = read_spectrum(fname, spec_type)

            # Fit the spectrum
            fit_result = analyser.fit_spectrum(spectrum=[x, y],
                                               update_params=update_flag,
                                               resid_limit=resid_limit,
                                               resid_type=resid_type,
                                               int_limit=int_limit,
                                               calc_od=graph_p,
                                               interp_method=interp_meth)

            # Add results to the dataframe
            row = [fname, info["spec_no"], info["time"]]
            for par in fit_result.params.values():
                row += [par.fit_val, par.fit_err]
            row += [fit_result.nerr, fit_result.int_lo, fit_result.int_hi,
                    fit_result.int_av]
            df.loc[loop] = row

            # Write the results to file
            w.write(f'{row[0]}')
            for i in row[1:]:
                w.write(f',{i}')
            w.write('\n')

            # Emit the progress
            progress_callback.emit(((loop+1)/len(spec_fnames))*100)

            # Emit graph data
            plotter.emit([fit_result, [x, y], df])

            # Check if analysis is finished
            if loop == len(spec_fnames)-1:
                logging.info('Analysis finished!')
                worker.kill()

            # Update loop counter
            loop += 1


# =============================================================================
# Connect to spectrometer
# =============================================================================

def connect_spectrometer(gui):
    """Connects or dissconnects to the spectrometer"""

    if not gui.connected_flag:

        # Connect to the spectrometer
        gui.spec = VSpectrometer(integration_time=gui.widgets.get("int_time"),
                                 coadds=gui.widgets.get("coadds"))

        # Update the GUI
        gui.spec_id.setText(gui.spec.serial_number)
        gui.connect_btn.setText('Disconnect')

        # Create a holder for the dark spectra
        gui.dark_spectrum = np.zeros(gui.spec.pixels)

        gui.connected_flag = True

    else:
        # Disconnect the spectrometer
        gui.spec.close()

        # Update the GUI
        gui.spec_id.setText('Not connected')
        gui.connect_btn.setText('Connect')

        gui.connected_flag = False


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
    spectrum, info = gui.spec.get_spectrum()

    # Plot the spectrum
    gui.plot_lines[2].setData(*spectrum)


# =============================================================================
# Widgets Object
# =============================================================================

class Widgets(dict):
    """Object to allow easy config/info transfer with PyQT Widgets"""

    def __init__(self):
        super().__init__()

    def get(self, key):
        """Get the value of a widget"""
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
        """Set the value of a widget"""
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
# Spinbox
# =============================================================================

# Create a Spinbox object for ease
class DSpinBox(QDoubleSpinBox):
    """Object for generating custom float spinboxes"""

    def __init__(self, value, range):
        super().__init__()
        self.setRange(*range)
        self.setValue(value)


class SpinBox(QSpinBox):
    """Object for generating custom integer spinboxes"""

    def __init__(self, value, range):
        super().__init__()
        self.setRange(*range)
        self.setValue(value)


# =============================================================================
# Parameter Table
# =============================================================================

class Table(QTableWidget):
    """Object to build parameter tables"""

    def __init__(self, parent, type, width, pname=None):
        super().__init__(parent)

        self._width = width
        self._type = type
        self._pname = pname

        if self._type == 'param':
            self._param_table()

        if self._type == 'poly':
            self._poly_table()

    def _param_table(self):
        self.setFixedWidth(self._width)
        self.setColumnCount(5)
        self.setRowCount(0)
        self.setHorizontalHeaderLabels(['Name', 'Value', 'Vary?', 'Xsec Path',
                                        ''])

    def _poly_table(self):
        self.setFixedWidth(self._width)
        self.setColumnCount(2)
        self.setRowCount(0)
        self.setHorizontalHeaderLabels(['Value', 'Vary?'])

    def add_row(self):
        """Add row to the bottom of the table"""
        n = self.rowCount()
        self.setRowCount(n+1)

        if self._type == 'param':
            cb = QCheckBox()
            self.setCellWidget(n, 2, cb)
            b = QPushButton('Browse')
            self.setItem(n, 1, QTableWidgetItem('0.0'))
            self.setCellWidget(n, 4, b)
            b.clicked.connect(partial(self.set_xsec, n))

        if self._type == 'poly':
            cb = QCheckBox()
            self.setCellWidget(n, 1, cb)

    def rem_row(self):
        """Remove the last row from the table"""
        n = self.rowCount()
        self.setRowCount(n-1)

    def contextMenuEvent(self, event):
        """Set up right click to add/remove rows"""
        menu = QMenu(self)
        addAction = menu.addAction('Add')
        remAction = menu.addAction('Remove')
        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == addAction:
            self.add_row()
        if action == remAction:
            self.rem_row()

    def set_xsec(self, n):
        """Command to set the cross-section in the table"""
        fname = QFileDialog.getOpenFileName()
        self.setItem(n, 3, QTableWidgetItem(fname[0]))

    def setData(self, data):
        """Method to populate the table using saved config"""

        for i in range(len(data)):
            self.add_row()
            line = data[i]
            if self._type == 'param':
                self.setItem(i, 0, QTableWidgetItem(line[0]))
                self.setItem(i, 1, QTableWidgetItem(str(line[1])))
                self.cellWidget(i, 2).setChecked(line[2])
                self.setItem(i, 3, QTableWidgetItem(line[3]))

            elif self._type == 'poly':
                self.setItem(i, 0, QTableWidgetItem(str(line[0])))
                self.cellWidget(i, 1).setChecked(line[1])

    def getData(self):
        """Method to extract the information from the table"""

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
                           self.item(i, 3).text()]
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
