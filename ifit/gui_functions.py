"""Assorted functions for the iFit GUI.

Contains functions for acquiring/analysing spectra in the background and
bespoke widgets for table features in the front end.
"""
import os
import sys
import time
import serial
import logging
import traceback
import numpy as np
import pandas as pd
from datetime import datetime
from functools import partial
import serial.tools.list_ports
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt, QObject, pyqtSignal
from PyQt5.QtWidgets import (QComboBox, QTextEdit, QLineEdit, QDoubleSpinBox,
                             QSpinBox, QCheckBox, QFileDialog, QPushButton,
                             QTableWidgetItem, QMenu, QTableWidget, QDialog,
                             QPlainTextEdit, QHeaderView, QFormLayout, QFrame,
                             QVBoxLayout, QHBoxLayout, QLabel)

from ifit.parameters import Parameters
from ifit.spectral_analysis import Analyser
from ifit.load_spectra import read_spectrum


logger = logging.getLogger(__name__)


# =============================================================================
# Logging text box
# =============================================================================

class QTextEditLogger(logging.Handler, QObject):
    """Record logs to the GUI."""

    appendPlainText = pyqtSignal(str)

    def __init__(self, parent):
        """Initialise."""
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


# =============================================================================
# Scope acquisition worker
# =============================================================================

class AcqScopeWorker(QObject):
    """Handle continuous scope acquisition."""

    # Define Signals
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    plotSpec = pyqtSignal(np.ndarray)

    def __init__(self, spectrometer):
        """Initialise."""
        super(QObject, self).__init__()
        self.spectrometer = spectrometer
        self.is_stopped = False

    def run(self):
        """Launch worker task."""
        while not self.is_stopped:
            spectrum, info = self.spectrometer.get_spectrum()
            self.plotSpec.emit(spectrum)

        self.finished.emit()

    def stop(self):
        """Terminate the acquisition."""
        self.is_stopped = True


# =============================================================================
# Measurement acquisition worker
# =============================================================================

class AcqSpecWorker(QObject):
    """Handle targeted acquisition and saving of spectra."""

    # Define Signals
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    status = pyqtSignal(str)
    error = pyqtSignal(tuple)
    plotSpec = pyqtSignal(np.ndarray)
    setDark = pyqtSignal(np.ndarray)
    setSpec = pyqtSignal(tuple)

    def __init__(self, spectrometer, gps, widgetData):
        """Initialise."""
        super(QObject, self).__init__()
        self.spectrometer = spectrometer
        self.gps = gps
        self.widgetData = widgetData
        self.is_paused = False
        self.is_stopped = False

#   Acquire Dark Spectra ======================================================

    def acquire_dark(self):
        """Acquire dark spectra."""
        try:
            # Log start of acquisition
            logger.info(f'Reading {self.widgetData["ndarks"]} dark spectra')

            # Set the save path location
            save_path = f'{self.widgetData["rt_save_path"]}/dark/'

            # Check the output folder exists and generate it if not
            if not os.path.isdir(save_path):
                os.makedirs(save_path)

            # If it exists create a new folder to avoid overwriting the old one
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
            dark_arr = np.full([self.widgetData["ndarks"],
                                self.spectrometer.pixels],
                               np.nan)

            for i in range(self.widgetData["ndarks"]):

                # Check if the loop is paused
                while self.is_paused:
                    time.sleep(0.01)

                # Read the spectrum and add to the array
                fname = f'{save_path}spectrum_{i:05d}.txt'
                spectrum, info = self.spectrometer.get_spectrum(fname=fname,
                                                                gps=self.gps)
                dark_arr[i] = spectrum[1]

                # Display the spectrum
                self.plotSpec.emit(spectrum)

                # Update the progress bar
                self.progress.emit(((i+1)/self.widgetData["ndarks"])*100)

                # Get if the worker is killed
                if self.is_stopped:
                    logger.info('Dark spectra acquisition interupted')
                    break

            # Record the dark spectrum in the GUI
            dark_spec = np.nanmean(dark_arr, axis=0)
            self.setDark.emit(dark_spec)

        except Exception:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.error.emit((exctype, value, traceback.format_exc()))

        self.finished.emit()

#   Acquire Measurement Spectra ===============================================

    def acquire_spec(self):
        """Acquire measurement spectra."""
        try:
            # Set the save location for the spectra
            save_path = self.widgetData['rt_save_path'] + '/spectra/'

            # Check the output folder exists and generate it if not
            if not os.path.isdir(save_path):
                os.makedirs(save_path)

            logger.info('Beginning spectra acquisition')

            # Create a spectrum number parameter
            nspec = 0

            while not self.is_stopped:

                # Check if the loop is paused
                while self.is_paused:
                    time.sleep(0.01)

                # Generate the spectrum filename
                spec_fname = f'{save_path}spectrum_{nspec:05d}.txt'

                # Check the file doesn't already exist
                while os.path.isfile(spec_fname):
                    nspec += 1
                    if nspec <= 99999:
                        spec_fname = f'{save_path}spectrum_{nspec:05d}.txt'

                    # If max spectrum number is readed (99999) then start a new
                    # folder
                    else:
                        logger.warning('Maximum spectrum number reached, '
                                       + 'creatinga new folder...')
                        count = 1
                        while os.path.isdir(save_path):
                            save_path = f'{self.widgetData["rt_save_path"]}' \
                                        + f'/spectra{count}/'
                            count += 1

                        # Restart the counter
                        nspec = 0
                        spec_fname = f'{save_path}spectrum_{nspec:05d}.txt'

                # Measure the spectrum
                spec, info = self.spectrometer.get_spectrum(fname=spec_fname,
                                                            gps=self.gps)

                # Add the spectrum number to the metadata
                info['spectrum_number'] = nspec

                # Display the spectrum
                self.plotSpec.emit(spec)
                self.setSpec.emit((spec, info, spec_fname))

        except Exception:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.error.emit((exctype, value, traceback.format_exc()))

        self.finished.emit()

#   Controls ==================================================================

    def pause(self):
        """Pause the analysis/acquisition."""
        if self.is_paused:
            self.is_paused = False
        else:
            self.is_paused = True

    def resume(self):
        """Resume the analysis/acquisition."""
        self.is_paused = False

    def stop(self):
        """Terminate the acquisition."""
        self.is_stopped = True


# =============================================================================
# Analysis worker
# =============================================================================

class AnalysisWorker(QObject):
    """Handle analysis of measured spectra."""

    # Define Signals
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    status = pyqtSignal(str)
    error = pyqtSignal(tuple)
    plotData = pyqtSignal(list)
    initializeTable = pyqtSignal(object)

    def __init__(self, mode, widgetData, dark_spec):
        """Initialise."""
        super(QObject, self).__init__()
        self.mode = mode
        self.widgetData = widgetData
        self.dark_spec = dark_spec
        self.is_paused = False
        self.is_stopped = False
        self.spectrum = None

#   Run analysis ==============================================================

    def run(self):
        """Launch the worker function."""
        try:
            self._run()
        except Exception:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.error.emit((exctype, value, traceback.format_exc()))
        self.finished.emit()

    def _run(self):
        """Run analysis."""
        # Create a loop counter
        loop = 0

        # Update the status
        self.status.emit('Loading')

        # Generate the anlyser
        self.analyser = self.generate_analyser()

        # Initialize results table
        self.initializeTable.emit(self.analyser.params)

        # Pull processing settings
        update_flag = self.widgetData['update_flag']
        resid_limit = self.widgetData['resid_limit']
        resid_type = self.widgetData['resid_type']
        int_limit = [self.widgetData['lo_int_limit'],
                     self.widgetData['hi_int_limit']]
        graph_p = [r[0] for r in self.widgetData['gas_params']]
        interp_meth = self.widgetData['interp_method']

        if self.widgetData['preshift_flag']:
            prefit_shift = self.widgetData['prefit_shift']
        else:
            prefit_shift = 0.0

        # Make a list of the fit parameters to form the buffer columns
        buffer_cols = ['Number']
        [buffer_cols.extend([par, f'{par}_err'])
         for par in self.analyser.params]
        buffer_cols += ['fit_quality', 'Lat', 'Lon']

        # Set write mode
        write_mode = 'w'

        # Pull analysis parameters for post analysis
        if self.mode == 'post_analyse':
            spec_fnames = self.widgetData['spec_fnames'].split('\n')
            nspec = len(spec_fnames)
            spec_type = self.widgetData['spec_type']
            wl_calib_file = self.widgetData['wl_calib']
            buffer = Buffer(len(spec_fnames), buffer_cols)
            save_path = self.widgetData['save_path']
            self.status.emit('Analysing')

        # Pull analysis parameters for real time analysis
        elif self.mode == 'rt_analyse':
            spec_fnames = None
            spec_type = 'iFit'
            wl_calib_file = None
            buffer = Buffer(2000, buffer_cols)
            date_str = datetime.strftime(datetime.now(), "%Y-%m-%d_%H%M%S")
            save_path = f'{self.widgetData["rt_save_path"]}/' \
                        + f'{date_str}_iFit_output.csv'
            if os.path.isfile(save_path):
                write_mode = 'a'
            self.status.emit('Acquiring')

        # Set status
        logger.info('Beginning analysis loop')

        # Open the output file
        with open(save_path, write_mode) as w:

            # If writing a new file, add the metadata and columns
            if write_mode == 'w':

                # Write the analysis metaddata
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                w.write('# iFit output file\n'
                        + f'# Analysis Time,{timestamp}\n'
                        + f'# Fit Window Low,{self.widgetData["fit_lo"]}\n'
                        + f'# Fit Window High,{self.widgetData["fit_hi"]}\n'
                        + f'# Flat corr.,{self.widgetData["flat_flag"]}\n'
                        + f'# Dark corr.,{self.widgetData["dark_flag"]}\n'
                        + f'# Stray corr.,{self.widgetData["stray_flag"]}\n')

                # Make a list of column names
                pre_cols = [['outp_lat', 'Lat'],
                            ['outp_lon', 'Lon'],
                            ['outp_alt', 'Alt']]
                post_cols = [['outp_intlo', 'int_lo'],
                             ['outp_inthi', 'int_hi'],
                             ['outp_intav', 'int_av'],
                             ['outp_resmax', 'max_resid'],
                             ['outp_resstd', 'std_resid'],
                             ['outp_fitqual', 'fit_quality']]

                # Add pre columns
                cols = ['File', 'Number', 'Time']
                [cols.append(name) for [key, name] in pre_cols
                 if self.widgetData[key]]

                # Add analysis parameter columns
                for par in self.analyser.params:
                    cols += [par, f'{par}_err']

                # Add post columns
                [cols.append(name) for [key, name] in post_cols
                 if self.widgetData[key]]

                # Write the header
                w.write(cols[0])
                [w.write(f',{c}') for c in cols[1:]]
                w.write('\n')

            # Begin the analysis loop
            while not self.is_stopped:

                # Check if analysis is paused
                while self.is_paused:
                    time.sleep(0.01)

                # Get the spectrum filename for post analysis
                if self.mode == 'post_analyse':
                    # Read in the spectrum
                    fname = spec_fnames[loop]
                    x, y, metadata, read_err = read_spectrum(fname, spec_type,
                                                             wl_calib_file)
                # Pull the spectrum for real time analysis
                elif self.mode == 'rt_analyse':
                    while self.spectrum is None:
                        time.sleep(0.01)
                    [x, y], metadata, fname = self.spectrum
                    self.spectrum = None

                # Fit the spectrum
                fit_result = self.analyser.fit_spectrum(
                    spectrum=[x, y],
                    update_params=update_flag,
                    resid_limit=resid_limit,
                    resid_type=resid_type,
                    int_limit=int_limit,
                    calc_od=graph_p,
                    interp_method=interp_meth,
                    prefit_shift=prefit_shift)

                # Write the fit results to file
                self._write_fit_results(w, fname, metadata, fit_result)

                # Add results to the buffer dataframe
                row = [metadata["spectrum_number"]]
                for par in fit_result.params.values():
                    row += [par.fit_val, par.fit_err]
                row += [fit_result.nerr, metadata['lat'], metadata['lon']]
                buffer.update(row)

                # Emit the progress
                self.progress.emit(((loop+1)/nspec)*100)

                # Emit graph data
                self.plotData.emit([fit_result, [x, y], buffer.df,
                                    metadata["spectrum_number"]])

                # Check if analysis is finished
                if self.mode == 'post_analyse' and loop == nspec-1:
                    logger.info('Analysis finished!')
                    self.stop()

                # Update loop counter
                loop += 1

#   Generate Analyser =========================================================

    def generate_analyser(self):
        """Generate iFit analyser."""
        # Initialise the Parameters
        logger.info('Generating model analyser')
        self.params = Parameters()

        # Pull the parameters from the parameter table
        for line in self.widgetData['gas_params']:
            self.params.add(name=line[0], value=line[1], vary=line[2],
                            xpath=line[4], plume_gas=line[3])
        for i, line in enumerate(self.widgetData['bgpoly_params']):
            self.params.add(name=f'bg_poly{i}', value=line[0], vary=line[1])
        for i, line in enumerate(self.widgetData['offset_params']):
            self.params.add(name=f'offset{i}', value=line[0], vary=line[1])
        for i, line in enumerate(self.widgetData['shift_params']):
            self.params.add(name=f'shift{i}', value=line[0], vary=line[1])

        # Check if ILS is in the fit
        if self.widgetData['ils_mode'] == 'Manual':
            self.params.add('fwem', value=float(self.widgetData['fwem']),
                            vary=self.widgetData['fwem_fit'])
            self.params.add('k', value=float(self.widgetData['k']),
                            vary=self.widgetData['k_fit'])
            self.params.add('a_w', value=float(self.widgetData['a_w']),
                            vary=self.widgetData['a_w_fit'])
            self.params.add('a_k', value=float(self.widgetData['a_k']),
                            vary=self.widgetData['a_k_fit'])

        # Add the light dilution factor
        if self.widgetData['ldf_fit'] or self.widgetData['ldf'] != 0.0:
            self.params.add('LDF', value=float(self.widgetData['ldf']),
                            vary=self.widgetData['ldf_fit'])

        # Report fitting parameters
        logger.info(self.params.pretty_print(cols=['name', 'value', 'vary',
                                                   'xpath']))

        # Get the bad pixels
        if self.widgetData['bad_pixels'] != '':
            bad_pixels = [int(i) for i
                          in self.widgetData['bad_pixels'].split(',')]
        else:
            bad_pixels = []

        # Generate the analyser
        analyser = Analyser(params=self.params,
                            fit_window=[self.widgetData['fit_lo'],
                                        self.widgetData['fit_hi']],
                            frs_path=self.widgetData['frs_path'],
                            model_padding=self.widgetData['model_padding'],
                            model_spacing=self.widgetData['model_spacing'],
                            flat_flag=self.widgetData['flat_flag'],
                            flat_path=self.widgetData['flat_path'],
                            stray_flag=self.widgetData['stray_flag'],
                            stray_window=[self.widgetData['stray_lo'],
                                          self.widgetData['stray_hi']],
                            dark_flag=self.widgetData['dark_flag'],
                            ils_type=self.widgetData['ils_mode'],
                            ils_path=self.widgetData['ils_path'],
                            despike_flag=self.widgetData['despike_flag'],
                            spike_limit=self.widgetData['spike_limit'],
                            bad_pixels=bad_pixels)

        # Add the dark spectrum
        if self.dark_spec is not None:
            analyser.dark_spec = self.dark_spec
        else:
            logger.warning('No dark spectrum, deactivating dark correction')
            analyser.dark_flag = False

        return analyser

#   Write results to file =====================================================

    def _write_fit_results(self, writer, spec_fname, metadata, fit_result):
        # Write spectrum info to file
        writer.write(f'{spec_fname},{metadata["spectrum_number"]},'
                     + f'{metadata["timestamp"]}')

        if self.widgetData['outp_lat']:
            writer.write(f',{metadata["lat"]}')
        if self.widgetData['outp_lon']:
            writer.write(f',{metadata["lon"]}')
        if self.widgetData['outp_alt']:
            writer.write(f',{metadata["alt"]}')

        # Write fit results to file
        for par in fit_result.params.values():
            writer.write(f',{par.fit_val},{par.fit_err}')

        # Write fit quality info and start a new line
        if self.widgetData['outp_intlo']:
            writer.write(f',{fit_result.int_lo}')
        if self.widgetData['outp_inthi']:
            writer.write(f',{fit_result.int_hi}')
        if self.widgetData['outp_intav']:
            writer.write(f',{fit_result.int_av}')
        if self.widgetData['outp_resmax']:
            writer.write(f',{fit_result.resid_max}')
        if self.widgetData['outp_resstd']:
            writer.write(f',{fit_result.resid_std}')
        if self.widgetData['outp_fitqual']:
            writer.write(f',{fit_result.nerr}')
        writer.write('\n')

#   Controls ==================================================================

    def pause(self):
        """Pause the analysis."""
        if self.is_paused:
            self.is_paused = False
        else:
            self.is_paused = True

    def resume(self):
        """Resume the analysis."""
        self.is_paused = False

    def stop(self):
        """Terminate the analysis."""
        self.is_stopped = True

    def set_spectrum(self, spectrum):
        """Set the spectrum file to be analysed in real time analysis."""
        self.spectrum = spectrum


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
        """Initialise."""
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


# ============================================================================
# GPS Worker
# ============================================================================

class GPSWorker(QObject):
    """Handle real-time GPS updates."""

    # Define signals
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    position = pyqtSignal(list)

    def __init__(self, gps):
        """Initialise."""
        super(QObject, self).__init__()
        self.gps = gps
        self.is_stopped = False

    def run(self):
        """Launch worker task."""
        while not self.is_stopped:
            ts = self.gps.timestamp
            lat = self.gps.lat
            lon = self.gps.lon
            alt = self.gps.alt
            self.position.emit([ts, lat, lon, alt])
            time.sleep(0.5)

        self.finished.emit()

    def stop(self):
        """Terminate the worker."""
        self.is_stopped = True


# =============================================================================
# Widgets Object
# =============================================================================

class Widgets(dict):
    """Object to allow easy config/info transfer with PyQT Widgets."""

    def __init__(self):
        """Initialise."""
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
        """Initialise."""
        super().__init__()
        if range is not None:
            self.setRange(*range)
        self.setValue(value)
        self.setSingleStep(step)


class SpinBox(QSpinBox):
    """Object for generating custom integer spinboxes."""

    def __init__(self, value, range):
        """Initialise."""
        super().__init__()
        self.setRange(*range)
        self.setValue(value)


# =============================================================================
# Parameter Table
# =============================================================================

class ParamTable(QTableWidget):
    """Object to build parameter tables."""

    def __init__(self, parent, type, width, height, pname=None):
        """Initialise."""
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

        return data


# =============================================================================
# GPS Combo Box
# =============================================================================

class GPSComboBox(QComboBox):
    def __init__(self, *args, **kwargs):
        super(GPSComboBox, self).__init__()
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.showMenu)

    def showMenu(self, pos):
        menu = QMenu()
        update_action = menu.addAction("Update")
        action = menu.exec_(self.mapToGlobal(pos))
        if action == update_action:
            self.update_com_list()

    def update_com_list(self):
        """Update available COM ports."""
        ports = serial.tools.list_ports.comports()
        self.clear()
        self.addItems([port.device for port in ports])


def browse(gui, widget, mode='single', filter=None):
    """Open native file dialogue."""
    # Check if specified file extensions
    if filter is not None:
        filter = filter + ';;All Files (*)'

    # Pick a single file to read
    if mode == 'single':
        fname, _ = QFileDialog.getOpenFileName(gui, 'Select File', '', filter)

    elif mode == 'multi':
        fname, _ = QFileDialog.getOpenFileNames(gui, 'Select Files', '',
                                                filter)

    elif mode == 'save':
        fname, _ = QFileDialog.getSaveFileName(gui, 'Save As', '', filter)

    elif mode == 'folder':
        fname = QFileDialog.getExistingDirectory(gui, 'Select Folder')

    # Get current working directory
    cwd = os.getcwd() + '/'
    cwd = cwd.replace("\\", "/")

    # Update the relavant widget for a single file
    if type(fname) == str and fname != '':
        if cwd in fname:
            fname = fname[len(cwd):]
        widget.setText(fname)

    # And for multiple files
    elif type(fname) == list and fname != []:
        for i, f in enumerate(fname):
            if cwd in f:
                fname[i] = f[len(cwd):]
        widget.setText('\n'.join(fname))


class GPSWizard(QDialog):
    """Opens a wizard to define a new station."""

    def __init__(self, parent=None):
        """Initialise the window."""
        super(GPSWizard, self).__init__(parent)

        self.parent = parent

        # Set the window properties
        self.setWindowTitle('Setup GPS')

        self._createApp()

    def _createApp(self):
        # Set the layout
        layout = QVBoxLayout()

        # Setup entry widgets
        self.widgets = {'COM Port': GPSComboBox(),
                        'Baudrate': QLineEdit('4800'),
                        'Stream to File': QCheckBox()}
        form_layout = QFormLayout()
        for key, item in self.widgets.items():
            form_layout.addRow(key + ':', item)

        # Add current COM ports and set file to true
        self.widgets['COM Port'].update_com_list()
        self.widgets['Stream to File'].setChecked(True)

        layout.addLayout(form_layout)

        # Add file name in its own layout
        fn_layout = QHBoxLayout()
        fname = f'{self.parent.widgets.get("rt_save_path")}/gps_output.txt'
        self.widgets['Filename'] = QLineEdit(fname)
        btn = QPushButton('Browse')
        # btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.widgets['Filename'],
                                    'save', None))
        fn_layout.addWidget(QLabel('GPS\nFilename:'))
        fn_layout.addWidget(self.widgets['Filename'])
        fn_layout.addWidget(btn)
        layout.addLayout(fn_layout)

        # Add cancel and accept buttons
        btn_layout = QHBoxLayout()
        cancel_btn = QPushButton('Cancel')
        cancel_btn.clicked.connect(self.cancel_action)
        accept_btn = QPushButton('Accept')
        accept_btn.clicked.connect(self.accept_action)
        btn_layout.addWidget(accept_btn)
        btn_layout.addWidget(cancel_btn)

        layout.addLayout(btn_layout)

        self.setLayout(layout)

    def accept_action(self):
        """Record the GPS connection data and exit."""
        self.gps_kwargs = {'comport': self.widgets['COM Port'].currentText(),
                           'baudrate': self.widgets['Baudrate'].text()}

        if self.widgets['Stream to File'].isChecked():
            self.gps_kwargs['filename'] = self.widgets['Filename'].text()
        else:
            self.gps_kwargs['filename'] = None
        self.accept()

    def cancel_action(self):
        """Close the window without connecting to the GPS."""
        self.close()


class QHLine(QFrame):
    """Horizontal line widget."""

    def __init__(self):
        """Initialise."""
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)


class QVLine(QFrame):
    """Vertical line widget."""

    def __init__(self):
        """Initialise."""
        super(QVLine, self).__init__()
        self.setFrameShape(QFrame.VLine)
        self.setFrameShadow(QFrame.Sunken)
