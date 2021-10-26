"""Contains scripts to launch GUI tools.

These functions laucnch new windows for tasks outside the standard GUI,
including:
    - Calculate light dilution
    - Characterise the flat spectrum and instrument line shape
    - Measure fluxes from traverse measurements
"""
import os
import sys
import yaml
import logging
import traceback
import numpy as np
import pandas as pd
import pyqtgraph as pg
from functools import partial
from scipy.optimize import curve_fit
from scipy.interpolate import griddata
from PyQt5.QtGui import QIcon, QPalette, QColor
from PyQt5.QtCore import Qt, QThreadPool, QObject, pyqtSignal, QThread
from PyQt5.QtWidgets import (QMainWindow, QWidget, QApplication, QGridLayout,
                             QLabel, QTextEdit, QLineEdit, QPushButton, QFrame,
                             QFileDialog, QScrollArea, QCheckBox, QSplitter,
                             QComboBox, QDoubleSpinBox, QTableWidget,
                             QTableWidgetItem, QTabWidget, QMessageBox)

from ifit.gui_functions import QTextEditLogger, DSpinBox, Widgets
from ifit.parameters import Parameters
from ifit.spectral_analysis import Analyser
from ifit.load_spectra import average_spectra
from ifit.light_dilution import generate_ld_curves

try:
    from .make_ils import super_gaussian
    from .haversine import haversine
except ImportError:
    from make_ils import super_gaussian
    from haversine import haversine


COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
          '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

logger = logging.getLogger(__name__)


# =============================================================================
# Calculate Flux
# =============================================================================

class CalcFlux(QMainWindow):
    """Open a window for flux analysis."""

    def __init__(self, parent=None):
        """Initialise the window."""
        super(CalcFlux, self).__init__(parent)

        # Set the window properties
        self.setWindowTitle('Calculate Flux')
        self.statusBar().showMessage('Ready')
        self.setGeometry(40, 40, 1000, 500)
        self.setWindowIcon(QIcon('bin/icons/flux.png'))

        # Set the window layout
        self.generalLayout = QGridLayout()
        self._centralWidget = QScrollArea()
        self.widget = QWidget()
        self.setCentralWidget(self._centralWidget)
        self.widget.setLayout(self.generalLayout)

        # Scroll Area Properties
        self._centralWidget.setWidgetResizable(True)
        self._centralWidget.setWidget(self.widget)

        self.save_flag = False

        self.widgets = Widgets()

        self._createApp()
        self.load_config()

    def _createApp(self):
        """Create the app widgets."""
        # Create a frame to hold program controls
        self.controlFrame = QFrame(self)
        self.controlFrame.setFrameShape(QFrame.StyledPanel)

        # Create a frame to hold volcano data
        self.volcanoFrame = QFrame(self)
        self.volcanoFrame.setFrameShape(QFrame.StyledPanel)

        # Create a frame to hold program outputs
        self.outputFrame = QFrame(self)
        self.outputFrame.setFrameShape(QFrame.StyledPanel)

        # Create a frame to hold graphs
        self.graphFrame = QFrame(self)
        self.graphFrame.setFrameShape(QFrame.StyledPanel)

        # Add splitters to allow for adjustment
        splitter1 = QSplitter(Qt.Vertical)
        splitter1.addWidget(self.controlFrame)
        splitter1.addWidget(self.volcanoFrame)
        splitter1.addWidget(self.outputFrame)

        splitter2 = QSplitter(Qt.Horizontal)
        splitter2.addWidget(splitter1)
        splitter2.addWidget(self.graphFrame)

        # Pack the Frames and splitters
        self.generalLayout.addWidget(splitter2)

        self._createControls()
        self._createVolcano()
        self._createOutput()
        self._createGraphs()

# =============================================================================
#   Program Controls
# =============================================================================

    def _createControls(self):
        """Build the main GUI controls."""
        # Generate main layout
        layout = QGridLayout(self.controlFrame)

        # Create input for iFit output
        layout.addWidget(QLabel('iFit File:'), 0, 0)
        self.widgets['so2_path'] = QLineEdit()
        self.widgets['so2_path'].setFixedSize(200, 25)
        layout.addWidget(self.widgets['so2_path'], 0, 1)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.widgets['so2_path'],
                                    'single', "Comma Separated (*.csv)"))
        layout.addWidget(btn, 0, 2)

        # Add checkbox to remove bad fits
        self.widgets['despike_flag'] = QCheckBox('Remove\nBad fits?')
        layout.addWidget(self.widgets['despike_flag'], 0, 3)

        # Create input for GPS intput
        layout.addWidget(QLabel('GPS File:'), 1, 0)
        self.widgets['gps_path'] = QLineEdit()
        self.widgets['gps_path'].setFixedSize(200, 25)
        layout.addWidget(self.widgets['gps_path'], 1, 1)
        gps_btn = QPushButton('Browse')
        gps_btn.setFixedSize(70, 25)
        gps_btn.clicked.connect(partial(browse, self, self.widgets['gps_path'],
                                        'single', "GPS File (*.txt)"))
        layout.addWidget(gps_btn, 1, 2)

        # Add checkbox to remove bad fits
        def gps_toggle():
            state = self.widgets['gps_file_flag'].isChecked()
            self.widgets['gps_path'].setReadOnly(state)
            gps_btn.setEnabled(state)
        self.widgets['gps_file_flag'] = QCheckBox('Use Separate\nGPS File?')
        self.widgets['gps_file_flag'].stateChanged.connect(gps_toggle)
        layout.addWidget(self.widgets['gps_file_flag'], 1, 3)
        gps_toggle()

        # Create inputs for the time difference
        layout.addWidget(QLabel('GPS Time\nDifference:'), 2, 0)
        self.widgets['tdiff'] = QDoubleSpinBox()
        self.widgets['tdiff'].setRange(-24, 24)
        self.widgets['tdiff'].setValue(0.0)
        layout.addWidget(self.widgets['tdiff'], 2, 1)

        # Make a button to read in the data
        btn = QPushButton('Import')
        btn.setFixedSize(100, 25)
        btn.clicked.connect(self.import_data)
        layout.addWidget(btn, 3, 3)

        # Create input for output
        layout.addWidget(QLabel('Output\nFolder:'), 3, 0)
        self.widgets['out_path'] = QLineEdit()
        self.widgets['out_path'].setFixedSize(200, 25)
        layout.addWidget(self.widgets['out_path'], 3, 1)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.widgets['out_path'],
                                    'folder'))
        layout.addWidget(btn, 3, 2)

# =============================================================================
#   Volcano Data
# =============================================================================

    def _createVolcano(self):
        """Create the volcano data inputs."""
        # Load volcano data
        self.volcano_data = {}
        if os.path.isfile('bin/volcano_data.yml'):
            with open('bin/volcano_data.yml', 'r') as ymlfile:
                self.volcano_data = yaml.load(ymlfile, Loader=yaml.FullLoader)

        # Create the layout
        layout = QGridLayout(self.volcanoFrame)

        # Create a combobox to hold the pre-saved volcano data
        layout.addWidget(QLabel('Volcano:'), 0, 0)
        self.widgets['volcano'] = QComboBox()
        self.widgets['volcano'].addItems(
            ['--select--'] + list(self.volcano_data.keys()))
        self.widgets['volcano'].currentIndexChanged.connect(
            self.update_volcano_data)
        layout.addWidget(self.widgets['volcano'], 0, 1)

        # Create inputs for the volcano latitude
        layout.addWidget(QLabel('Volcano\nLatitude:'), 1, 0)
        self.vlat = QLineEdit()
        self.vlat.textChanged.connect(self._on_volcano_change)
        layout.addWidget(self.vlat, 1, 1)

        # Create inputs for the volcano longitude
        layout.addWidget(QLabel('Volcano\nLongitutde:'), 2, 0)
        self.vlon = QLineEdit()
        self.vlon.textChanged.connect(self._on_volcano_change)
        layout.addWidget(self.vlon, 2, 1)

        # Create inputs for the wind speed
        layout.addWidget(QLabel('Wind\nSpeed:'), 1, 2)
        self.widgets['wind_speed'] = QDoubleSpinBox()
        self.widgets['wind_speed'].setRange(0, 100)
        self.widgets['wind_speed'].setValue(1.0)
        layout.addWidget(self.widgets['wind_speed'], 1, 3)

        # Create input for wind units
        self.widgets['wind_units'] = QComboBox()
        self.widgets['wind_units'].addItems(['m/s', 'knots'])
        layout.addWidget(self.widgets['wind_units'], 1, 4)

        # Create inputs for the wind speed
        layout.addWidget(QLabel('Wind\nError:'), 2, 2)
        self.widgets['wind_error'] = QDoubleSpinBox()
        self.widgets['wind_error'].setRange(0, 1000)
        layout.addWidget(self.widgets['wind_error'], 2, 3)

        # Create input for wind units
        self.widgets['err_units'] = QComboBox()
        self.widgets['err_units'].addItems(['%', 'abs'])
        layout.addWidget(self.widgets['err_units'], 2, 4)

# =============================================================================
#   Create outputs
# =============================================================================

    def _createOutput(self):
        """Create the program outputs."""
        # Generate the layout
        layout = QGridLayout(self.outputFrame)

        # Add a button to calculate the flux
        calc_btn = QPushButton('Calculate Flux')
        calc_btn.setFixedSize(90, 25)
        calc_btn.clicked.connect(self.calc_flux)
        layout.addWidget(calc_btn, 0, 0)

        # Add a button to remove the last flux
        rem_btn = QPushButton('Remove Last')
        rem_btn.setFixedSize(90, 25)
        rem_btn.clicked.connect(self.del_trav)
        layout.addWidget(rem_btn, 0, 1)

        # Add a button to save the fluxes
        sav_btn = QPushButton('Save Fluxes')
        sav_btn.setFixedSize(90, 25)
        sav_btn.clicked.connect(self.save_fluxes)
        layout.addWidget(sav_btn, 0, 2)

        # Create a combobox to hold the pre-saved volcano data
        label = QLabel('Flux\nUnits:')
        label.setAlignment(Qt.AlignRight)
        layout.addWidget(label, 0, 3)
        self.widgets['flux_units'] = QComboBox()
        self.widgets['flux_units'].addItems(['kg/s', 't/day'])
        layout.addWidget(self.widgets['flux_units'], 0, 4)

        # Add a table to hold the flux results
        self.fluxTable = QTableWidget()
        self.fluxTable.setColumnCount(2)
        self.fluxTable.setRowCount(0)
        self.fluxTable.setHorizontalHeaderLabels(['Flux', 'Error'])
        layout.addWidget(self.fluxTable, 1, 0, 1, 5)

        # Create a textbox to display the program messages
        self.logBox = QTextEditLogger(self)
        self.logBox.setFormatter(logging.Formatter('%(message)s'))
        logging.getLogger().addHandler(self.logBox)
        logging.getLogger().setLevel(logging.INFO)
        layout.addWidget(self.logBox.widget, 2, 0, 1, 5)

# =============================================================================
#   Graphs
# =============================================================================

    def _createGraphs(self):
        """Generate the graphs."""
        layout = QGridLayout(self.graphFrame)
        pg.setConfigOptions(antialias=True)

        # Generate tabs for the graphs and settings
        tab1 = QWidget()
        tab2 = QWidget()

        # Form the tab widget
        tabwidget = QTabWidget()
        tabwidget.addTab(tab1, 'Traverse')
        tabwidget.addTab(tab2, 'Map')
        layout.addWidget(tabwidget, 0, 0)

        g1layout = QGridLayout(tab1)
        graphwin = pg.GraphicsLayoutWidget(show=True)

        # Make the graphs
        self.ax = graphwin.addPlot()

        self.ax.setDownsampling(mode='peak')
        self.ax.setClipToView(True)
        self.ax.showGrid(x=True, y=True)

        # Add axis labels
        self.ax.setLabel('left', 'SO<sub>2</sub> SCD (molec cm<sup>-2</sup>)')
        self.ax.setLabel('bottom', 'Time')

        # Add the plot elements
        b0 = QColor('#1f77b4')
        b0.setAlpha(120)
        l1 = pg.PlotCurveItem(pen=None)
        l2 = pg.PlotCurveItem(pen=None)
        pfill = pg.FillBetweenItem(l1, l2, pen=None, brush=b0)
        self.ax.addItem(pfill)
        line = pg.PlotCurveItem(pen=pg.mkPen(COLORS[0], width=2))
        sline = pg.PlotCurveItem(pen=pg.mkPen(COLORS[2], width=2))

        # Add the region selector
        region_select = pg.LinearRegionItem()
        region_select.setZValue(-1)
        region_select.sigRegionChangeFinished.connect(self._on_region_change)

        # Add the baseline
        baseline = pg.InfiniteLine(0, angle=0, movable=True)

        self.ax_elements = {'line': line, 'sline': sline, 'hi_err': l1,
                            'lo_err': l2, 'region_select': region_select,
                            'baseline': baseline}
        for el in self.ax_elements.values():
            self.ax.addItem(el)

        # Add the graphs to the layout
        g1layout.addWidget(graphwin, 0, 0, 0, 5)

        # Create a combobox to determine x axis as time or number
        self.widgets['x_axis'] = QComboBox()
        self.widgets['x_axis'].addItems(['Time', 'Number'])
        self.widgets['x_axis'].currentIndexChanged.connect(self.switch_x_axis)
        g1layout.addWidget(self.widgets['x_axis'], 1, 2)

        g2layout = QGridLayout(tab2)
        graphwin = pg.GraphicsLayoutWidget(show=True)

        # Make the maps
        self.mapax = graphwin.addPlot()
        self.mapax.setDownsampling(mode='peak')
        self.mapax.setClipToView(True)
        self.mapax.showGrid(x=True, y=True)
        self.mapax.setAspectLocked(True)
        self.mapax.addLegend()

        # Add plot elements
        full_trav_line = pg.PlotCurveItem(pen=pg.mkPen(COLORS[0], width=2))
        trav_scat = pg.ScatterPlotItem()
        trav_scat.setZValue(-1)
        trav_line = pg.PlotCurveItem(pen=pg.mkPen(COLORS[2], width=2))
        volcano = pg.ScatterPlotItem(symbol='o', size=10, name='Volcano',
                                     brush=pg.mkBrush(color='r'))
        plume_start = pg.ScatterPlotItem(
            symbol='t1', pen=pg.mkPen(color='w'), size=10,
            brush=pg.mkBrush(color=COLORS[2]), name='Plume Start')
        plume_stop = pg.ScatterPlotItem(
            symbol='t', pen=pg.mkPen(color='w'), size=10,
            brush=pg.mkBrush(color=COLORS[2]), name='Plume Centre')
        plume_cent = pg.ScatterPlotItem(
            symbol='o', pen=pg.mkPen(color='w'), size=10,
            brush=pg.mkBrush(color=COLORS[2]), name='Plume End')
        self.map_elements = {'full_trav_line': full_trav_line,
                             'trav_line': trav_line,
                             'trav_scat': trav_scat,
                             'plume_start': plume_start,
                             'plume_stop': plume_stop,
                             'plume_cent': plume_cent,
                             'volcano': volcano}
        for el in self.map_elements.values():
            self.mapax.addItem(el)

        # Add axis labels
        self.mapax.setLabel('left', 'Latitude')
        self.mapax.setLabel('bottom', 'Longitude')

        # Generate the colormap
        self.cm = pg.colormap.get('magma', source='matplotlib')

        # Make pens for plotting
        self.p0 = pg.mkPen(color='#1f77b4', width=1.5)
        self.p1 = pg.mkPen(color='#ff7f0e', width=2.0)
        self.b0 = QColor('#1f77b4')
        self.b0.setAlpha(120)

        # Add the graphs to the layout
        g2layout.addWidget(graphwin, 0, 0, 0, 0)

# =============================================================================
#   Update Volcano Data
# =============================================================================

    def update_volcano_data(self):
        """Update the volcano data on combobox change."""
        volc = str(self.widgets.get('volcano'))

        if volc != '--select--':
            data = self.volcano_data[volc]
            self.vlat.setText(str(data[0]))
            self.vlon.setText(str(data[1]))

# =============================================================================
#   Slot functions
# =============================================================================

    def _on_region_change(self):
        try:
            # Find the bounds of the selected area
            i0, i1 = self.ax_elements['region_select'].getRegion()
            if self.widgets.get('x_axis') == 'Time':
                idx = np.where(np.logical_and(self.so2_time >= i0,
                                              self.so2_time <= i1))
            else:
                idx = np.where(np.logical_and(self.so2_num >= i0,
                                              self.so2_num <= i1))

            # Extract the relevant gps data and interpolate onto the so2 grid
            so2_time = self.so2_time[idx]
            so2_num = self.so2_num[idx]
            lat = self.lat[idx]
            lon = self.lon[idx]

            # Update the graph
            ploty = self.so2_scd[idx] / 10**self.order
            if self.widgets.get('x_axis') == 'Time':
                plotx = so2_time
            else:
                plotx = so2_num
            self.ax_elements['sline'].setData(plotx, ploty)

            # Update the map
            self.map_elements['trav_line'].setData(lon, lat)
            self.map_elements['plume_start'].setData([lon[0]], [lat[0]])
            self.map_elements['plume_stop'].setData([lon[-1]], [lat[-1]])

        except AttributeError:
            pass

    def _on_volcano_change(self):
        try:
            vlat = float(self.vlat.text())
            vlon = float(self.vlon.text())
            self.map_elements['volcano'].setData([vlon], [vlat])
        except ValueError:
            pass

# =============================================================================
#   Import data
# =============================================================================

    def import_data(self):
        """Import the traverse data."""
        logger.info('Importing traverse data...')

        # Ask to save any outstanding fluxes
        if self.save_flag:
            options = QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
            reply = QMessageBox.question(self, 'Message',
                                         "Would you like to save the fluxes?",
                                         options, QMessageBox.Cancel)

            if reply == QMessageBox.Yes:
                self.save_fluxes()
            elif reply == QMessageBox.No:
                pass
            else:
                return

        # Read in the SO2 results
        logger.info('Importing gas data...')
        so2_df = pd.read_csv(self.widgets.get('so2_path'),
                             parse_dates=['Time'], comment='#')

        self.so2_time = np.array([t.hour + t.minute/60 + t.second/3600
                                  for t in so2_df['Time']])
        self.so2_num = so2_df['Number'].to_numpy()
        self.so2_scd = so2_df['SO2'].to_numpy()
        self.so2_err = so2_df['SO2_err'].to_numpy()

        if self.widgets['despike_flag'].isChecked():
            try:
                mask = so2_df['fit_quality'] != 1
                self.so2_scd = np.ma.masked_where(mask, self.so2_scd)
                self.so2_err = np.ma.masked_where(mask, self.so2_err)
            except KeyError:
                logger.warning('Fit quality not found in iFit file!')

        if self.widgets['gps_file_flag'].isChecked():

            # Read in the GPS data
            logger.info('Importing GPS data...')
            gps_df = pd.read_table(self.widgets.get('gps_path'), sep='\t',
                                   parse_dates=['time'])
            lat = gps_df['latitude'].to_numpy()
            lon = gps_df['longitude'].to_numpy()

            # Correct the GPS time for the time difference
            tdiff = float(self.widgets.get('tdiff'))
            gps_time = np.array([t.hour + t.minute/60 + t.second/3600 + tdiff
                                 for t in gps_df['time']])

            # Interpolate the GPS locations onto the spectra times
            self.lat = griddata(gps_time, lat, self.so2_time)
            self.lon = griddata(gps_time, lon, self.so2_time)

        else:
            # Pull from the iFit results file
            self.lat = so2_df['Lat'].to_numpy()
            self.lon = so2_df['Lon'].to_numpy()

        ploty = self.so2_scd

        if np.nanmax(self.so2_scd) > 1e6:
            self.order = int(np.ceil(np.log10(np.nanmax(ploty)))) - 1
            ploty = ploty / 10**self.order
            plot_err = [ploty - self.so2_err/10**self.order,
                        ploty + self.so2_err/10**self.order]
            self.ax.setLabel('left', f'SO2 SCD (1e{self.order})')

        if self.widgets.get('x_axis') == 'Time':
            plotx = self.so2_time
        else:
            plotx = self.so2_num

        # Update graph
        self.ax_elements['line'].setData(plotx, ploty)
        self.ax_elements['hi_err'].setData(plotx, plot_err[0])
        self.ax_elements['lo_err'].setData(plotx, plot_err[1])
        self.ax_elements['region_select'].setRegion([plotx[0], plotx[-1]])
        self.ax_elements['baseline'].setValue(0)

        # Update map
        norm = (ploty-np.nanmin(ploty)) / np.nanmax(ploty-np.nanmin(ploty))
        pens = [pg.mkPen(color=self.cm.map(val)) for val in norm]
        brushes = [pg.mkBrush(color=self.cm.map(val))
                   for val in norm]
        self.map_elements['trav_scat'].setData(
            x=self.lon, y=self.lat, pen=pens, brush=brushes)
        self.map_elements['full_trav_line'].setData(x=self.lon, y=self.lat)

        # Create a traverse counter and dictionary to hold the results
        self.flux_data = {}
        self.trav_no = 0

        # Reset the results table
        self.fluxTable.setRowCount(0)

        # Turn off the svae flag
        self.save_flag = False

        logger.info('Traverses imported!')

    def switch_x_axis(self):
        """Switch traverse x axis."""
        try:
            ploty = self.so2_scd

            if np.nanmax(self.so2_scd) > 1e6:
                self.order = int(np.ceil(np.log10(np.nanmax(ploty)))) - 1
                ploty = ploty / 10**self.order
                plot_err = [ploty - self.so2_err/10**self.order,
                            ploty + self.so2_err/10**self.order]
                self.ax.setLabel('left', f'Fit value (1e{self.order})')

            if self.widgets.get('x_axis') == 'Time':
                plotx = self.so2_time
            else:
                plotx = self.so2_num

            # Update the graph
            self.ax_elements['line'].setData(plotx, ploty)
            self.ax_elements['hi_err'].setData(plotx, plot_err[0])
            self.ax_elements['lo_err'].setData(plotx, plot_err[1])
            self.ax_elements['region_select'].setRegion([plotx[0], plotx[-1]])

        except AttributeError:
            pass

# =============================================================================
# Calculate Flux
# =============================================================================

    def calc_flux(self):
        """Calculate the flux from the selected traverse."""
        logger.info('Calculating flux:')

        # Pull the relavant data from the GUI
        vlat = float(self.vlat.text())
        vlon = float(self.vlon.text())
        wind_speed = float(self.widgets.get('wind_speed'))
        wind_error = float(self.widgets.get('wind_error'))
        flux_units = self.widgets.get('flux_units')
        wind_units = self.widgets.get('wind_units')
        err_units = self.widgets.get('err_units')

        # Change units if required
        if err_units == 'abs':
            wind_error = wind_error / wind_speed
        else:
            wind_error = wind_error / 100

        # If wind speed is in knotts, convert to m/s
        if wind_units == 'knots':
            wind_speed = wind_speed * 0.5144444444
        else:
            wind_speed = wind_speed

        logger.info(f'Wind speed: {wind_speed:.02f} m/s')

        # Find the bounds of the selected area
        i0, i1 = self.ax_elements['region_select'].getRegion()
        if self.widgets.get('x_axis') == 'Time':
            idx = np.where(np.logical_and(self.so2_time >= i0,
                                          self.so2_time <= i1))
        else:
            idx = np.where(np.logical_and(self.so2_num >= i0,
                                          self.so2_num <= i1))

        # Extract the travese from the SO2 data
        so2_time = self.so2_time[idx]
        so2_scd = self.so2_scd[idx]
        so2_err = self.so2_err[idx]
        lat = self.lat[idx]
        lon = self.lon[idx]

        # Correct the baseline in SO2
        baseline = self.ax_elements['baseline'].value()*(10**self.order)
        so2_scd = np.subtract(so2_scd, baseline)

        # Find the centre of mass of the plume
        cum_so2_scd = np.nancumsum(so2_scd)
        peak_idx = np.abs(cum_so2_scd - cum_so2_scd[-1]/2).argmin()

        # Calculate the plume bearing
        volc_loc = [vlat, vlon]
        peak_loc = [lat[peak_idx], lon[peak_idx]]
        plume_dist, plume_bearing = haversine(volc_loc, peak_loc)

        logger.info(f'Distance from vent: {plume_dist:.02f} m')
        logger.info(f'Azimuth from vent: {np.degrees(plume_bearing):.02f} deg')

        # Calculate the distance and bearing of each measurement vector
        vect = [haversine([lat[i-1], lon[i-1]], [lat[i], lon[i]])
                for i in range(1, len(lat))]

        # Unpack the distance and bearing from the vectors
        trav_dist, trav_bearing = np.asarray(vect).T

        # Correct the distance for the angle between the travel and plume
        # vectors
        corr_factor = np.sin(plume_bearing-trav_bearing)
        corr_dist = np.multiply(trav_dist, corr_factor)

        # Convert so2 amounts from molecules.cm-2 to molecules.m-2
        so2_molec_per_m2 = so2_scd * 1.0e4

        # Multiply by the distance moved and sum
        so2_molec_per_m = np.nansum(np.multiply(so2_molec_per_m2[1:],
                                                corr_dist))

        # Multiply by the wind speed
        so2_molec_per_s = so2_molec_per_m * wind_speed

        # Convert to moles
        so2_moles_per_s = so2_molec_per_s / 6.022e23

        # Convert to kg/s. Molar mass of SO2 is 64.066 g/mole
        so2_kg_per_s = so2_moles_per_s * 0.064066
        logger.info(f'Total SO2 mass: {abs(so2_kg_per_s/wind_speed):.02f} kg')

        # Convert to t/day if required
        if flux_units == 't/day':
            flux = abs(so2_kg_per_s * 60*60*24 / 1000.0)
        else:
            flux = abs(so2_kg_per_s)

        logger.info(f'SO2 Flux: {flux:.02f} {flux_units}')

        # Calculate the Flux Error
        so2_err[abs(so2_err) == np.inf] = np.nan
        tot_so2_err = np.nansum(so2_err)
        frac_so2_err = tot_so2_err / np.nansum(so2_scd)

        # Combine with the wind speed error
        frac_err = ((frac_so2_err)**2 + (wind_error)**2)**0.5
        flux_err = flux * frac_err

        # Add the flux result to the table
        self.fluxTable.setRowCount(self.fluxTable.rowCount()+1)
        self.fluxTable.setItem(self.trav_no, 0,
                               QTableWidgetItem(f'{flux:.02f}'))
        self.fluxTable.setItem(self.trav_no, 1,
                               QTableWidgetItem(f'{flux_err:.02f}'))

        # Plot the traverse graphs
        self.map_elements['trav_line'].setData(lon, lat)
        self.map_elements['plume_start'].setData([lon[0]], [lat[0]])
        self.map_elements['plume_cent'].setData([peak_loc[1]], [peak_loc[0]])
        self.map_elements['plume_stop'].setData([lon[-1]], [lat[-1]])
        self.map_elements['volcano'].setData([vlon], [vlat])

        # Collate the results
        output_data = np.column_stack([so2_time[1:], so2_scd[1:], so2_err[1:],
                                       lat[1:], lon[1:], trav_dist,
                                       trav_bearing, corr_factor])
        self.flux_data[f'trav{self.trav_no}'] = {'flux': flux,
                                                 'flux_err': flux_err,
                                                 'flux_units': flux_units,
                                                 'wind_speed': wind_speed,
                                                 'vlat': vlat,
                                                 'vlon': vlon,
                                                 'peak_loc': peak_loc,
                                                 'output_data': output_data}

        # Increment the counter
        self.trav_no += 1

        self.save_flag = True

    def del_trav(self):
        """Delete the last traverse."""
        if self.trav_no > 0:
            logger.info('Removing last traverse')
            self.trav_no -= 1
            self.flux_data.pop(f'trav{self.trav_no}')
            self.fluxTable.setRowCount(self.fluxTable.rowCount()-1)
        if self.trav_no == 0:
            self.save_flag = False

    def save_fluxes(self):
        """Output the flux results."""
        # Make sure the output directory exists, and create if not
        out_path = self.widgets.get('out_path')
        if not os.path.isdir(out_path):
            os.makedirs(out_path)

        for i, [key, data] in enumerate(self.flux_data.items()):

            # Write the detailed output
            with open(out_path + 'flux_results.csv', 'w') as w:
                w.write(f'Traverse Number,{i+1}\n'
                        + f'Flux ({data["flux_units"]}),{data["flux"]}\n'
                        + f'Error ({data["flux_units"]}),{data["flux_err"]}\n'
                        + f'Wind Speed (m/s), {data["wind_speed"]}\n'
                        + f'Volcano Lat/Lon,{data["vlat"]},{data["vlon"]}\n'
                        + f'Plume Center,{data["peak_loc"][0]},'
                        + f'{data["peak_loc"][1]}\n'
                        + 'Time (Decimal hours),SO2 SCD,SO2 Err,Lat,Lon,'
                        + 'Distance,Bearing,Correction Factor\n')

                for i, line in enumerate(data["output_data"]):
                    w.write(str(line[0]))
                    for j in line[1:]:
                        w.write(f',{j}')
                    w.write('\n')
                w.write('\n')

        # Write summary output
        flux_units = data["flux_units"]
        fluxes = [data['flux'] for data in self.flux_data.values()]
        errors = [data['flux_err'] for data in self.flux_data.values()]
        w_mean_flux = np.average(fluxes, weights=np.power(errors, -2))
        w_mean_error = (1 / np.sum(np.power(errors, -2)))**0.5
        with open(out_path + 'flux_summary.txt', 'w') as w:
            w.write(f'Flux Summary\nFlux ({flux_units}), Error\n')
            for i in range(len(fluxes)):
                w.write(f'{fluxes[i]}, {errors[i]}\n')
            w.write(f'Av Flux = {np.average(fluxes):.02f} {flux_units}\n')
            w.write(f'Stdev = {np.std(fluxes):.02f} {flux_units}\n')
            w.write(f'Weighted Mean Flux = {w_mean_flux:.02f}'
                    + f' (+/- {w_mean_error:.02f}) {flux_units}')

        logger.info('Fluxes saved')
        self.save_flag = False

    def closeEvent(self, event):
        """Handle GUI closure."""
        if self.save_flag:
            options = QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
            reply = QMessageBox.question(self, 'Message',
                                         "Would you like to save the fluxes?",
                                         options, QMessageBox.Cancel)

            if reply == QMessageBox.Yes:
                self.save_fluxes()
                self.save_config()
                event.accept()
            elif reply == QMessageBox.No:
                event.accept()
                self.save_config()
            else:
                event.ignore()
        else:
            self.save_config()

    def load_config(self):
        """Load previous session settings."""
        # Open the config file
        try:
            with open('bin/.calc_flux_config.yml', 'r') as ymlfile:
                config = yaml.load(ymlfile, Loader=yaml.FullLoader)

            for key, item in config.items():
                try:
                    self.widgets.set(key, item)
                except Exception:
                    logger.warning(f'Failed to load {key} from config file')
        except FileNotFoundError:
            logger.warning('Unable to load config file!')

    def save_config(self):
        """Save session settings."""
        config = {}
        for key in self.widgets:
            config[key] = self.widgets.get(key)

        # Write the config
        with open('bin/.calc_flux_config.yml', 'w') as outfile:
            yaml.dump(config, outfile)


# =============================================================================
# Measure ILS
# =============================================================================

class ILSWindow(QMainWindow):
    """Opens ILS analysis window."""

    def __init__(self, parent=None):
        """Initialise the ILS window."""
        super(ILSWindow, self).__init__(parent)

        # Set the window properties
        self.setWindowTitle('Measure ILS')
        self.statusBar().showMessage('Ready')
        self.setGeometry(40, 40, 1000, 600)
        self.setWindowIcon(QIcon('bin/icons/ils.png'))

        # Set the window layout
        self.generalLayout = QGridLayout()
        self._centralWidget = QScrollArea()
        self.widget = QWidget()
        self.setCentralWidget(self._centralWidget)
        self.widget.setLayout(self.generalLayout)

        # Scroll Area Properties
        self._centralWidget.setWidgetResizable(True)
        self._centralWidget.setWidget(self.widget)

        self._createApp()

    def _createApp(self):

        layout = QGridLayout(self._centralWidget)

        # Add option for spectra type
        layout.addWidget(QLabel('Format:'), 0, 0)
        self.spec_type = QComboBox()
        self.spec_type.setToolTip('Choose spectrum format')
        self.spec_type.addItems(['iFit',
                                 'Master.Scope',
                                 'Spectrasuite',
                                 'mobileDOAS',
                                 'Basic'])
        self.spec_type.setFixedSize(100, 20)
        layout.addWidget(self.spec_type, 0, 1)

        # Add an input for the spectra selection
        layout.addWidget(QLabel('Spectra:'), 1, 0)
        self.spec_fnames = QTextEdit()
        self.spec_fnames.setFixedSize(300, 150)
        layout.addWidget(self.spec_fnames, 1, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.spec_fnames, 'multi'))
        layout.addWidget(btn, 1, 4)

        # Add an input for the dark spectra selection
        layout.addWidget(QLabel('Dark\nSpectra:'), 2, 0)
        self.dark_fnames = QTextEdit()
        self.dark_fnames.setFixedSize(300, 150)
        layout.addWidget(self.dark_fnames, 2, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.dark_fnames, 'multi'))
        layout.addWidget(btn, 2, 4)

        # Add an input for the save selection
        layout.addWidget(QLabel('Save:'), 3, 0)
        self.save_path = QLineEdit()
        layout.addWidget(self.save_path, 3, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.save_path,
                                    'save'))
        layout.addWidget(btn, 3, 4)

        btn = QPushButton('Import')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(self.import_spectra)
        layout.addWidget(btn, 4, 1)

        btn = QPushButton('Fit')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(self.measure_ils)
        layout.addWidget(btn, 4, 2)

        btn = QPushButton('Save')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(self.save_fit)
        layout.addWidget(btn, 4, 3)

        graphwin = pg.GraphicsLayoutWidget(show=True)
        pg.setConfigOptions(antialias=True)
        # pg.setConfigOptions(useOpenGL=True)

        w = QWidget()
        glayout = QGridLayout(w)

        # Make the graphs
        self.ax0 = graphwin.addPlot(row=0, col=0)
        self.ax1 = graphwin.addPlot(row=1, col=0)
        self.plot_axes = [self.ax0, self.ax1]

        for ax in self.plot_axes:
            ax.setDownsampling(mode='peak')
            ax.setClipToView(True)
            ax.showGrid(x=True, y=True)

        # Add axis labels
        self.ax0.setLabel('left', 'Intensity (counts)')
        self.ax1.setLabel('left', 'Intensity (arb)')
        self.ax1.setLabel('bottom', 'Wavelength (nm)')

        # Add the graphs to the layout
        glayout.addWidget(graphwin, 0, 0, 1, 7)
        layout.addWidget(w, 0, 5, 5, 1)

    def import_spectra(self):
        """Read in ILS spectra and display."""
        # Read in the dark spectra
        files = self.dark_fnames.toPlainText()
        if files == '':
            dark = 0
        else:
            x, dark = average_spectra(files.split('\n'),
                                      self.spec_type.currentText())
            # for i, fname in enumerate(files.split('\n')):
            #     x, y = np.loadtxt(fname, unpack=True)
            #     if i == 0:
            #         spec = y
            #     else:
            #         spec += y
            #
            #     dark = np.divide(spec, i+1)

        # Read in the measurement spectra
        for i, fname in enumerate(self.spec_fnames.toPlainText().split('\n')):
            x, y = np.loadtxt(fname, unpack=True)
            if i == 0:
                spec = y
            else:
                spec += y

        self.x = x
        self.spec = np.divide(spec, i+1) - dark

        # Add to the graph
        self.ax0.clear()
        self.l0 = self.ax0.plot(self.x, self.spec,
                                pen=pg.mkPen(color='#1f77b4', width=1.0))

        # Add the region selector
        self.lr = pg.LinearRegionItem([x[0], x[-1]])
        self.lr.setZValue(-10)
        self.ax0.addItem(self.lr)

    def measure_ils(self):
        """Measure the ILS on the selected line."""
        # Get the highlighted region
        i0, i1 = self.lr.getRegion()
        idx = np.where(np.logical_and(self.x >= i0, self.x <= i1))
        grid = self.x[idx]
        line = self.spec[idx]

        # Find the line
        center = grid[line.argmax()]

        # Center grid on zero
        grid = grid - center

        # Remove the offset and normalise
        line = line - min(line)
        line = np.divide(line, max(line))

        # Fit super gaussian
        p0 = [0.34, 2, 0.0, 0.0, 0.0, 1.0, 0.0]
        self.popt, pcov = curve_fit(super_gaussian, grid, line, p0=p0)
        ngrid = np.arange(grid[0], grid[-1]+0.01, 0.01)
        fit = super_gaussian(ngrid, *self.popt)

        # Plot
        self.ax1.clear()
        self.ax1.plot(grid, line, pen=None, symbol='+', symbolPen=None,
                      symbolSize=10, symbolBrush='#1f77b4')
        self.ax1.plot(ngrid, fit, pen=pg.mkPen(color='#ff7f0e', width=1.0))

    def save_fit(self):
        """Save the ILS fit parameters."""
        w, k, a_w, a_k, shift, amp, offset = self.popt

        fwem = 2*w

        np.savetxt(self.save_path.text(), [fwem, k, a_w, a_k],
                   header='ILS parameters (FWEM, k, a_w, a_k)')


# =============================================================================
# Measure Flat Spectrum
# =============================================================================

class FLATWindow(QMainWindow):
    """Open a window for flat-field analysis."""

    def __init__(self, parent=None):
        """Initialise the window."""
        super(FLATWindow, self).__init__(parent)

        # Set the window properties
        self.setWindowTitle('Measure Flat Field')
        self.statusBar().showMessage('Ready')
        self.setGeometry(40, 40, 1000, 500)
        self.setWindowIcon(QIcon('bin/icons/flat.png'))

        # Set the window layout
        self.generalLayout = QGridLayout()
        self._centralWidget = QScrollArea()
        self.widget = QWidget()
        self.setCentralWidget(self._centralWidget)
        self.widget.setLayout(self.generalLayout)

        # Scroll Area Properties
        self._centralWidget.setWidgetResizable(True)
        self._centralWidget.setWidget(self.widget)

        self._createApp()

    def _createApp(self):

        layout = QGridLayout(self._centralWidget)

        # Add an input for the spectra selection
        layout.addWidget(QLabel('Spectra:'), 0, 0)
        self.spec_fnames = QTextEdit()
        self.spec_fnames.setFixedSize(300, 150)
        layout.addWidget(self.spec_fnames, 0, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.spec_fnames, 'multi'))
        layout.addWidget(btn, 0, 4)

        # Add an input for the dark spectra selection
        layout.addWidget(QLabel('Dark\nSpectra:'), 1, 0)
        self.dark_fnames = QTextEdit()
        self.dark_fnames.setFixedSize(300, 150)
        layout.addWidget(self.dark_fnames, 1, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.dark_fnames, 'multi'))
        layout.addWidget(btn, 1, 4)

        # Add an input for the save selection
        layout.addWidget(QLabel('Save:'), 2, 0)
        self.save_path = QLineEdit()
        layout.addWidget(self.save_path, 2, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.save_path,
                                    'save'))
        layout.addWidget(btn, 2, 4)

        btn = QPushButton('Import')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(self.import_spectra)
        layout.addWidget(btn, 3, 1)

        btn = QPushButton('Fit')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(self.measure_flat)
        layout.addWidget(btn, 3, 2)

        btn = QPushButton('Save')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(self.save_fit)
        layout.addWidget(btn, 3, 3)

        graphwin = pg.GraphicsLayoutWidget(show=True)
        pg.setConfigOptions(antialias=True)
        # pg.setConfigOptions(useOpenGL=True)

        w = QWidget()
        glayout = QGridLayout(w)

        # Make the graphs
        self.ax0 = graphwin.addPlot(row=0, col=0)
        self.ax1 = graphwin.addPlot(row=1, col=0)
        self.ax2 = graphwin.addPlot(row=2, col=0)
        self.plot_axes = [self.ax0, self.ax1, self.ax2]

        for ax in self.plot_axes:
            ax.setDownsampling(mode='peak')
            ax.setClipToView(True)
            ax.showGrid(x=True, y=True)

        # Add axis labels
        self.ax0.setLabel('left', 'Intensity (counts)')
        self.ax1.setLabel('left', 'Intensity (counts)')
        self.ax2.setLabel('left', 'Flat Response')
        self.ax2.setLabel('bottom', 'Wavelength (nm)')

        # Add the graphs to the layout
        glayout.addWidget(graphwin, 0, 0, 1, 7)
        layout.addWidget(w, 0, 5, 5, 1)

    def import_spectra(self):
        """Read in ILS spectra and display."""
        # Read in the dark spectra
        files = self.dark_fnames.toPlainText()
        if files == '':
            dark = 0
        else:
            for i, fname in enumerate(files.split('\n')):
                x, y = np.loadtxt(fname, unpack=True)
                if i == 0:
                    spec = y
                else:
                    spec += y

                dark = np.divide(spec, i+1)

        # Read in the measurement spectra
        for i, fname in enumerate(self.spec_fnames.toPlainText().split('\n')):
            x, y = np.loadtxt(fname, unpack=True)
            if i == 0:
                spec = y
            else:
                spec += y

        self.x = x
        self.y = np.divide(spec, i+1) - dark

        # Add to the graph
        self.ax0.clear()
        self.ax1.clear()
        self.ax2.clear()
        self.l0 = self.ax0.plot(self.x, self.y,
                                pen=pg.mkPen(color='#1f77b4', width=1.0))

        # Add the region selector
        self.lr = pg.LinearRegionItem([x[0], x[-1]])
        self.lr.setZValue(-10)
        self.ax0.addItem(self.lr)

    def measure_flat(self):
        """Measure the flat spectrum across the selected region."""
        width = 5

        # Get the highlighted region
        i0, i1 = self.lr.getRegion()
        idx = np.where(np.logical_and(self.x >= i0, self.x <= i1))
        self.grid = self.x[idx]
        self.spec = self.y[idx]

        # Create the boxcar window
        window = np.ones(width) / width

        # Pad the array with values to avoid edge effects
        pre_array = np.ones(width-1) * self.spec[0]
        post_array = np.ones(width-1) * self.spec[-1]

        # Add padding to the origional array
        spec = np.append(pre_array, self.spec)
        spec = np.append(spec, post_array)

        # Convolve with boxcar to smooth
        smooth_spec = np.convolve(spec, window, 'same')

        # Cut array to origional size
        smooth_spec = smooth_spec[width-1:-(width-1)]

        self.ax1.plot(self.grid, self.spec, pen=pg.mkPen(color='#1f77b4'))
        self.ax1.plot(self.grid, smooth_spec, pen=pg.mkPen(color='#ff7f0e'))

        self.flat = np.divide(self.spec, smooth_spec)

        self.ax2.plot(self.grid, self.flat, pen=pg.mkPen(color='#1f77b4'))

    def save_fit(self):
        """Save the flat response."""
        data = np.column_stack([self.grid, self.flat])
        header = 'Flat spectrum\nWavelength (nm),       Flat Response'
        np.savetxt(self.save_path.text(), data, header=header)


# =============================================================================
# Measure Light Dilution
# =============================================================================

class LDFWindow(QMainWindow):
    """Open a window for light dilution analysis."""

    def __init__(self, widgetData, parent=None):
        """Initialise the window."""
        super(LDFWindow, self).__init__(parent)

        # Set the window properties
        self.setWindowTitle('Light Dilution Analysis')
        self.statusBar().showMessage('Ready')
        self.setGeometry(40, 40, 1200, 700)
        self.setWindowIcon(QIcon('bin/icons/ldf.png'))

        # Generate the threadpool for launching background processes
        self.threadpool = QThreadPool()

        # Set the window layout
        self.generalLayout = QGridLayout()
        self._centralWidget = QScrollArea()
        self.widget = QWidget()
        self.setCentralWidget(self._centralWidget)
        self.widget.setLayout(self.generalLayout)

        # Scroll Area Properties
        self._centralWidget.setWidgetResizable(True)
        self._centralWidget.setWidget(self.widget)

        # Save the main GUI widget data
        self.widgetData = widgetData

        self._createApp()

    def _createApp(self):
        """Create the main app widgets."""
        # Create a frame to hold program controls
        self.controlFrame = QFrame(self)
        self.controlFrame.setFrameShape(QFrame.StyledPanel)

        # Create a frame to hold program outputs
        self.outputFrame = QFrame(self)
        self.outputFrame.setFrameShape(QFrame.StyledPanel)

        # Create a frame to hold graphs
        self.graphFrame = QFrame(self)
        self.graphFrame.setFrameShape(QFrame.StyledPanel)

        # Add splitters to allow for adjustment
        splitter1 = QSplitter(Qt.Vertical)
        splitter1.addWidget(self.controlFrame)
        splitter1.addWidget(self.outputFrame)

        splitter2 = QSplitter(Qt.Horizontal)
        splitter2.addWidget(splitter1)
        splitter2.addWidget(self.graphFrame)

        # Pack the Frames and splitters
        self.generalLayout.addWidget(splitter2)

        # Create the individual widgets
        self._createControls()
        self._createOutput()
        self._createGraphs()

# =============================================================================
#   Program controls
# =============================================================================

    def _createControls(self):
        """Create main analysis controls."""
        # Setup tab layout
        tablayout = QGridLayout(self.controlFrame)

        # Generate tabs for the gaphs and settings
        tab1 = QWidget()
        tab2 = QWidget()

        # Form the tab widget
        tabwidget = QTabWidget()
        tabwidget.addTab(tab1, 'Generate Curves')
        tabwidget.addTab(tab2, 'Load Data')
        tablayout.addWidget(tabwidget, 0, 0)

        # Generate Curves =====================================================

        # Setup the main layout
        layout = QGridLayout(tab1)
        layout.setAlignment(Qt.AlignTop)

        # Create an option menu for the spectra format
        layout.addWidget(QLabel('Format:'), 0, 0)
        self.spec_type = QComboBox()
        self.spec_type.addItems(['iFit',
                                 'Master.Scope',
                                 'Spectrasuite',
                                 'mobileDOAS',
                                 'Basic'])
        self.spec_type.setFixedSize(100, 20)
        layout.addWidget(self.spec_type, 0, 1)
        index = self.spec_type.findText(self.widgetData['spec_type'],
                                        Qt.MatchFixedString)
        if index >= 0:
            self.spec_type.setCurrentIndex(index)

        # Add an input for the spectra selection
        layout.addWidget(QLabel('Spectra:'), 1, 0)
        self.spec_fnames = QTextEdit()
        self.spec_fnames.setFixedHeight(75)
        layout.addWidget(self.spec_fnames, 1, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.spec_fnames, 'multi'))
        layout.addWidget(btn, 1, 4)

        # Add an input for the dark selection
        layout.addWidget(QLabel('Dark Spectra:'), 2, 0)
        self.dark_fnames = QTextEdit()
        self.dark_fnames.setFixedHeight(75)
        layout.addWidget(self.dark_fnames, 2, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.dark_fnames, 'multi'))
        layout.addWidget(btn, 2, 4)

        # Add spinboxs for the fit windows
        layout.addWidget(QLabel('Fit Window 1:\n    (nm)'), 3, 0, 2, 1)
        self.wb1_lo = DSpinBox(306, [0, 10000], 0.1)
        self.wb1_lo.setFixedSize(70, 20)
        layout.addWidget(self.wb1_lo, 3, 1)
        self.wb1_hi = DSpinBox(316, [0, 10000], 0.1)
        self.wb1_hi.setFixedSize(70, 20)
        layout.addWidget(self.wb1_hi, 4, 1)

        layout.addWidget(QLabel('Fit Window 2:\n    (nm)'), 3, 2, 2, 1)
        self.wb2_lo = DSpinBox(312, [0, 10000], 0.1)
        self.wb2_lo.setFixedSize(70, 20)
        layout.addWidget(self.wb2_lo, 3, 3)
        self.wb2_hi = DSpinBox(322, [0, 10000], 0.1)
        self.wb2_hi.setFixedSize(70, 20)
        layout.addWidget(self.wb2_hi, 4, 3)

        # Add entries for the SO2 grid
        layout.addWidget(QLabel('SO<sub>2</sub> Grid Low:'), 5, 0)
        self.so2_grid_lo = QLineEdit()
        self.so2_grid_lo.setText('0.0')
        self.so2_grid_lo.setFixedSize(100, 20)
        layout.addWidget(self.so2_grid_lo, 5, 1)
        layout.addWidget(QLabel('SO<sub>2</sub> Grid High:'), 6, 0)
        self.so2_grid_hi = QLineEdit()
        self.so2_grid_hi.setText('1.0e19')
        self.so2_grid_hi.setFixedSize(100, 20)
        layout.addWidget(self.so2_grid_hi, 6, 1)
        layout.addWidget(QLabel('SO<sub>2</sub> Grid Step:'), 7, 0)
        self.so2_grid_step = QLineEdit()
        self.so2_grid_step.setText('5.0e17')
        self.so2_grid_step.setFixedSize(100, 20)
        layout.addWidget(self.so2_grid_step, 7, 1)

        # Add entries for the LDF grid
        layout.addWidget(QLabel('LDF Grid Low:'), 5, 2)
        self.ldf_grid_lo = QLineEdit()
        self.ldf_grid_lo.setText('0.0')
        self.ldf_grid_lo.setFixedSize(100, 20)
        layout.addWidget(self.ldf_grid_lo, 5, 3)
        layout.addWidget(QLabel('LDF Grid High:'), 6, 2)
        self.ldf_grid_hi = QLineEdit()
        self.ldf_grid_hi.setText('0.9')
        self.ldf_grid_hi.setFixedSize(100, 20)
        layout.addWidget(self.ldf_grid_hi, 6, 3)
        layout.addWidget(QLabel('LDF Grid Step:'), 7, 2)
        self.ldf_grid_step = QLineEdit()
        self.ldf_grid_step.setText('0.1')
        self.ldf_grid_step.setFixedSize(100, 20)
        layout.addWidget(self.ldf_grid_step, 7, 3)

        # Add an input for the save selection
        layout.addWidget(QLabel('Output File:'), 8, 0)
        self.save_path = QLineEdit()
        self.save_path.setFixedHeight(25)
        layout.addWidget(self.save_path, 8, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.save_path,
                                    'save', "Text (*.txt)"))
        layout.addWidget(btn, 8, 4)

        # Add button to begin analysis
        self.start_btn = QPushButton('Begin!')
        self.start_btn.clicked.connect(partial(self.calc_ld_curves))
        self.start_btn.setFixedSize(90, 25)
        layout.addWidget(self.start_btn, 9, 1)

        # Add button to pause analysis
        self.save_btn = QPushButton('Save')
        self.save_btn.clicked.connect(partial(self.save_ld_curves))
        self.save_btn.setFixedSize(90, 25)
        self.save_btn.setEnabled(False)
        layout.addWidget(self.save_btn, 9, 2)

        # Load Curves =========================================================

        # Setup the main layout
        layout = QGridLayout(tab2)
        layout.setAlignment(Qt.AlignTop)

        # Add an input for the loading selection
        layout.addWidget(QLabel('Light Dilution\nCurve File:'), 0, 0)
        self.load_path = QLineEdit()
        self.load_path.setFixedHeight(25)
        layout.addWidget(self.load_path, 0, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.load_path,
                                    'single', "Text (*.txt)"))
        layout.addWidget(btn, 0, 4)
        btn = QPushButton('Load')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(self.load_ld_curves)
        layout.addWidget(btn, 0, 5)

        layout.addWidget(QHLine(), 1, 0, 1, 5)

        # Add an input for the iFit output files
        layout.addWidget(QLabel('iFit Output:\n(W1):'), 2, 0)
        self.ifit_1_path = QLineEdit()
        self.ifit_1_path.setFixedHeight(25)
        layout.addWidget(self.ifit_1_path, 2, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.ifit_1_path,
                                    'single', "Comma separated (*.csv)"))
        layout.addWidget(btn, 2, 4)

        layout.addWidget(QLabel('iFit Output:\n(W2):'), 3, 0)
        self.ifit_2_path = QLineEdit()
        self.ifit_2_path.setFixedHeight(25)
        layout.addWidget(self.ifit_2_path, 3, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.ifit_2_path,
                                    'single', "Comma separated (*.csv)"))
        layout.addWidget(btn, 3, 4)

        btn = QPushButton('Load')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(self.load_ifit_data)
        layout.addWidget(btn, 2, 5, 2, 1)

# =============================================================================
#   Program outputs
# =============================================================================

    def _createOutput(self):
        """Create program output widgets."""
        # Generate the layout
        layout = QGridLayout(self.outputFrame)

        # Create a textbox to display the program messages
        self.logBox = QTextEditLogger(self)
        self.logBox.setFormatter(logging.Formatter('%(message)s'))
        logging.getLogger().addHandler(self.logBox)
        logging.getLogger().setLevel(logging.INFO)
        layout.addWidget(self.logBox.widget, 0, 0)

# =============================================================================
#   Graphs
# =============================================================================

    def _createGraphs(self):
        """Generate the graphs."""
        # Generate the graph layout and window
        layout = QGridLayout(self.graphFrame)
        pg.setConfigOptions(antialias=True)
        graphwin = pg.GraphicsLayoutWidget(show=True)

        # Make the graphs
        self.ax = graphwin.addPlot()

        self.ax.setDownsampling(mode='peak')
        self.ax.setClipToView(True)
        self.ax.showGrid(x=True, y=True)

        # Add axis labels
        self.ax.setLabel('left', 'SO<sub>2</sub> W2 (molec cm<sup>-2</sup>)')
        self.ax.setLabel('bottom', 'SO<sub>2</sub> W1 (molec cm<sup>-2</sup>)')

        # Create an order parameter to scale the graphs
        self.order = None

        # Add the graphs to the layout
        layout.addWidget(graphwin, 0, 0, 0, 0)

# =============================================================================
#   Program functions
# =============================================================================

    def calc_ld_curves(self):
        """Calculate Light Dilution curve data."""
        # Disable buttons
        self.start_btn.setEnabled(False)
        self.save_btn.setEnabled(False)

        # Pull the waveband data from the LD gui and pad for the initial fit
        self.widgetData['fit_lo'] = self.wb1_lo.value() - 1
        self.widgetData['fit_hi'] = self.wb2_hi.value() + 1

        # Grab the spectra type from the GUI
        self.widgetData['spec_type'] = self.spec_type.currentText()

        # Pull the waveband, SO2 and LDF grid info from the GUI
        ld_kwargs = {'wb1': [self.wb1_lo.value(), self.wb1_hi.value()],
                     'wb2': [self.wb2_lo.value(), self.wb2_hi.value()],
                     'so2_lims': [float(self.so2_grid_lo.text()),
                                  float(self.so2_grid_hi.text())],
                     'so2_step': float(self.so2_grid_step.text()),
                     'ldf_lims': [float(self.ldf_grid_lo.text()),
                                  float(self.ldf_grid_hi.text())],
                     'ldf_step': float(self.ldf_grid_step.text())}

        # Load the spectra to fit
        spec_fnames = self.spec_fnames.toPlainText().split('\n')
        dark_fnames = self.dark_fnames.toPlainText().split('\n')

        # Initialise the analysis worker
        self.ldThread = QThread()
        self.ldWorker = LDWorker(spec_fnames, dark_fnames, self.widgetData,
                                 ld_kwargs)
        self.ldWorker.moveToThread(self.ldThread)
        self.ldThread.started.connect(self.ldWorker.run)
        self.ldWorker.finished.connect(self.analysis_complete)
        self.ldWorker.data.connect(self.handle_data)
        self.ldWorker.finished.connect(self.ldThread.quit)
        self.ldWorker.finished.connect(self.ldWorker.deleteLater)
        self.ldThread.finished.connect(self.ldThread.deleteLater)
        self.ldThread.start()

    def handle_data(self, ld_results):
        """Handle the LD results."""
        self.ld_results = ld_results
        self.plot_ld_curves(ld_results)

    def analysis_complete(self):
        """Run when LD analysis is complete."""
        self.start_btn.setEnabled(True)
        self.save_btn.setEnabled(True)

    def save_ld_curves(self):
        """Save Light Dilution curve data."""
        if self.save_path.text() == '':
            filter = "Text (*.txt)"
            fname, _ = QFileDialog.getSaveFileName(self, 'Save As', '', filter)
        else:
            fname = self.save_path.text()

        header = 'LDF Curves generated by iFit. All SCDs in [molec/cm^2]\n' \
                 + f'W1: {self.wb1_lo.value()} - {self.wb1_hi.value()} nm\n' \
                 + f'W2: {self.wb2_lo.value()} - {self.wb2_hi.value()} nm\n' \
                 + 'LDF\tModel_SO2_SCD\tW1_SO2_SCD\tW1_SO2_Err\t' \
                 + 'W2_SO2_SCD\tW2_SO2_Err'
        np.savetxt(fname, self.ld_results, header=header)

        logger.info('Light dilution curve data saved!')

    def load_ld_curves(self):
        """Load the light dilution curve data."""
        if self.load_path.text() == '':
            filter = "Text (*.txt)"
            fname, _ = QFileDialog.getOpenFileName(self, 'Load', '', filter)
        else:
            fname = self.load_path.text()

        self.ld_results = np.loadtxt(fname)
        self.plot_ld_curves(self.ld_results)
        logger.info('Light dilution curve data loaded!')

    def load_ifit_data(self):
        """Load the light dilution curve data."""
        # Read in the iFit analysis results
        df1 = pd.read_csv(self.ifit_1_path.text(), comment='#')
        df2 = pd.read_csv(self.ifit_2_path.text(), comment='#')

        if self.order is None:
            max_val = np.max([np.nanmax(df1['SO2']), np.nanmax(df2['SO2'])])

            if ~np.isnan(max_val) and max_val > 1e6:
                self.order = int(np.ceil(np.log10(max_val))) - 1
                self.ax.setLabel(
                    'left',
                    'SO<sub>2</sub> W2 (molec cm<sup>-2</sup>) '
                    + f'(1e{self.order})'
                )
                self.ax.setLabel(
                    'bottom',
                    'SO<sub>2</sub> W1 (molec cm<sup>-2</sup>) '
                    + f'(1e{self.order})'
                )

            else:
                self.order = 0
                self.ax.setLabel(
                    'left', 'SO<sub>2</sub> W2 (molec cm<sup>-2</sup>)')
                self.ax.setLabel(
                    'bottom', 'SO<sub>2</sub> W1 (molec cm<sup>-2</sup>)')

        plotx = df1['SO2'] / 10**self.order
        ploty = df2['SO2'] / 10**self.order

        # Plot on the graph
        p0 = pg.mkPen(color='#1f77b4', width=0.5)
        line = self.ax.plot(plotx, ploty, pen=None, symbolPen=p0,
                            symbol='o', brush=None)
        line.setAlpha(0.75, False)

        logger.info('iFit data loaded!')

    def plot_ld_curves(self, ld_results):
        """Plot the light dilution curves."""
        # Clear the plot
        self.ax.clear()
        self.order = None

        # Pull out the unique LDF values
        ldf_grid = np.unique(ld_results[:, 0])

        # Check for large number in the time series. This is due to
        # a bug in pyqtgraph not displaying large numbers
        if self.order is None:
            max_val = np.max([np.nanmax(np.abs(ld_results[:, 2])),
                              np.nanmax(np.abs(ld_results[:, 4]))])

            # Calculate the required order and update the axes labels
            if ~np.isnan(max_val) and max_val > 1e6:
                self.order = int(np.ceil(np.log10(max_val))) - 1
                self.ax.setLabel(
                    'left',
                    'SO<sub>2</sub> W2 (molec cm<sup>-2</sup>) '
                    + f'(1e{self.order})'
                )
                self.ax.setLabel(
                    'bottom',
                    'SO<sub>2</sub> W1 (molec cm<sup>-2</sup>) '
                    + f'(1e{self.order})'
                )
            else:
                self.order = 0
                self.ax.setLabel(
                    'left', 'SO<sub>2</sub> W2 (molec cm<sup>-2</sup>)')
                self.ax.setLabel(
                    'bottom', 'SO<sub>2</sub> W1 (molec cm<sup>-2</sup>)')

        # Get the curve data for each LDF value
        legend = self.ax.addLegend()
        for i, ldf in enumerate(ldf_grid):
            row_idx = np.where(ld_results[:, 0] == ldf)[0]
            so2_scd_1 = ld_results[row_idx, 2]
            # so2_err_1 = ld_results[row_idx, 3]
            so2_scd_2 = ld_results[row_idx, 4]
            # so2_err_2 = ld_results[row_idx, 5]

            # Adjust for large numbers if neccissary
            if self.order is not None:
                so2_scd_1 = so2_scd_1 / 10**self.order
                so2_scd_2 = so2_scd_2 / 10**self.order

            # Plot!
            line = self.ax.plot(so2_scd_1, so2_scd_2,
                                pen=pg.mkPen(color=COLORS[i % 10], width=2.0))
            legend.addItem(line, name=f'LDF={ldf:.02f}')


# =============================================================================
# Useful functions
# =============================================================================

def browse(gui, widget, mode='single', filter=False):
    """Open file dialouge."""
    if not filter:
        filter = None
    else:
        filter = filter + ';;All Files (*)'

    if mode == 'single':
        fname, _ = QFileDialog.getOpenFileName(gui, 'Select File', '',
                                               filter)
        if fname != '':
            widget.setText(fname)

    elif mode == 'multi':
        fnames, _ = QFileDialog.getOpenFileNames(gui, 'Select Files', '',
                                                 filter)
        if fnames != []:
            widget.setText('\n'.join(fnames))

    elif mode == 'save':
        fname, _ = QFileDialog.getSaveFileName(gui, 'Save As', '', filter)
        if fname != '':
            widget.setText(fname)

    elif mode == 'folder':
        fname = QFileDialog.getExistingDirectory(gui, 'Select Foler')
        if fname != '':
            widget.setText(fname + '/')


# Create a worker to handle QThreads for light dilution analysis
class LDWorker(QObject):
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

    # Define signals
    finished = pyqtSignal()
    data = pyqtSignal(np.ndarray)
    error = pyqtSignal(tuple)

    def __init__(self, spec_fnames, dark_fnames, widgetData, ld_kwargs):
        """Initialise."""
        super(QObject, self).__init__()

        # Create stopped and paused flags
        self.is_paused = False
        self.is_killed = False

        # Create a holder for the spectrum filepath
        self.spec_fnames = spec_fnames
        self.dark_fnames = dark_fnames
        self.widgetData = widgetData
        self.ld_kwargs = ld_kwargs

    def run(self):
        """Launch worker function."""
        try:
            self._run()
        except Exception:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.error.emit((exctype, value, traceback.format_exc()))
        self.finished.emit()

    def _run(self):
        """Launch LD analysis from the GUI."""
        # Read in the spectra
        spectrum = average_spectra(self.spec_fnames,
                                   self.widgetData['spec_type'],
                                   self.widgetData['wl_calib'])

        if self.widgetData['dark_flag']:
            x, dark = average_spectra(self.dark_fnames,
                                      self.widgetData['spec_type'],
                                      self.widgetData['wl_calib'])
        else:
            dark = 0

        # Generate the analyser
        logger.info('Generating the iFit analyser...')
        self.analyser = self.generate_analyser()
        self.analyser.dark_spec = dark

        # Generate the light dilution curves
        logger.info('Beginning light dilution calculations')
        ld_results = generate_ld_curves(self.analyser, spectrum,
                                        **self.ld_kwargs)

        self.data.emit(ld_results)

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

        return analyser

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
        """Terminate the analysis/acquisition."""
        if self.is_paused:
            self.is_paused = False
        self.is_stopped = True


class QHLine(QFrame):
    """Horizontal line widget."""

    def __init__(self):
        """Initialise."""
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)


# Cliet Code
def main():
    """Start the main function."""
    # Create an instance of QApplication
    app = QApplication(sys.argv)

    app.setStyle("Fusion")

    # Use a palette to switch to dark colors:
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(53, 53, 53))
    palette.setColor(QPalette.WindowText, Qt.white)
    palette.setColor(QPalette.Base, QColor(25, 25, 25))
    palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    palette.setColor(QPalette.ToolTipBase, Qt.black)
    palette.setColor(QPalette.ToolTipText, Qt.white)
    palette.setColor(QPalette.Text, Qt.white)
    palette.setColor(QPalette.Button, QColor(53, 53, 53))
    palette.setColor(QPalette.Active, QPalette.Button, QColor(53, 53, 53))
    palette.setColor(QPalette.ButtonText, Qt.white)
    palette.setColor(QPalette.BrightText, Qt.red)
    palette.setColor(QPalette.Link, QColor(42, 130, 218))
    palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    palette.setColor(QPalette.HighlightedText, Qt.black)
    palette.setColor(QPalette.Disabled, QPalette.ButtonText, Qt.darkGray)
    app.setPalette(palette)

    # Show the GUI
    view = CalcFlux()
    view.show()

    # Execute the main loop
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
    main()
    main()
    main()
