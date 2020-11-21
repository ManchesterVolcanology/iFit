import os
import sys
import yaml
import numpy as np
import pandas as pd
import pyqtgraph as pg
from functools import partial
from scipy.optimize import curve_fit
from scipy.interpolate import griddata
from PyQt5.QtGui import QIcon, QPalette, QColor
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QMainWindow, QWidget, QApplication, QGridLayout,
                             QLabel, QTextEdit, QLineEdit, QPushButton, QFrame,
                             QFileDialog, QScrollArea, QCheckBox, QSplitter,
                             QComboBox, QDoubleSpinBox, QTableWidget,
                             QTableWidgetItem, QTabWidget)

try:
    from .make_ils import super_gaussian
    from .haversine import haversine
except ImportError:
    from make_ils import super_gaussian
    from haversine import haversine


# =============================================================================
# Calculate Flux
# =============================================================================

class CalcFlux(QMainWindow):
    def __init__(self, parent=None):
        super(CalcFlux, self).__init__(parent)

        # Set the window properties
        self.setWindowTitle('Calculate Flux')
        self.statusBar().showMessage('Ready')
        self.setGeometry(40, 40, 1000, 500)
        self.setWindowIcon(QIcon('bin/icon.ico'))

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
        """Create the app widgets"""

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
        """Build the main GUI controls"""

        # Generate main layout
        layout = QGridLayout(self.controlFrame)

        # Create input for iFit output
        layout.addWidget(QLabel('iFit File:'), 0, 0)
        self.so2_path = QLineEdit()
        self.so2_path.setFixedSize(200, 25)
        layout.addWidget(self.so2_path, 0, 1)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.so2_path, 'single',
                                    "Comma Separated (*.csv)"))
        layout.addWidget(btn, 0, 2)

        # Add checkbox to remove bad fits
        self.despike_flag = QCheckBox('Remove\nBad fits?')
        layout.addWidget(self.despike_flag, 0, 3, 1, 2)

        # Create input for GPS intput
        layout.addWidget(QLabel('GPS File:'), 1, 0)
        self.gps_path = QLineEdit()
        self.gps_path.setFixedSize(200, 25)
        layout.addWidget(self.gps_path, 1, 1)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.gps_path, 'single',
                                    "GPS File (*.txt)"))
        layout.addWidget(btn, 1, 2)

        # Make a button to read in the data
        btn = QPushButton('Import')
        btn.setFixedSize(100, 25)
        btn.clicked.connect(self.import_data)
        layout.addWidget(btn, 1, 3, 1, 2)

        # Create input for output
        layout.addWidget(QLabel('Output\nFolder:'), 2, 0)
        self.out_path = QLineEdit()
        self.out_path.setFixedSize(200, 25)
        layout.addWidget(self.out_path, 2, 1)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.out_path, 'folder'))
        layout.addWidget(btn, 2, 2)

        # Create a combobox to hold the pre-saved volcano data
        layout.addWidget(QLabel('Flux Units:'), 2, 3)
        self.flux_units = QComboBox()
        self.flux_units.addItems(['kg/s', 't/day'])
        layout.addWidget(self.flux_units, 2, 4)

# =============================================================================
#   Volcano Data
# =============================================================================

    def _createVolcano(self):
        """Inputs for the volcano data"""

        # Load volcano data
        self.volcano_data = {}
        if os.path.isfile('bin/volcano_data.yaml'):
            with open('bin/volcano_data.yaml', 'r') as ymlfile:
                self.volcano_data = yaml.load(ymlfile, Loader=yaml.FullLoader)

        # Create the layout
        layout = QGridLayout(self.volcanoFrame)

        # Create a combobox to hold the pre-saved volcano data
        layout.addWidget(QLabel('Volcano:'), 0, 0)
        self.volcano = QComboBox()
        self.volcano.addItems(['--select--'] + list(self.volcano_data.keys()))
        self.volcano.currentIndexChanged.connect(self.update_volcano_data)
        layout.addWidget(self.volcano, 0, 1)

        # Create inputs for the volcano latitude
        layout.addWidget(QLabel('Volcano\nLatitude:'), 1, 0)
        self.vlat = QLineEdit()
        layout.addWidget(self.vlat, 1, 1)

        # Create inputs for the volcano longitude
        layout.addWidget(QLabel('Volcano\nLongitutde:'), 2, 0)
        self.vlon = QLineEdit()
        layout.addWidget(self.vlon, 2, 1)

        # Create inputs for the time difference
        layout.addWidget(QLabel('Time\nDifference:'), 3, 0)
        self.tdiff = QLineEdit()
        layout.addWidget(self.tdiff, 3, 1)

        # Create inputs for the wind speed
        layout.addWidget(QLabel('Wind\nSpeed:'), 1, 2)
        self.wind_speed = QDoubleSpinBox()
        self.wind_speed.setRange(0, 100)
        self.wind_speed.setValue(1.0)
        layout.addWidget(self.wind_speed, 1, 3)

        # Create input for wind units
        self.wind_units = QComboBox()
        self.wind_units.addItems(['m/s', 'knots'])
        layout.addWidget(self.wind_units, 1, 4)

        # Create inputs for the wind speed
        layout.addWidget(QLabel('Wind\nError:'), 2, 2)
        self.wind_error = QDoubleSpinBox()
        self.wind_error.setRange(0, 1000)
        layout.addWidget(self.wind_error, 2, 3)

        # Create input for wind units
        self.err_units = QComboBox()
        self.err_units.addItems(['%', 'abs'])
        layout.addWidget(self.err_units, 2, 4)

# =============================================================================
#   Create outputs
# =============================================================================

    def _createOutput(self):
        """Program Outputs"""

        # Generate the layout
        layout = QGridLayout(self.outputFrame)

        # Add a button to calculate the flux
        btn = QPushButton('Calculate Flux')
        btn.setFixedSize(90, 25)
        btn.clicked.connect(self.calc_flux)
        layout.addWidget(btn, 0, 0)

        # Add a table to hold the flux results
        self.fluxTable = QTableWidget()
        self.fluxTable.setColumnCount(2)
        self.fluxTable.setRowCount(0)
        self.fluxTable.setHorizontalHeaderLabels(['Flux', 'Error'])

        layout.addWidget(self.fluxTable, 1, 0)

# =============================================================================
#   Graphs
# =============================================================================

    def _createGraphs(self):
        """Generate the graphs"""

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
        graphwin = pg.GraphicsWindow(show=True)

        # Make the graphs
        self.ax = graphwin.addPlot()

        self.ax.setDownsampling(mode='peak')
        self.ax.setClipToView(True)
        self.ax.showGrid(x=True, y=True)

        # Add axis labels
        self.ax.setLabel('left', 'SO<sub>2</sub> SCD (molec cm<sup>-2</sup>)')
        self.ax.setLabel('bottom', 'Time')

        # Add the graphs to the layout
        g1layout.addWidget(graphwin, 0, 0, 0, 0)

        g2layout = QGridLayout(tab2)
        graphwin = pg.GraphicsWindow(show=True)

        # Make the graphs
        self.mapax = graphwin.addPlot()
        self.mapax.setDownsampling(mode='peak')
        self.mapax.setClipToView(True)
        self.mapax.showGrid(x=True, y=True)
        self.mapax.setAspectLocked(True)

        # Add axis labels
        self.mapax.setLabel('left', 'Latitude')
        self.mapax.setLabel('bottom', 'Longitude')

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
        """Update the volcano data on combobox change"""

        volc = str(self.volcano.currentText())

        if volc != '--select--':
            data = self.volcano_data[volc]
            self.vlat.setText(str(data[0]))
            self.vlon.setText(str(data[1]))
            self.tdiff.setText(str(data[2]))

# =============================================================================
#   Import data
# =============================================================================

    def import_data(self):
        """Import the traverse data"""

        # Read in the SO2 results
        so2_df = pd.read_csv(self.so2_path.text(), parse_dates=['Time'])

        self.so2_time = np.array([t.hour + t.minute/60 + t.second/3600
                                  for t in so2_df['Time']])
        self.so2_scd = so2_df['SO2'].to_numpy()
        self.so2_err = so2_df['SO2_err'].to_numpy()

        if self.despike_flag.isChecked():
            mask = so2_df['fit_quality'] != 1
            self.so2_scd = np.ma.masked_where(mask, self.so2_scd)
            self.so2_err = np.ma.masked_where(mask, self.so2_err)

        # Read in the GPS data
        gps_df = pd.read_table(self.gps_path.text(), sep='\t',
                               parse_dates=['time'])

        self.gps_time = np.array([t.hour + t.minute/60 + t.second/3600
                                  for t in gps_df['time']])
        self.lat = gps_df['latitude'].to_numpy()
        self.lon = gps_df['longitude'].to_numpy()
        self.alt = gps_df['altitude (m)'].to_numpy()

        ploty = self.so2_scd

        if np.nanmax(self.so2_scd) > 1e6:
            order = int(np.ceil(np.log10(np.nanmax(ploty)))) - 1
            ploty = ploty / 10**order
            plot_err = [ploty - self.so2_err/10**order,
                        ploty + self.so2_err/10**order]
            self.ax.setLabel('left', f'Fit value (1e{order})')

        # Update the graph
        self.ax.clear()
        l1 = pg.PlotCurveItem(self.so2_time, plot_err[0], pen=None)
        l2 = pg.PlotCurveItem(self.so2_time, plot_err[1], pen=None)
        pfill = pg.FillBetweenItem(l1, l2, pen=None, brush=self.b0)
        self.ax.addItem(pfill)
        self.ax.plot(self.so2_time, ploty, pen=self.p0)

        # Add the region selector
        self.lr = pg.LinearRegionItem([self.so2_time[0], self.so2_time[-1]])
        self.lr.setZValue(-10)
        self.ax.addItem(self.lr)

        # Create a traverse counter
        self.trav_no = 0

        # Reset the results table
        self.fluxTable.setRowCount(0)

        # Initialise holder arrays for flux results
        self.fluxes = []
        self.errors = []

# =============================================================================
# Calculate Flux
# =============================================================================

    def calc_flux(self):
        """Calculate the flux from the selected traverse"""

        # Pull the relavant data from the GUI
        vlat = float(self.vlat.text())
        vlon = float(self.vlon.text())
        tdiff = float(self.tdiff.text())
        wind_speed = float(self.wind_speed.value())
        wind_error = float(self.wind_error.value())
        flux_units = self.flux_units.currentText()
        wind_units = self.wind_units.currentText()
        err_units = self.err_units.currentText()

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

        # Find the bounds of the selected area
        i0, i1 = self.lr.getRegion()
        idx = np.where(np.logical_and(self.so2_time >= i0,
                                      self.so2_time <= i1))

        # Extract the travese from the SO2 data
        so2_time = self.so2_time[idx]
        so2_scd = self.so2_scd[idx]
        so2_err = self.so2_err[idx]

        # Find the centre of mass of the plume
        cum_so2_scd = np.cumsum(so2_scd)
        peak_idx = np.abs(cum_so2_scd - cum_so2_scd[-1]/2).argmin()

        # Correct for the time diffeence between the SO2 and GPS data
        gps_time = self.gps_time + tdiff

        # Extract the relevant gps data and interpolate onto the so2 grid
        lat = griddata(gps_time, self.lat, so2_time)
        lon = griddata(gps_time, self.lon, so2_time)

        # Calculate the plume bearing
        volc_loc = [vlat, vlon]
        peak_loc = [lat[peak_idx], lon[peak_idx]]
        plume_dist, plume_bearing = haversine(volc_loc, peak_loc)

        # Calculate the distance and bearing of each measurement vector
        vect = [haversine([lat[i-1], lon[i-1]], [lat[i], lon[i]])
                for i in range(len(lat) - 1)]

        # Unpack the distance and bearing from the vectors
        trav_dist, trav_bearing = np.asarray(vect).T

        # Correct the distance for the angle between the travel and plume
        # vectors
        corr_factor = np.sin(plume_bearing-trav_bearing)
        corr_dist = np.multiply(trav_dist, corr_factor)

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
        if flux_units == 't/day':
            flux = abs(so2_kg_per_s * 60*60*24 / 1000.0)
        else:
            flux = abs(so2_kg_per_s)

        # Calculate the Flux Error
        tot_so2_err = np.sum(np.power(so2_err, 2)) ** 0.5
        frac_so2_err = tot_so2_err / np.sum(so2_scd)

        # Combine with the wind speed error
        frac_err = ((frac_so2_err)**2 + (wind_error)**2)**0.5
        flux_err = flux * frac_err

        # Add the flux result to the table
        self.fluxTable.setRowCount(self.fluxTable.rowCount()+1)
        self.fluxTable.setItem(self.trav_no, 0,
                               QTableWidgetItem(f'{flux:.02f}'))
        self.fluxTable.setItem(self.trav_no, 1,
                               QTableWidgetItem(f'{flux_err:.02f}'))

        # Make sure the output directory exists, and create if not
        out_path = self.out_path.text()
        if not os.path.isdir(out_path):
            os.makedirs(out_path)

        # Plot the traverse graphs
        self.mapax.clear()
        self.mapax.plot(self.lon, self.lat, pen=self.p0)
        self.mapax.plot(lon, lat, pen=self.p1)
        self.mapax.plot([vlon], [vlat], pen=None, symbol='o', symbolPen=None,
                        symbolSize=10, symbolBrush=(255, 255, 255))
        self.mapax.plot([peak_loc[1]], [peak_loc[0]], pen=None, symbol='o',
                        symbolPen=None, symbolSize=10,
                        symbolBrush=(0, 255, 0))

        # Collate the output data
        self.fluxes.append(flux)
        self.errors.append(flux_err)
        output_data = np.column_stack([so2_time[1:], so2_scd[1:], so2_err[1:],
                                       lat[1:], lon[1:], trav_dist,
                                       trav_bearing, corr_factor])

        # Write detailed output
        with open(out_path + 'flux_results.csv', 'a') as w:
            w.write(f'Traverse Number,{self.trav_no}\n')
            w.write(f'Flux ({flux_units}),{flux}\n')
            w.write(f'Error ({flux_units}),{flux_err}\n')
            w.write(f'Wind Speed (m/s), {wind_speed}\n')
            w.write(f'Volcano Lat/Lon,{vlat},{vlon}\n')
            w.write(f'Plume Center,{peak_loc[0]},{peak_loc[1]}\n')
            w.write('Time (Decimal hours),SO2 SCD,SO2 Err,Lat,Lon,Distance,'
                    + 'Bearing,Correction Factor\n')

            for i, line in enumerate(output_data):
                w.write(str(line[0]))
                for j in line[1:]:
                    w.write(f',{j}')
                w.write('\n')
            w.write('\n')

        # Write summary output
        w_mean_flux = np.average(self.fluxes,
                                 weights=np.power(self.errors, -2))
        w_mean_error = (1 / np.sum(np.power(self.errors, -2)))**0.5
        with open(out_path + 'flux_summary.txt', 'w') as w:
            w.write(f'Flux Summary\nFlux ({flux_units}), Error\n')
            for i in range(len(self.fluxes)):
                w.write(f'{self.fluxes[i]}, {self.errors[i]}\n')
            w.write(f'Av Flux = {np.average(self.fluxes):.02f} {flux_units}\n')
            w.write(f'Stdev = {np.std(self.fluxes):.02f} {flux_units}\n')
            w.write(f'Weighted Mean Flux = {w_mean_flux:.02f}'
                    + f' (+/- {w_mean_error:.02f}) {flux_units}')

        self.trav_no += 1


# =============================================================================
# Measure ILS
# =============================================================================

class ILSWindow(QMainWindow):
    def __init__(self, parent=None):
        super(ILSWindow, self).__init__(parent)

        # Set the window properties
        self.setWindowTitle('Measure ILS')
        self.statusBar().showMessage('Ready')
        self.setGeometry(40, 40, 1000, 600)
        self.setWindowIcon(QIcon('bin/icon.ico'))

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
        btn.clicked.connect(self.measure_ils)
        layout.addWidget(btn, 3, 2)

        btn = QPushButton('Save')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(self.save_fit)
        layout.addWidget(btn, 3, 3)

        graphwin = pg.GraphicsWindow(show=True)
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
        """Read in ILS spectra and display"""

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
        """Measure the ILS on the selected line"""

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
        w, k, a_w, a_k, shift, amp, offset = self.popt

        fwem = 2*w

        np.savetxt(self.save_path.text(), [fwem, k, a_w, a_k],
                   header='ILS parameters (FWEM, k, a_w, a_k)')


# =============================================================================
# Measure Flat Spectrum
# =============================================================================

class FLATWindow(QMainWindow):
    def __init__(self, parent=None):
        super(FLATWindow, self).__init__(parent)

        # Set the window properties
        self.setWindowTitle('Measure Flat Field')
        self.statusBar().showMessage('Ready')
        self.setGeometry(40, 40, 1000, 500)
        self.setWindowIcon(QIcon('bin/icon.ico'))

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

        graphwin = pg.GraphicsWindow(show=True)
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
        """Read in ILS spectra and display"""
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
        """Measure the flat spectrum across the selected region"""

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

        data = np.column_stack([self.grid, self.flat])
        header = 'Flat spectrum\nWavelength (nm),       Flat Response'
        np.savetxt(self.save_path.text(), data, header=header)


def browse(gui, widget, mode='single', filter=False):

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


# Cliet Code
def main():
    """Main function"""
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
