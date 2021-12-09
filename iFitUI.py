"""Main script for iFit graphical user interface."""
import os
import sys
import yaml
import logging
import numpy as np
import pyqtgraph as pg
from datetime import datetime
from functools import partial
from logging.handlers import RotatingFileHandler
from PyQt5.QtGui import QIcon, QPalette, QColor, QFont
from PyQt5.QtCore import Qt, pyqtSlot, QThread
from PyQt5.QtWidgets import (QMainWindow, QWidget, QApplication, QGridLayout,
                             QMessageBox, QLabel, QComboBox, QTextEdit,
                             QLineEdit, QPushButton, QProgressBar, QFrame,
                             QSplitter, QCheckBox, QSizePolicy, QSpacerItem,
                             QTabWidget, QAction, QFileDialog, QScrollArea,
                             QToolBar, QTableWidget, QHeaderView,
                             QTableWidgetItem)

from ifit.gui_functions import (Widgets, SpinBox, DSpinBox, ParamTable,
                                AnalysisWorker, QTextEditLogger, GPSWorker,
                                AcqScopeWorker, AcqSpecWorker, QHLine, QVLine,
                                browse, GPSWizard)
from ifit.gui_tools import ILSWindow, FLATWindow, CalcFlux, LDFWindow
from ifit.spectrometers import Spectrometer
from ifit.load_spectra import average_spectra
from ifit.gps import GPS

__version__ = '3.4'
__author__ = 'Ben Esse'

# Set up logging
logger = logging.getLogger()
if not os.path.isdir('bin/'):
    os.makedirs('bin/')
fh = RotatingFileHandler('bin/iFit.log', maxBytes=20000, backupCount=5)
fh.setLevel(logging.INFO)
fmt = '%(asctime)s %(levelname)s %(module)s %(funcName)s %(message)s'
fh.setFormatter(logging.Formatter(fmt))
logger.addHandler(fh)


class MainWindow(QMainWindow):
    """View for the iFit GUI."""

    def __init__(self):
        """View initialiser."""
        super().__init__()

        # Set the window properties
        self.setWindowTitle(f'iFit {__version__}')
        self.statusBar().showMessage('Ready')
        self.setGeometry(40, 40, 1210, 700)
        self.setWindowIcon(QIcon('bin/icons/main.ico'))

        # Set the window layout
        self.generalLayout = QGridLayout()
        self._centralWidget = QScrollArea()
        self.widget = QWidget()
        self.setCentralWidget(self._centralWidget)
        self.widget.setLayout(self.generalLayout)

        # Scroll Area Properties
        self._centralWidget.setWidgetResizable(True)
        self._centralWidget.setWidget(self.widget)

        # Setup widget stylesheets
        QTabWidget().setStyleSheet('QTabWidget { font-size: 18pt; }')

        # Create an empty dictionary to hold the GUI widgets
        self.widgets = Widgets()

        # Set the default theme
        self.theme = 'Dark'

        # Build the GUI
        self._createApp()

        # Update widgets from loaded config file
        self.config = {}
        self.config_fname = None
        if os.path.isfile('bin/.config'):
            with open('bin/.config', 'r') as r:
                self.config_fname = r.readline().strip()
            self.load_config(fname=self.config_fname)

        # Update GUI theme
        if self.theme == 'Dark':
            self.changeThemeDark()
        elif self.theme == 'Light':
            self.changeThemeLight()

    def _createApp(self):
        """Build the main GUI."""
        # Generate actions
        # Save action
        saveAct = QAction(QIcon('bin/icons/save.png'), '&Save', self)
        saveAct.setShortcut('Ctrl+S')
        saveAct.triggered.connect(partial(self.save_config, False))

        # Save As action
        saveasAct = QAction(QIcon('bin/icons/saveas.png'), '&Save As', self)
        saveasAct.setShortcut('Ctrl+Shift+S')
        saveasAct.triggered.connect(partial(self.save_config, True))

        # Load action
        loadAct = QAction(QIcon('bin/icons/open.png'), '&Load', self)
        loadAct.triggered.connect(partial(self.load_config, None))

        # Change theme action
        themeAct = QAction(QIcon('bin/icons/theme.png'), '&Change Theme', self)
        themeAct.triggered.connect(self.change_theme)

        # ILS GUI action
        ilsAct = QAction(QIcon('bin/icons/ils.png'), '&Measure ILS', self)
        ilsAct.triggered.connect(self.open_ils_window)
        flatAct = QAction(QIcon('bin/icons/flat.png'), '&Measure\nFlat Field',
                          self)
        # Flat GUI action
        flatAct.triggered.connect(self.open_flat_window)
        fluxAct = QAction(QIcon('bin/icons/flux.png'), '&Calculate flux', self)
        fluxAct.triggered.connect(self.open_flux_window)

        # LDF GUI action
        ldfAct = QAction(QIcon('bin/icons/ldf.png'),
                         '&Light Dilution\nAnalysis', self)
        ldfAct.triggered.connect(self.open_ldf_window)

        # Add menubar
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(saveAct)
        fileMenu.addAction(saveasAct)
        fileMenu.addAction(loadAct)
        toolMenu = menubar.addMenu('&Tools')
        toolMenu.addAction(ilsAct)
        toolMenu.addAction(flatAct)
        toolMenu.addAction(fluxAct)
        toolMenu.addAction(ldfAct)
        toolMenu = menubar.addMenu('&View')
        toolMenu.addAction(themeAct)

        # Create a toolbar
        toolbar = QToolBar("Main toolbar")
        self.addToolBar(toolbar)
        toolbar.addAction(saveAct)
        toolbar.addAction(saveasAct)
        toolbar.addAction(loadAct)
        toolbar.addAction(themeAct)
        toolbar.addSeparator()
        toolbar.addAction(ilsAct)
        toolbar.addAction(flatAct)
        toolbar.addAction(fluxAct)
        toolbar.addAction(ldfAct)

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

        # Generate the GUI widgets
        self._createControls()
        self._createOuput()
        self._createGraphs()

# =============================================================================
#   Generate the program controls
# =============================================================================

    def _createControls(self):
        """Generate the main GUI controls."""
        # Setup tab layout
        tablayout = QGridLayout(self.controlFrame)

        # Generate tabs for the gaphs and settings
        tab1 = QWidget()
        tab2 = QWidget()

        # Form the tab widget
        tabwidget = QTabWidget()
        tabwidget.addTab(tab1, 'Real-Time')
        tabwidget.addTab(tab2, 'Post-Analysis')
        tablayout.addWidget(tabwidget, 0, 0)

# =============================================================================
#       Real time controls
# =============================================================================

        # Setup the layout
        layout = QGridLayout(tab1)
        layout.setAlignment(Qt.AlignTop)
        nrow = 0

        # Add an input for the save selection
        layout.addWidget(QLabel('Output\nFolder:'), nrow, 0)
        self.widgets['rt_save_path'] = QLineEdit('Results')
        self.widgets['rt_save_path'].setToolTip('Folder to hold real time '
                                                + 'results')
        layout.addWidget(self.widgets['rt_save_path'], nrow, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.widgets['rt_save_path'],
                                    'folder', None))
        layout.addWidget(btn, nrow, 4)
        nrow += 1

        layout.addWidget(QHLine(), nrow, 0, 1, 5)
        nrow += 1

        # Set a label for the spectrometer ID
        self.spectro_connected_flag = False
        layout.addWidget(QLabel('Spectrometer:'), nrow, 0)
        self.spec_id = QLabel('Not connected')
        self.spec_id.setToolTip('Spectrometer Serial Number')
        layout.addWidget(self.spec_id, nrow, 1)

        # Create a button to connect to a spectrometer
        self.spec_connect_btn = QPushButton('Connect')
        self.spec_connect_btn.setToolTip(
            'Connect or disconnect the spectrometer')
        self.spec_connect_btn.clicked.connect(self.connect_spectrometer)
        layout.addWidget(self.spec_connect_btn, nrow, 2)

        nrow += 1

        # Create a control for the spectrometer integration time
        layout.addWidget(QLabel('Integration\nTime (ms):'), nrow, 0)
        self.widgets['int_time'] = SpinBox(100, [10, 1000000])
        self.widgets['int_time'].setToolTip('Spectrometer integration time')
        layout.addWidget(self.widgets['int_time'], nrow, 1)

        # Create a button to update the integration time
        self.update_inttime_btn = QPushButton('Update')
        self.update_inttime_btn.setToolTip('Update Spectrometer Integration'
                                           + ' Time')
        self.update_inttime_btn.clicked.connect(self.update_int_time)
        self.update_inttime_btn.setEnabled(False)
        layout.addWidget(self.update_inttime_btn, nrow, 2)

        # Add stereo button for non-liniarity correction
        self.widgets['nonlin_flag'] = QCheckBox('Correct\nNon-Linearity?')
        self.widgets['nonlin_flag'].setToolTip('Turn on correction for non-'
                                               + 'linear intensity response'
                                               + ' correction')
        layout.addWidget(self.widgets['nonlin_flag'], nrow, 3)

        nrow += 1

        # Create a control for the spectrometer coadds
        layout.addWidget(QLabel('Coadds:'), nrow, 0)
        self.widgets['coadds'] = SpinBox(10, [1, 1000000])
        self.widgets['coadds'].setToolTip('No. spectra to average')
        layout.addWidget(self.widgets['coadds'], nrow, 1)

        # Create a button to update the coadds
        self.update_coadds_btn = QPushButton('Update')
        self.update_coadds_btn.setToolTip('Update Spectrometer Coadds')
        self.update_coadds_btn.clicked.connect(self.update_coadds)
        self.update_coadds_btn.setEnabled(False)
        layout.addWidget(self.update_coadds_btn, nrow, 2)

        # Add stereo button for non-liniarity correction
        self.widgets['eldark_flag'] = QCheckBox('Correct\nElectronic dark?')
        self.widgets['eldark_flag'].setToolTip('Turn on electronic dark '
                                               + 'correction')
        layout.addWidget(self.widgets['eldark_flag'], nrow, 3)

        nrow += 1

        # Create a control for the number of dark spectra
        layout.addWidget(QLabel('No. Dark\nSpectra:'), nrow, 0)
        self.widgets['ndarks'] = SpinBox(10, [1, 1000000])
        self.widgets['ndarks'].setToolTip('Set number of dark spectra to '
                                          + 'measure')
        layout.addWidget(self.widgets['ndarks'], nrow, 1)

        # Create a button to acquire the dark spectra
        self.acquire_darks_btn = QPushButton('Acquire')
        self.acquire_darks_btn.setToolTip('Measure Dark Spectra')
        self.acquire_darks_btn.clicked.connect(partial(self.begin_acquisition,
                                                       'acquire_darks'))
        self.acquire_darks_btn.setEnabled(False)
        layout.addWidget(self.acquire_darks_btn, nrow, 2)

        # Create a button to toggle real-time analysis
        self.rt_fitting_flag = False
        self.rt_flag_btn = QPushButton('Fitting OFF')
        self.rt_flag_btn.setToolTip('Toggle Real Time Fitting')
        self.rt_flag_btn.clicked.connect(self.toggle_fitting)
        self.rt_flag_btn.setEnabled(False)
        self.rt_flag_btn.setStyleSheet("background-color: red")
        layout.addWidget(self.rt_flag_btn, nrow, 3)

        nrow += 1

        layout.addWidget(QHLine(), nrow, 0, 1, 5)
        nrow += 1

        # Add GPS connection settings
        self.gps_connected_flag = False
        self.gps = None
        layout.addWidget(QLabel('GPS Status:'), nrow, 0)
        self.gps_status = QLabel('Not connected')
        layout.addWidget(self.gps_status, nrow, 1)
        self.gps_connect_btn = QPushButton('Connect')
        self.gps_connect_btn.clicked.connect(self.connect_gps)
        layout.addWidget(self.gps_connect_btn, nrow, 2)
        nrow += 1

        layout.addWidget(QHLine(), nrow, 0, 1, 5)
        nrow += 1

        # Add button to begin analysis
        self.rt_start_btn = QPushButton('Begin!')
        self.rt_start_btn.setToolTip('Begin Spectra Acquisition')
        self.rt_start_btn.clicked.connect(partial(self.begin_acquisition,
                                                  'acquire_cont'))
        self.rt_start_btn.setFixedSize(90, 25)
        self.rt_start_btn.setEnabled(False)
        layout.addWidget(self.rt_start_btn, nrow, 1)

        # Add button to pause analysis
        self.rt_pause_btn = QPushButton('Pause')
        self.rt_pause_btn.setToolTip('Pause/Play Spectra Acquisition')
        self.rt_pause_btn.clicked.connect(partial(self.pause))
        self.rt_pause_btn.setFixedSize(90, 25)
        self.rt_pause_btn.setEnabled(False)
        layout.addWidget(self.rt_pause_btn, nrow, 2)

        # Add button to stop analysis
        self.rt_stop_btn = QPushButton('Stop')
        self.rt_stop_btn.setToolTip('Stop Spectra Acquisition')
        self.rt_stop_btn.clicked.connect(partial(self.stop))
        self.rt_stop_btn.setFixedSize(90, 25)
        self.rt_stop_btn.setEnabled(False)
        layout.addWidget(self.rt_stop_btn, nrow, 3)

# =============================================================================
#       Post-procesing controls
# =============================================================================

        # Setup the layout
        layout = QGridLayout(tab2)
        layout.setAlignment(Qt.AlignTop)

        # Create an option menu for the spectra format
        layout.addWidget(QLabel('Format:'), 0, 0)
        self.widgets['spec_type'] = QComboBox()
        self.widgets['spec_type'].setToolTip('Choose spectrum format')
        self.widgets['spec_type'].addItems(['iFit',
                                            'iFit (old)',
                                            'Master.Scope',
                                            'Spectrasuite',
                                            'mobileDOAS',
                                            'Basic'])
        self.widgets['spec_type'].setFixedSize(100, 20)
        layout.addWidget(self.widgets['spec_type'], 0, 1)

        # Add an input for the spectra selection
        layout.addWidget(QLabel('Spectra:'), 1, 0)
        self.widgets['spec_fnames'] = QTextEdit()
        self.widgets['spec_fnames'].setToolTip('Measurement spectrum files')
        layout.addWidget(self.widgets['spec_fnames'], 1, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.widgets['spec_fnames'],
                                    'multi', None))
        layout.addWidget(btn, 1, 4)

        # Add an input for the dark selection
        layout.addWidget(QLabel('Dark\nSpectra:'), 2, 0)
        self.widgets['dark_fnames'] = QTextEdit()
        self.widgets['dark_fnames'].setToolTip('Dark spectrum files')
        layout.addWidget(self.widgets['dark_fnames'], 2, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.widgets['dark_fnames'],
                                    'multi', None))
        layout.addWidget(btn, 2, 4)

        # Add an input for the save selection
        layout.addWidget(QLabel('Output\nFile:'), 3, 0)
        self.widgets['save_path'] = QLineEdit()
        self.widgets['save_path'].setToolTip('Output .csv file')
        layout.addWidget(self.widgets['save_path'], 3, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.widgets['save_path'],
                                    'save', "Comma Separated (*.csv)"))
        layout.addWidget(btn, 3, 4)

        # Add button to begin analysis
        self.start_btn = QPushButton('Begin!')
        self.start_btn.setToolTip('Begin spectra analysis')
        self.start_btn.clicked.connect(partial(self.begin_analysis,
                                               'post_analyse'))
        self.start_btn.setFixedSize(90, 25)
        layout.addWidget(self.start_btn, 4, 1)

        # Add button to pause analysis
        self.pause_btn = QPushButton('Pause')
        self.pause_btn.setToolTip('Pause/play spectra analysis')
        self.pause_btn.clicked.connect(partial(self.pause))
        self.pause_btn.setFixedSize(90, 25)
        self.pause_btn.setEnabled(False)
        layout.addWidget(self.pause_btn, 4, 2)

        # Add button to stop analysis
        self.stop_btn = QPushButton('Stop')
        self.stop_btn.setToolTip('Stop spectra analysis')
        self.stop_btn.clicked.connect(partial(self.stop))
        self.stop_btn.setFixedSize(90, 25)
        self.stop_btn.setEnabled(False)
        layout.addWidget(self.stop_btn, 4, 3)

# =============================================================================
#   Generate the program outputs
# =============================================================================

    def _createOuput(self):
        """Build the main GUI visual ouputs."""
        layout = QGridLayout(self.outputFrame)
        layout.setAlignment(Qt.AlignTop)

        # Add a progress bar
        self.progress = QProgressBar(self)
        self.progress.setFixedSize(480, 25)
        layout.addWidget(self.progress, 0, 0, 1, 6)

        # Add spectrum number indicator
        layout.addWidget(QLabel('Spectrum Number:'), 1, 0)
        self.spectrum_number = QLabel('-')
        layout.addWidget(self.spectrum_number, 1, 1)

        # Add numerical outputs
        layout.addWidget(QLabel('Last amt:'), 2, 0)
        self.last_amt = QLabel('-')
        layout.addWidget(self.last_amt, 2, 1)
        layout.addWidget(QLabel('+/-'), 2, 2)
        self.last_err = QLabel('-')
        layout.addWidget(self.last_err, 2, 3)

        # Create a textbox to display the program logs
        self.logBox = QTextEditLogger(self)
        fmt = logging.Formatter('%(asctime)s - %(message)s', '%H:%M:%S')
        self.logBox.setFormatter(fmt)
        logger.addHandler(self.logBox)
        logger.setLevel(logging.INFO)
        layout.addWidget(self.logBox.widget, 3, 0, 1, 6)
        msg = 'Welcome to iFit! Written by Ben Esse'
        self.logBox.widget.appendPlainText(msg)

# =============================================================================
#   Set up graphs and settings
# =============================================================================

    def _createGraphs(self):
        """Build the graphical display and program settings."""
        layout = QGridLayout(self.graphFrame)

        # Generate tabs for the graphs and settings
        tab1 = QWidget()
        tab2 = QWidget()
        tab3 = QWidget()

        # Form the tab widget
        tabwidget = QTabWidget()
        tabwidget.addTab(tab1, 'Analysis')
        tabwidget.addTab(tab3, 'Scope')
        tabwidget.addTab(tab2, 'Settings')
        layout.addWidget(tabwidget, 0, 0)

        # Make analysis tab holder
        analysis_tabs = QTabWidget()
        graphtab = QWidget()
        tabletab = QWidget()
        maptab = QWidget()
        analysis_tabs.addTab(graphtab, 'Graphs')
        analysis_tabs.addTab(tabletab, 'Table')
        analysis_tabs.addTab(maptab, 'Map')

        analysis_layout = QGridLayout(tab1)
        analysis_layout.addWidget(analysis_tabs, 0, 0, 1, 8)

# =============================================================================
#       Set up the analysis graphs
# =============================================================================

        # Generate the graph window
        self.graphwin = pg.GraphicsLayoutWidget(show=True)
        pg.setConfigOptions(antialias=True)
        # pg.setConfigOptions(useOpenGL=True)

        # Make the graphs
        ax0 = self.graphwin.addPlot(row=0, col=0)
        ax1 = self.graphwin.addPlot(row=0, col=1)
        ax2 = self.graphwin.addPlot(row=1, col=0)
        ax3 = self.graphwin.addPlot(row=1, col=1)
        ax4 = self.graphwin.addPlot(row=2, col=0, colspan=2)
        self.plot_axes = [ax0, ax1, ax2, ax3, ax4]

        for ax in self.plot_axes:
            ax.setDownsampling(mode='peak')
            ax.setClipToView(True)
            ax.showGrid(x=True, y=True)

        # Add axis labels
        ax0.setLabel('left', 'Intensity (counts)')
        ax1.setLabel('left', 'Intensity (counts)')
        ax2.setLabel('left', 'Residual')
        ax3.setLabel('left', 'Optical Depth')
        ax2.setLabel('bottom', 'Wavelength (nm)')
        ax3.setLabel('bottom', 'Wavelength (nm)')
        ax4.setLabel('left', 'Fitted Value')
        ax4.setLabel('bottom', 'Spectrum Number')

        # Initialise the lines
        p0 = pg.mkPen(color='#1f77b4', width=1.5)
        p1 = pg.mkPen(color='#ff7f0e', width=1.0)
        l0 = ax0.plot(pen=p0)
        l1 = ax0.plot(pen=p1)
        l2 = ax1.plot(pen=p0)
        l3 = ax2.plot(pen=p0)
        l4 = ax3.plot(pen=p0)
        l5 = ax3.plot(pen=p1)
        l6 = ax4.plot(pen=pg.mkPen(color='#1f77b4', width=2.0))

        # Initialise the error plot
        self.err_plot = pg.ErrorBarItem(pen=p0)
        ax4.addItem(self.err_plot)

        # Add legend to first axis
        legend = ax0.addLegend()
        legend.addItem(l0, 'Spectrum')
        legend.addItem(l1, 'Fit')

        # Colate lines into a dictionary
        self.plot_lines = {'spectrum': l0,
                           'fit': l1,
                           'measurement': l2,
                           'residual': l3,
                           'meas_od': l4,
                           'fit_od': l5,
                           'series': l6}

        # Add the graphs to the layout
        glayout = QGridLayout(graphtab)
        glayout.addWidget(self.graphwin, 0, 0)

# =============================================================================
#      Analysis results table
# =============================================================================

        # Generate results table
        tlayout = QGridLayout(tabletab)
        self.results_table = QTableWidget(0, 4)
        self.results_table.setHorizontalHeaderLabels(
            ['Parameter', 'Vary?', 'Fit Value', 'Fit Error'])
        self.results_table.horizontalHeader().setStretchLastSection(True)
        self.results_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch)

        tlayout.addWidget(self.results_table, 0, 0)

# =============================================================================
#      Analysis results map
# =============================================================================

        # Make the layout
        mlayout = QGridLayout(maptab)

        # Add markers for the current lat/lon/alt values
        self.gps_timestamp = QLabel('-')
        self.gps_lat = QLabel('-')
        self.gps_lon = QLabel('-')
        self.gps_alt = QLabel('-')
        mlayout.addWidget(QLabel('GPS Time:'), 0, 0)
        mlayout.addWidget(self.gps_timestamp, 0, 1)
        mlayout.addWidget(QLabel('Latitude:'), 0, 2)
        mlayout.addWidget(self.gps_lat, 0, 3)
        mlayout.addWidget(QLabel('Longitude:'), 0, 4)
        mlayout.addWidget(self.gps_lon, 0, 5)
        mlayout.addWidget(QLabel('Altitude:'), 0, 6)
        mlayout.addWidget(self.gps_alt, 0, 7)

        # Add checkbox to control if the colorlimit is fixed or automatic
        self.widgets['auto_map_scale'] = QCheckBox('Auto Scale Map?')
        mlayout.addWidget(self.widgets['auto_map_scale'], 0, 8)

        # Make the graphics window
        self.mapwin = pg.GraphicsLayoutWidget(show=True)
        mlayout.addWidget(self.mapwin, 2, 0, 1, 9)

        # Generate the axes
        self.map_ax = self.mapwin.addPlot(row=0, col=0)
        self.map_ax.setDownsampling(mode='peak')
        self.map_ax.setClipToView(True)
        self.map_ax.showGrid(x=True, y=True)
        self.map_ax.setAspectLocked(True)

        # Generate the colorbar
        self.cmap = pg.colormap.get('viridis')
        im = pg.ImageItem()
        self.cbar = pg.ColorBarItem(values=(0, 100), colorMap=self.cmap,
                                    label='Fit Value')
        self.cbar.setImageItem(im)
        self.cbar.sigLevelsChangeFinished.connect(self.update_map)
        self.mapwin.addItem(self.cbar, 0, 2)

        # Add the axis labels
        self.map_ax.setLabel('left', 'Latitude')
        self.map_ax.setLabel('bottom', 'Longitude')

        # Initialise the plots for the GPS data
        self.gps_line = self.map_ax.plot(pen=p0)
        self.gps_scatter = pg.ScatterPlotItem()
        self.map_ax.addItem(self.gps_scatter)

# =============================================================================
#      Graph settings
# =============================================================================

        # Create a checkbox to turn plotting on or off
        self.widgets['graph_flag'] = QCheckBox('Update Graphs?')
        self.widgets['graph_flag'].setChecked(True)
        self.widgets['graph_flag'].setToolTip('Display graph plots (can slow '
                                              + 'analysis)')
        analysis_layout.addWidget(self.widgets['graph_flag'], 1, 0)

        # Create a checkbox to only display good fits
        self.widgets['good_fit_flag'] = QCheckBox('Only Show\nGood Fits?')
        self.widgets['good_fit_flag'].setToolTip('Only display results for '
                                                 + 'fits that pass the '
                                                 + 'quality checks')
        analysis_layout.addWidget(self.widgets['good_fit_flag'], 1, 1)

        # Add combo box for the graph parameter
        analysis_layout.addWidget(QLabel('Target\nParameter:'), 1, 2)
        self.widgets['graph_param'] = QComboBox()
        self.widgets['graph_param'].addItems([''])
        analysis_layout.addWidget(self.widgets['graph_param'], 1, 3)

        # Create a checkbox to plot the fit error
        self.widgets['graph_err_flag'] = QCheckBox('Show error?')
        self.widgets['graph_err_flag'].setToolTip('Plot fit errors for the '
                                                  + 'tarjet parameter')
        analysis_layout.addWidget(self.widgets['graph_err_flag'], 1, 4)

        # Create a checkbox to turn scrolling on or off
        self.widgets['scroll_flag'] = QCheckBox('Scroll Graphs?')
        self.widgets['scroll_flag'].setToolTip('Allow graphs to scroll\n'
                                               + '(limits no. spectra '
                                               + 'displayed)')
        analysis_layout.addWidget(self.widgets['scroll_flag'], 1, 5)

        # Add spinbox for the graph scroll amount
        analysis_layout.addWidget(QLabel('No. Spectra\nTo Display:'), 1, 6)
        self.widgets['scroll_amt'] = SpinBox(100, [1, 10000])
        # self.widgets['scroll_amt'].setFixedSize(70, 20)
        analysis_layout.addWidget(self.widgets['scroll_amt'], 1, 7)

# =============================================================================
#      Set up the scope plot
# =============================================================================

        self.scopewin = pg.GraphicsLayoutWidget(show=True)
        slayout = QGridLayout(tab3)

        # Make the graph
        self.scope_ax = self.scopewin.addPlot(row=0, col=0)
        self.scope_ax.setDownsampling(mode='peak')
        self.scope_ax.setClipToView(True)
        self.scope_ax.showGrid(x=True, y=True)

        self.scope_line = self.scope_ax.plot([], [], pen=p0)

        # Add the graphs to the layout
        slayout.addWidget(self.scopewin, 0, 0)

# =============================================================================
#       Create settings
# =============================================================================

        # Create tabs for settings
        slayout = QGridLayout(tab2)

        stab1 = QWidget()
        stab2 = QWidget()
        stab3 = QWidget()
        stab4 = QWidget()

        tabwidget = QTabWidget()
        tabwidget.addTab(stab1, 'Model')
        tabwidget.addTab(stab2, 'Spectrometer')
        tabwidget.addTab(stab3, 'Parameters')
        tabwidget.addTab(stab4, 'Outputs')
        slayout.addWidget(tabwidget, 0, 0)

# =============================================================================
#       Model Settings
# =============================================================================

        # Setup the layout
        layout = QGridLayout(stab1)
        layout.setAlignment(Qt.AlignTop)
        nrow = 1
        ncol = 0

        # Add column header
        header = QLabel('Model Settings')
        header.setAlignment(Qt.AlignCenter)
        header.setFont(QFont('Arial', 14))
        layout.addWidget(header, 0, ncol, 1, 2)

        # Add spinboxs for the fit window
        layout.addWidget(QLabel('Fit Window:\n    (nm)'), nrow, ncol, 2, 1)
        self.widgets['fit_lo'] = DSpinBox(310, [0, 10000])
        self.widgets['fit_lo'].setToolTip('Fit window lower limit')
        self.widgets['fit_lo'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['fit_lo'], nrow, ncol+1)
        nrow += 1
        self.widgets['fit_hi'] = DSpinBox(320, [0, 10000])
        self.widgets['fit_hi'].setToolTip('Fit window upper limit')
        self.widgets['fit_hi'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['fit_hi'], nrow, ncol+1)
        nrow += 1

        # Add spinbox for the model grid padding
        layout.addWidget(QLabel('Model Grid\nPadding (nm):'), nrow, ncol)
        self.widgets['model_padding'] = DSpinBox(1.0, [0, 10000])
        self.widgets['model_padding'].setToolTip('Pad the analysis window to '
                                                 + 'avoid edge effects in fit')
        self.widgets['model_padding'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['model_padding'], nrow, ncol+1)
        nrow += 1

        # Add spinbox for the model grid spacing
        layout.addWidget(QLabel('Model Grid\nSpacing (nm):'), nrow, ncol)
        self.widgets['model_spacing'] = DSpinBox(0.01, [0, 10])
        self.widgets['model_spacing'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['model_spacing'], nrow, ncol+1)
        nrow += 1

        # Add combo box for interpolation method
        layout.addWidget(QLabel('Interpolation\nMethod:'), nrow, ncol)
        self.widgets['interp_method'] = QComboBox()
        self.widgets['interp_method'].setToolTip('How the model spectrum is '
                                                 + 'interpolated onto the '
                                                 + 'model grid')
        self.widgets['interp_method'].addItems(['cubic', 'linear'])
        self.widgets['interp_method'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['interp_method'], nrow, ncol+1)

        # New column
        layout.addWidget(QVLine(), 0, ncol+2, 10, 1)
        nrow = 1
        ncol += 3

        # Add column header
        header = QLabel('Spectrum Preprocessing')
        header.setAlignment(Qt.AlignCenter)
        header.setFont(QFont('Ariel', 14))
        layout.addWidget(header, 0, ncol, 1, 3)

        # Add sterio button for dark correction
        self.widgets['dark_flag'] = QCheckBox('Correct Dark\nSpectrum?')
        self.widgets['dark_flag'].setToolTip('Remove averaged dark spectrum '
                                             + 'before fitting')
        layout.addWidget(self.widgets['dark_flag'], nrow, ncol, 1, 2)
        # nrow += 1

        # Add sterio button for flat correction
        self.widgets['flat_flag'] = QCheckBox('Correct Flat\nSpectrum?')
        self.widgets['flat_flag'].setToolTip('Remove flat-field spectrum '
                                             + 'before fitting')
        layout.addWidget(self.widgets['flat_flag'], nrow, ncol+2, 1, 2)
        nrow += 1

        # Add option to apply pre-fit wavelength shift
        layout.addWidget(QLabel('Pre-Fit\nShift (nm):'), nrow, ncol)
        self.widgets['prefit_shift'] = DSpinBox(0.0, [-1000, 1000], 0.1)
        self.widgets['prefit_shift'].setFixedSize(70, 20)
        nrow += 1

        # Add spinboxs for the stray light window
        layout.addWidget(QLabel('Stray Light\nWindow (nm):'), nrow, ncol, 2, 1)
        self.widgets['stray_lo'] = DSpinBox(280, [0, 10000])
        self.widgets['stray_lo'].setToolTip('Stray window lower limit')
        self.widgets['stray_lo'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['stray_lo'], nrow, ncol+1)
        nrow += 1
        self.widgets['stray_hi'] = DSpinBox(290, [0, 10000])
        self.widgets['stray_hi'].setToolTip('Stray window upper limit')
        self.widgets['stray_hi'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['stray_hi'], nrow, ncol+1)
        nrow += 1

        # Add sterio button for stray light correction
        self.widgets['stray_flag'] = QCheckBox('Remove?')
        self.widgets['stray_flag'].setToolTip('Apply stray light correction '
                                              + 'before fitting')
        layout.addWidget(self.widgets['stray_flag'], nrow-2, ncol+2, 2, 1)

        # Add spinbox to control spike removal
        layout.addWidget(QLabel('Spike Limit\n(counts):'), nrow, ncol)
        self.widgets['spike_limit'] = SpinBox(1000, [0, 10000000])
        self.widgets['spike_limit'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['spike_limit'], nrow, ncol+1)
        self.widgets['despike_flag'] = QCheckBox('Remove?')
        self.widgets['despike_flag'].setToolTip('Remove spikey pixels, based '
                                                + 'on the Spike Limit')
        layout.addWidget(self.widgets['despike_flag'], nrow, ncol+2)
        nrow += 1

        # New column
        layout.addWidget(QVLine(), 0, ncol+3, 10, 1)
        nrow = 1
        ncol += 4

        # Add column header
        header = QLabel('Quality Control')
        header.setAlignment(Qt.AlignCenter)
        header.setFont(QFont('Ariel', 14))
        layout.addWidget(header, 0, ncol, 1, 2)

        # Add combo box for residual display
        layout.addWidget(QLabel('Residual Display:'), nrow, ncol)
        self.widgets['resid_type'] = QComboBox()
        self.widgets['resid_type'].addItems(['Percentage', 'Absolute'])
        self.widgets['resid_type'].setFixedSize(100, 20)
        layout.addWidget(self.widgets['resid_type'], nrow, ncol+1)
        nrow += 1

        # Add sterio button for auto-update of fit params
        self.widgets['update_flag'] = QCheckBox('Auto-update\nFit Parameters?')
        self.widgets['update_flag'].setToolTip('Use previous fit results as '
                                               + 'next first guess\n(if fit '
                                               + 'is good)')
        layout.addWidget(self.widgets['update_flag'], nrow, ncol, 1, 2)
        nrow += 1

        # Add spinbox for the residual limit
        layout.addWidget(QLabel('Residual Limit:'), nrow, ncol)
        self.widgets['resid_limit'] = DSpinBox(10.0, [0, 10000])
        self.widgets['resid_limit'].setToolTip('Residual limit for a good fit')
        self.widgets['resid_limit'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['resid_limit'], nrow, ncol+1)
        nrow += 1

        # Add spinboxs for the intensity limits
        layout.addWidget(QLabel('Intensity Limits:'), nrow, ncol, 2, 1)
        self.widgets['lo_int_limit'] = SpinBox(0, [0, 100000])
        self.widgets['lo_int_limit'].setToolTip('Low intensity limit for a '
                                                + 'good fit')
        self.widgets['lo_int_limit'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['lo_int_limit'], nrow, ncol+1)
        nrow += 1
        self.widgets['hi_int_limit'] = SpinBox(70000, [0, 100000])
        self.widgets['hi_int_limit'].setToolTip('High intensity limit for a '
                                                + 'good fit')
        self.widgets['hi_int_limit'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['hi_int_limit'], nrow, ncol+1)
        nrow += 1

        vspacer = QSpacerItem(QSizePolicy.Minimum, QSizePolicy.Expanding)
        layout.addItem(vspacer, nrow, 0, 1, -1)

        hspacer = QSpacerItem(QSizePolicy.Expanding, QSizePolicy.Minimum)
        layout.addItem(hspacer, 0, 10, -1, 1)

# =============================================================================
#       Spectrometer Settings
# =============================================================================

        # Setup the layout
        layout = QGridLayout(stab2)
        layout.setAlignment(Qt.AlignTop)
        nrow = 0
        ncol = 0

        # Create entries for the ILS input
        layout.addWidget(QLabel('Generate ILS:'), nrow, ncol)
        self.widgets['ils_mode'] = QComboBox()
        self.widgets['ils_mode'].addItems(['Manual', 'Params', 'File'])
        self.widgets['ils_mode'].setFixedSize(100, 20)
        layout.addWidget(self.widgets['ils_mode'], nrow, ncol+1)
        nrow += 1

        # Add an input for the ILS parameters
        layout.addWidget(QLabel('ILS Parameter\nFile:'), nrow, ncol)
        self.widgets['ils_path'] = QLineEdit()
        self.widgets['ils_path'].setFixedSize(300, 25)
        layout.addWidget(self.widgets['ils_path'], nrow, ncol+1, 1, 2)
        btn = QPushButton('Browse')
        btn.setFixedSize(100, 25)
        btn.clicked.connect(partial(browse, self, self.widgets['ils_path'],
                                    'single', None))
        layout.addWidget(btn, nrow, ncol+3)
        nrow += 1

        # Add an input for the flat spectrum
        layout.addWidget(QLabel('Flat Spectrum:'), nrow, ncol)
        self.widgets['flat_path'] = QLineEdit()
        self.widgets['flat_path'].setFixedSize(300, 25)
        layout.addWidget(self.widgets['flat_path'], nrow, ncol+1, 1, 2)
        btn = QPushButton('Browse')
        btn.setFixedSize(100, 25)
        btn.clicked.connect(partial(browse, self, self.widgets['flat_path'],
                                    'single', None))
        layout.addWidget(btn, nrow, ncol+3)
        nrow += 1

        # Add an input for the wavelength calibration
        layout.addWidget(QLabel('Wavelength\nCalibration:'), nrow, ncol)
        self.widgets['wl_calib'] = QLineEdit()
        self.widgets['wl_calib'].setFixedSize(300, 25)
        layout.addWidget(self.widgets['wl_calib'], nrow, ncol+1, 1, 2)
        btn = QPushButton('Browse')
        btn.setFixedSize(100, 25)
        btn.clicked.connect(partial(browse, self, self.widgets['wl_calib'],
                                    'single', None))
        layout.addWidget(btn, nrow, ncol+3)
        nrow += 1

        layout.addWidget(QHLine(), nrow, 0, 1, 10)
        nrow += 1

        # Add inputs for the manual ILS parameters
        layout.addWidget(QLabel('FWEM:'), nrow, ncol)
        self.widgets['fwem'] = QLineEdit("0.6")
        self.widgets['fwem'].setFixedSize(100, 20)
        layout.addWidget(self.widgets['fwem'], nrow, ncol+1)
        self.widgets['fwem_fit'] = QCheckBox('Fit?')
        self.widgets['fwem_fit'].setChecked(True)
        layout.addWidget(self.widgets['fwem_fit'], nrow, ncol+2)
        nrow += 1

        layout.addWidget(QLabel('k:'), nrow, ncol)
        self.widgets['k'] = QLineEdit("2.0")
        self.widgets['k'].setFixedSize(100, 20)
        layout.addWidget(self.widgets['k'], nrow, ncol+1)
        self.widgets['k_fit'] = QCheckBox('Fit?')
        self.widgets['k_fit'].setChecked(True)
        layout.addWidget(self.widgets['k_fit'], nrow, ncol+2)
        nrow += 1

        layout.addWidget(QLabel('a<sub>w</sub>:'), nrow, ncol)
        self.widgets['a_w'] = QLineEdit("0.0")
        self.widgets['a_w'].setFixedSize(100, 20)
        layout.addWidget(self.widgets['a_w'], nrow, ncol+1)
        self.widgets['a_w_fit'] = QCheckBox('Fit?')
        self.widgets['a_w_fit'].setChecked(True)
        layout.addWidget(self.widgets['a_w_fit'], nrow, ncol+2)
        nrow += 1

        layout.addWidget(QLabel('a<sub>k</sub>:'), nrow, ncol)
        self.widgets['a_k'] = QLineEdit("0.0")
        self.widgets['a_k'].setFixedSize(100, 20)
        layout.addWidget(self.widgets['a_k'], nrow, ncol+1)
        self.widgets['a_k_fit'] = QCheckBox('Fit?')
        self.widgets['a_k_fit'].setChecked(True)
        layout.addWidget(self.widgets['a_k_fit'], nrow, ncol+2)
        nrow += 1

        layout.addWidget(QHLine(), nrow, 0, 1, 10)
        nrow += 1

        # Add an input for the bad pixel numbers
        layout.addWidget(QLabel('Bad pixels:'), nrow, ncol)
        self.widgets['bad_pixels'] = QLineEdit()
        self.widgets['bad_pixels'].setFixedSize(300, 20)
        layout.addWidget(self.widgets['bad_pixels'], nrow, ncol+1, 1, 2)

        vspacer = QSpacerItem(QSizePolicy.Minimum, QSizePolicy.Expanding)
        layout.addItem(vspacer, nrow, 0, 1, -1)

        hspacer = QSpacerItem(QSizePolicy.Expanding, QSizePolicy.Minimum)
        layout.addItem(hspacer, 0, 4, -1, 1)

# =============================================================================
#       Parameter Settings
# =============================================================================

        # Setup the layout
        layout = QGridLayout(stab3)
        layout.setAlignment(Qt.AlignTop)
        nrow = 0
        ncol = 0

        # Add an input for the FRS file
        layout.addWidget(QLabel('FRS File:'), nrow, ncol)
        self.widgets['frs_path'] = QLineEdit()
        self.widgets['frs_path'].setFixedSize(300, 25)
        layout.addWidget(self.widgets['frs_path'], nrow, ncol+1)
        btn = QPushButton('Browse')
        btn.setFixedSize(100, 25)
        btn.clicked.connect(partial(browse, self, self.widgets['frs_path'],
                                    'single', None))
        layout.addWidget(btn, nrow, ncol+2)
        nrow += 1

        # Create the tabs for the parameters
        ptab1 = QWidget()
        ptab2 = QWidget()
        ptab3 = QWidget()
        ptab4 = QWidget()
        ptab5 = QWidget()

        tabwidget = QTabWidget()
        tabwidget.setFixedWidth(620)
        tabwidget.addTab(ptab1, 'Absorbers')
        tabwidget.addTab(ptab2, 'Polynomial')
        tabwidget.addTab(ptab3, 'Offset')
        tabwidget.addTab(ptab4, 'Shift')
        tabwidget.addTab(ptab5, 'Light Dilution')
        layout.addWidget(tabwidget, 1, 0, 1, 6)

        vspacer = QSpacerItem(QSizePolicy.Minimum, QSizePolicy.Expanding)
        layout.addItem(vspacer, nrow, 0, 1, -1)

        hspacer = QSpacerItem(QSizePolicy.Expanding, QSizePolicy.Minimum)
        layout.addItem(hspacer, 0, 5, -1, 1)

        # Create the absorber and polynomial tables
        self.gas_table = ParamTable(ptab1, 'param', 620, 400)
        self.bgpoly_table = ParamTable(ptab2, 'poly', 550, 400, 'bg_poly')
        self.offset_table = ParamTable(ptab3, 'poly', 550, 400, 'offset')
        self.shift_table = ParamTable(ptab4, 'poly', 550, 400, 'shift')

        # Link the parameter table to the plot parameter combobox
        self.gas_table.cellChanged.connect(self.update_plot_params)

        # Add controls for light dilution
        ld_layout = QGridLayout(ptab5)
        ld_layout.setAlignment(Qt.AlignTop)

        # Add manual or fitted LDF controls
        ld_layout.addWidget(QLabel('LDF:'), 1, 0)
        self.widgets['ldf'] = DSpinBox(0, [0, 1], 0.1)
        self.widgets['ldf'].setFixedSize(100, 20)
        ld_layout.addWidget(self.widgets['ldf'], 1, 1)
        self.widgets['ldf_fit'] = QCheckBox('Fit?')
        ld_layout.addWidget(self.widgets['ldf_fit'], 1, 2)

# =============================================================================
#   Output Settings
# =============================================================================

        # Setup the layout
        layout = QGridLayout(stab4)
        layout.setAlignment(Qt.AlignTop)
        nrow = 0
        ncol = 0

        # Add a header
        header = QLabel('Output File Parameters')
        header.setAlignment(Qt.AlignLeft)
        header.setFont(QFont('Ariel', 12))
        layout.addWidget(header, nrow, 0, 1, 4)
        nrow += 1

        # Set the possible parameters
        output_params = {'outp_lat':     ['Latitude', 1, 0],
                         'outp_lon':     ['Longitude', 2, 0],
                         'outp_alt':     ['Altitude', 3, 0],
                         'outp_intlo':   ['Low Fit\nIntensity', 1, 1],
                         'outp_inthi':   ['High Fit\nIntensity', 2, 1],
                         'outp_intav':   ['Average Fit\nIntensity', 3, 1],
                         'outp_resmax':  ['Maximum\nResidual', 1, 2],
                         'outp_resstd':  ['Stdev.\nResidual', 2, 2],
                         'outp_fitqual': ['Fit Quality', 3, 2]}

        # Genertae the widgets
        for key, [label, row, col] in output_params.items():
            self.widgets[key] = QCheckBox(label)
            layout.addWidget(self.widgets[key], row, col)

# =============================================================================
# Tool Windows
# =============================================================================

    def open_ils_window(self):
        """Open ILS analysis."""
        win = ILSWindow(self)
        win.show()

    def open_flat_window(self):
        """Open flat analysis."""
        win = FLATWindow(self)
        win.show()

    def open_flux_window(self):
        """Open flux analysis."""
        win = CalcFlux(self)
        win.show()

    def open_ldf_window(self):
        """Open LDF analysis."""
        # Pull the plotting data from the GUI
        widgetData = {'gas_params':    self.gas_table.getData(),
                      'bgpoly_params': self.bgpoly_table.getData(),
                      'offset_params': self.offset_table.getData(),
                      'shift_params':  self.shift_table.getData()}

        for label in self.widgets:
            widgetData[label] = self.widgets.get(label)

        win = LDFWindow(widgetData, self)
        win.show()

    def update_plot_params(self):
        """Update plot parameter options."""
        rows = self.gas_table.getData()
        params = [r[0] for r in rows]
        self.widgets['graph_param'].clear()
        self.widgets['graph_param'].addItems(params)

    def closeEvent(self, event):
        """Handle GUI closure."""
        # Pull widget values
        config = {'gas_params':    self.gas_table.getData(),
                  'bgpoly_params': self.bgpoly_table.getData(),
                  'offset_params': self.offset_table.getData(),
                  'shift_params':  self.shift_table.getData()}
        for label in self.widgets:
            config[label] = self.widgets.get(label)
        config['theme'] = self.theme

        # Check if the config matches the current widget states
        save_flag = True
        for k in config.keys():
            if k not in self.config.keys() or config[k] != self.config[k]:
                save_flag = False

        # If there have been no changes, ask if want to quit
        if save_flag:
            reply = QMessageBox.question(self, 'Message',
                                         "Are you sure to quit?",
                                         QMessageBox.Yes | QMessageBox.No,
                                         QMessageBox.No)

            if reply == QMessageBox.Yes:
                event.accept()
            else:
                event.ignore()

        # If there have been changes, ask if want to save
        else:
            msg = "Would you like to save before quitting?"
            options = QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
            reply = QMessageBox.question(self, 'Message', msg, options,
                                         QMessageBox.Cancel)

            if reply == QMessageBox.Yes:
                self.save_config(asksavepath=False)
                event.accept()
            if reply == QMessageBox.No:
                event.accept()
            else:
                event.ignore()

# =============================================================================
# Save config
# =============================================================================

    def save_config(self, asksavepath=True):
        """Save the config file."""
        # Get the GUI configuration
        config = {'gas_params':    self.gas_table.getData(),
                  'bgpoly_params': self.bgpoly_table.getData(),
                  'offset_params': self.offset_table.getData(),
                  'shift_params':  self.shift_table.getData()}
        for label in self.widgets:
            config[label] = self.widgets.get(label)
        config['theme'] = self.theme

        # Get save filename if required
        if asksavepath or self.config_fname is None:
            filter = 'YAML (*.yml *.yaml);;All Files (*)'
            fname, s = QFileDialog.getSaveFileName(self, 'Save Config', '',
                                                   filter)
            # If valid, proceed. If not, return
            if fname != '' and fname is not None:
                self.config_fname = fname
            else:
                return

        # Write the config
        with open(self.config_fname, 'w') as outfile:
            yaml.dump(config, outfile)

        # Log the update
        logger.info(f'Config file saved to {self.config_fname}')

        #
        with open('bin/.config', 'w') as w:
            w.write(self.config_fname)

        self.config = config

# =============================================================================
# Load config
# =============================================================================

    def load_config(self, fname=None):
        """Read the config file."""
        if fname is None:
            filter = 'YAML (*.yml *.yaml);;All Files (*)'
            fname, tfile = QFileDialog.getOpenFileName(self, 'Load Config', '',
                                                       filter)

        # Open the config file
        try:
            with open(fname, 'r') as ymlfile:
                config = yaml.load(ymlfile, Loader=yaml.FullLoader)

            for label in config:
                try:
                    if label == 'gas_params':
                        self.gas_table.setRowCount(0)
                        self.gas_table.setData(config['gas_params'])

                    elif label == 'bgpoly_params':
                        self.bgpoly_table.setRowCount(0)
                        self.bgpoly_table.setData(config['bgpoly_params'])

                    elif label == 'offset_params':
                        self.offset_table.setRowCount(0)
                        self.offset_table.setData(config['offset_params'])

                    elif label == 'shift_params':
                        self.shift_table.setRowCount(0)
                        self.shift_table.setData(config['shift_params'])

                    elif label == 'theme':
                        self.theme = config[label]

                    else:
                        self.widgets.set(label, config[label])
                except Exception:
                    logger.warning(f'Failed to load {label} from config file')

            logger.info(f'Config file loaded from {self.config_fname}')

            # Update the config file settings
            self.config_fname = fname
            with open('bin/.config', 'w') as w:
                w.write(self.config_fname)

        except FileNotFoundError:
            logger.warning(f'Unable to load config file {self.config_fname}')
            config = {}
        self.config = config
        return config

# =============================================================================
#   Global Slots
# =============================================================================

    def update_progress(self, prog):
        """Slot to update the progress bar."""
        self.progress.setValue(prog)

    def update_status(self, status):
        """Update the status."""
        self.statusBar().showMessage(status)

    def update_error(self, error):
        """Slot to update error messages from the worker."""
        exctype, value, trace = error
        logger.warning(f'Uncaught exception!\n{trace}')

    def pause(self):
        """Pause the worker loop."""
        # Update button text
        if self.rt_pause_btn.text() == 'Pause':
            self.rt_pause_btn.setText('Continue')
        else:
            self.rt_pause_btn.setText('Pause')
        if self.pause_btn.text() == 'Pause':
            self.pause_btn.setText('Continue')
        else:
            self.pause_btn.setText('Pause')

        try:
            self.analysisWorker.pause()
        except (AttributeError, RuntimeError):
            pass
        try:
            self.specWorker.pause()
        except (AttributeError, RuntimeError):
            pass

    def stop(self):
        """Kill the worker loop."""
        try:
            self.analysisWorker.stop()
            logger.info('Analysis stopped')
        except (AttributeError, RuntimeError):
            pass
        try:
            self.specWorker.stop()
            logger.info('Acquisition stopped')
        except (AttributeError, RuntimeError):
            pass

# =============================================================================
#   Analysis Loop Setup
# =============================================================================

    def analysis_complete(self):
        """Slot to run once the analysis worker is finished."""
        # Renable the start button
        self.start_btn.setEnabled(True)
        self.pause_btn.setEnabled(False)
        self.stop_btn.setEnabled(False)
        self.pause_btn.setText('Pause')

        # Set the status bar
        self.update_status('Ready')

    def get_plot_data(self, plot_info):
        """Catch plot info emitted by the analysis loop."""
        # Unpack the data
        self.fit_result, self.spectrum, self.df, spec_no = plot_info

        # Get the parameter to plot
        self.key = self.widgets.get('graph_param')

        # Update the numerical output
        amt = self.fit_result.params[self.key].fit_val
        err = self.fit_result.params[self.key].fit_err
        self.spectrum_number.setText(str(int(spec_no)))
        self.last_amt.setText(f'{amt:.03g}')
        self.last_err.setText(f'{err:.03g}')

        self.update_plots()
        self.update_results_table()
        self.update_map()

    def update_plots(self):
        """Update the plots."""
        # Plot the data
        t1 = self.widgets.get('graph_flag')
        t2 = self.analysisWorker.is_paused
        if t1 and not t2:

            # Get the time series data
            plotx = np.array(self.df['Number'].to_numpy(), dtype=float)
            ploty = np.array(self.df[self.key].to_numpy(), dtype=float)
            erry = np.array(self.df[self.key + '_err'].to_numpy(),
                            dtype=float)

            # Remove failed fits
            fail_idx = np.where(~np.isnan(ploty))
            plotx = plotx[fail_idx]
            ploty = ploty[fail_idx]
            erry = erry[fail_idx]

            # Remove bad points if desired
            if self.widgets.get('good_fit_flag'):
                fit_quality = self.df['fit_quality'].to_numpy()[fail_idx]
                qual_idx = np.where(fit_quality == 1)[0]
                plotx = plotx[qual_idx]
                ploty = ploty[qual_idx]
                erry = erry[qual_idx]

            # Scroll the graphs if desired
            if self.widgets.get('scroll_flag'):
                npts = self.widgets.get('scroll_amt')
                if len(plotx) > npts:
                    plotx = plotx[-npts:]
                    ploty = ploty[-npts:]
                    erry = erry[-npts:]

            # Plot the data
            plot_data = {'spectrum': [self.fit_result.grid,
                                      self.fit_result.spec],
                         'fit': [self.fit_result.grid,
                                 self.fit_result.fit],
                         'measurement': self.spectrum,
                         'residual': [self.fit_result.grid,
                                      self.fit_result.resid],
                         'meas_od': [self.fit_result.grid,
                                     self.fit_result.meas_od[self.key]],
                         'fit_od': [self.fit_result.grid,
                                    self.fit_result.synth_od[self.key]],
                         'series': [plotx, ploty]}
            for key, data in plot_data.items():
                self.plot_lines[key].setData(*data)

            # Show error plot if selected
            if self.widgets.get('graph_err_flag'):
                self.err_plot.setVisible(True)
                self.err_plot.setData(x=plotx, y=ploty, top=erry, bottom=erry)
            else:
                self.err_plot.setVisible(False)

    def update_map(self):
        """Update data map."""
        if self.widgets.get('graph_flag'):
            # Get non-null values from the dataframe
            try:
                df = self.df.dropna(subset=['Lat', 'Lon'])
            except AttributeError:
                return

            # Remove bad points if desired
            if self.widgets.get('good_fit_flag'):
                df = df[df['fit_quality'] == 1]

            # Scroll the graphs if desired
            if self.widgets.get('scroll_flag'):
                npts = self.widgets.get('scroll_amt')
                df = df.tail(npts)

            # Pull the plot values
            lat = df['Lat'].to_numpy()
            lon = df['Lon'].to_numpy()
            values = df[self.key].to_numpy()

            # Get the colormap limits and normalise the data
            map_lo_lim, map_hi_lim = self.cbar.levels()
            norm_values = (values - map_lo_lim) / (map_hi_lim - map_lo_lim)

            # Convert to colors
            pens = [pg.mkPen(color=self.cmap.map(val)) for val in norm_values]
            brushes = [pg.mkBrush(color=self.cmap.map(val))
                       for val in norm_values]

            # Update the map data
            self.gps_scatter.setData(x=lon, y=lat, pen=pens, brush=brushes)

    def toggle_map_auto_scale(self):
        """Toggle manual control over map colorbar."""
        if self.widgets.get('auto_map_scale'):
            self.cbar.sigLevelsChangeFinished.connect(self.update_map)
        else:
            self.cbar.sigLevelsChangeFinished.connect(None)

    def initialize_results_table(self, params):
        """Initialize table rows."""
        # Clear all current rows
        self.results_table.clearContents()

        # Make the rows
        self.results_table.setRowCount(len(params.keys()))

        for i, [name, param] in enumerate(params.items()):
            self.results_table.setItem(i, 0, QTableWidgetItem(name))
            self.results_table.setItem(i, 1, QTableWidgetItem(str(param.vary)))

    def update_results_table(self):
        """Update results table with fit output."""
        # Pull the fitted parameters
        params = self.fit_result.params

        # Update the table items
        for i, [name, p] in enumerate(params.items()):
            self.results_table.setItem(i, 2, QTableWidgetItem(str(p.fit_val)))
            self.results_table.setItem(i, 3, QTableWidgetItem(str(p.fit_err)))

    def begin_analysis(self, analysis_mode):
        """Run spectral analysis."""
        # Pull the plotting data from the GUI
        widgetData = {'gas_params':    self.gas_table.getData(),
                      'bgpoly_params': self.bgpoly_table.getData(),
                      'offset_params': self.offset_table.getData(),
                      'shift_params':  self.shift_table.getData()}
        for label in self.widgets:
            widgetData[label] = self.widgets.get(label)

        # Read the dark spectra or assign from the GUI
        if widgetData['dark_flag'] and analysis_mode == 'rt_analyse':
            dark_spec = self.dark_spectrum
        elif widgetData['dark_flag'] and analysis_mode == 'post_analyse':
            dark_fnames = widgetData['dark_fnames'].split('\n')
            if len(dark_fnames) == 0:
                logger.warning('No dark spectra selected, disabling dark '
                               + 'correction')

            else:
                x, dark_spec = average_spectra(dark_fnames,
                                               widgetData['spec_type'],
                                               widgetData['wl_calib'])

        else:
            dark_spec = None

        # Initialise the analysis worker
        self.analysisThread = QThread()
        self.analysisWorker = AnalysisWorker(analysis_mode, widgetData,
                                             dark_spec)
        self.analysisWorker.moveToThread(self.analysisThread)
        self.analysisThread.started.connect(self.analysisWorker.run)
        self.analysisWorker.progress.connect(self.update_progress)
        self.analysisWorker.error.connect(self.update_error)
        self.analysisWorker.initializeTable.connect(
            self.initialize_results_table)
        self.analysisWorker.plotData.connect(self.get_plot_data)
        self.analysisWorker.finished.connect(self.analysis_complete)
        self.analysisWorker.finished.connect(self.analysisThread.quit)
        self.analysisWorker.finished.connect(self.analysisWorker.deleteLater)
        self.analysisThread.finished.connect(self.analysisThread.deleteLater)
        self.analysisThread.start()

        # Disable the start button and enable the pause/stop buttons
        self.start_btn.setEnabled(False)
        self.pause_btn.setEnabled(True)
        self.stop_btn.setEnabled(True)

# =============================================================================
# Connect to spectrometer
# =============================================================================

    def connect_spectrometer(self):
        """Connect or dissconnect the spectrometer."""
        if not self.spectro_connected_flag:

            # Connect to the spectrometer
            w = self.widgets
            spec = Spectrometer(integration_time=w.get("int_time"),
                                coadds=w.get("coadds"),
                                correct_dark_counts=w.get("nonlin_flag"),
                                correct_nonlinearity=w.get("eldark_flag")
                                )

            # Check if connection was successful
            if spec.serial_number is not None:

                # Add the spectrometer to the parent GUI
                self.spectrometer = spec

                # Update the GUI
                self.spec_id.setText(self.spectrometer.serial_number)
                self.spec_connect_btn.setText('Disconnect')

                # Create a holder for the dark spectra
                self.dark_spectrum = np.zeros(self.spectrometer.pixels)

                # Update GUI features
                self.spectro_connected_flag = True
                self.rt_flag_btn.setEnabled(True)
                self.acquire_darks_btn.setEnabled(True)
                self.update_inttime_btn.setEnabled(True)
                self.update_coadds_btn.setEnabled(True)
                self.rt_start_btn.setEnabled(True)
                for k in ["nonlin_flag", "eldark_flag"]:
                    self.widgets[k].setEnabled(False)
                    self.widgets[k].setStyleSheet("color: darkGray")

                # Begin scope acquisition
                self.start_scope_acquisition()

        else:

            # Kill scope acquisition
            self.scopeWorker.stop()
            self.scopeThread.quit()
            self.scopeThread.wait()

            # Disconnect the spectrometer
            self.spectrometer.close()

            # Update the GUI
            self.spec_id.setText('Not connected')
            self.spec_connect_btn.setText('Connect')

            # Update GUI features
            self.spectro_connected_flag = False
            self.rt_flag_btn.setEnabled(False)
            self.acquire_darks_btn.setEnabled(False)
            self.update_inttime_btn.setEnabled(False)
            self.update_coadds_btn.setEnabled(False)
            self.rt_start_btn.setEnabled(False)
            for k in ["nonlin_flag", "eldark_flag"]:
                self.widgets[k].setEnabled(True)
                self.widgets[k].setStyleSheet("color: white")

    def update_int_time(self):
        """Update the spectrometer integration time."""
        self.spectrometer.update_integration_time(self.widgets.get('int_time'))

    def update_coadds(self):
        """Update the spectrometer coadds."""
        self.spectrometer.update_coadds(self.widgets.get('coadds'))

    def toggle_fitting(self):
        """Toggle real time fitting on and off."""
        if self.rt_fitting_flag:
            self.rt_fitting_flag = False
            self.rt_flag_btn.setStyleSheet("background-color: red")
            self.rt_flag_btn.setText('Fitting OFF')
            logger.info('Fitting turned off')

        else:
            self.rt_fitting_flag = True
            self.rt_flag_btn.setStyleSheet("background-color: green")
            self.rt_flag_btn.setText('Fitting ON')
            logger.info('Fitting turned on')

# =============================================================================
#   Connect to GPS
# =============================================================================

    def connect_gps(self):
        """Input information for a GPS connection."""
        if not self.gps_connected_flag:
            dialog = GPSWizard(self)
            if dialog.exec_():
                self.gps = GPS(**dialog.gps_kwargs)
                self.gps_status.setText('Connected')
                self.gps_connect_btn.setText('Disconnect')
                self.gps_connected_flag = True
                logger.info('GPS connected')

                # Start GPS aquisition
                self.gpsThread = QThread()
                self.gpsWorker = GPSWorker(self.gps)
                self.gpsWorker.moveToThread(self.gpsThread)
                self.gpsThread.started.connect(self.gpsWorker.run)
                self.gpsWorker.position.connect(self.update_position)
                self.gpsWorker.error.connect(self.update_error)
                self.gpsWorker.finished.connect(self.gpsThread.quit)
                self.gpsWorker.finished.connect(self.gpsWorker.deleteLater)
                self.gpsThread.finished.connect(self.gpsThread.deleteLater)
                self.gpsThread.start()

        else:
            self.gpsWorker.stop()
            self.gpsThread.quit()
            self.gpsThread.wait()
            self.gps.close()
            self.gps = None
            self.gps_status.setText('Not connected')
            self.gps_connect_btn.setText('Connect')
            self.gps_connected_flag = False
            logger.info('GPS disconnected')

    def update_position(self, location):
        """Update the current position on the GUI."""
        # Unpack location data
        ts, lat, lon, alt = location

        # Determine the directions
        if lat >= 0:
            lat_dir = 'N'
        else:
            lat_dir = 'S'
        if lon >= 0:
            lon_dir = 'E'
        else:
            lon_dir = 'W'

        # Update GUI labels
        self.gps_timestamp.setText(f'{ts}')
        self.gps_lat.setText(f'{abs(lat):.6f} {lat_dir}')
        self.gps_lon.setText(f'{abs(lon):.6f} {lon_dir}')
        self.gps_alt.setText(f'{alt:.1f} m')

        # Update GUI graph
        xdata, ydata = self.gps_line.getData()

        if xdata is None or ydata is None:
            self.gps_line.setData([lon], [lat])

        elif lat != ydata[-1] and lon != xdata[-1]:
            new_xdata = np.append(xdata, lon)
            new_ydata = np.append(ydata, lat)
            self.gps_line.setData(new_xdata, new_ydata)

# =============================================================================
#   Acquisition Loop Setup
# =============================================================================

    def acquisition_complete(self):
        """Slot to run once the acquisition worker is finished."""
        # Renable the start button
        self.rt_start_btn.setEnabled(True)
        self.rt_pause_btn.setEnabled(False)
        self.rt_stop_btn.setEnabled(False)
        self.spec_connect_btn.setEnabled(True)
        self.acquire_darks_btn.setEnabled(True)
        self.rt_flag_btn.setEnabled(True)
        self.rt_pause_btn.setText('Pause')

        # Reset the range on the progress bar
        self.progress.setRange(0, 100)

        # Turn off the logger
        logger.removeHandler(self.acquisition_logger)

        # Set the status bar
        self.statusBar().showMessage('Ready')

        # Begin background measurements
        self.start_scope_acquisition()

    def plot_spectrum(self, spectrum):
        """Display measured spectrum."""
        self.scope_line.setData(*spectrum)

    def set_meas_spectrum(self, spectrum):
        """Update the spectrum to fit."""
        try:
            self.analysisWorker.set_spectrum(spectrum)
        except AttributeError:
            pass

    def set_dark_spectrum(self, dark_spec):
        """Update program dark spectrum."""
        self.dark_spectrum = dark_spec

    def start_scope_acquisition(self):
        """Set up scope acquisition."""
        # This section is for testing with a virtual spectrometer
        #######################################################################
        # self.spectrometer.fpath = 'Example/spectrum_00000.txt'
        #######################################################################

        # Initialise the acquisition worker
        self.scopeThread = QThread()
        self.scopeWorker = AcqScopeWorker(self.spectrometer)
        self.scopeWorker.moveToThread(self.scopeThread)
        self.scopeThread.started.connect(self.scopeWorker.run)
        self.scopeWorker.plotSpec.connect(self.plot_spectrum)
        self.scopeWorker.error.connect(self.update_error)
        self.scopeWorker.finished.connect(self.scopeThread.quit)
        self.scopeWorker.finished.connect(self.scopeWorker.deleteLater)
        self.scopeThread.finished.connect(self.scopeThread.deleteLater)
        self.scopeThread.start()

    def begin_acquisition(self, mode):
        """Set up spectra acquisition."""
        # Check a results folder has been chosen
        if self.widgets.get("rt_save_path") == '':
            logger.error('Please select an output folder!')
            return

        # Ensure the results folder exists
        if not os.path.exists(self.widgets.get("rt_save_path")):
            os.makedirs(self.widgets.get("rt_save_path"))

        # Create a log handler
        date_str = datetime.strftime(datetime.now(), "%Y-%m-%d")
        log_fname = f'{self.widgets.get("rt_save_path")}/{date_str}_iFit.log'
        self.acquisition_logger = logging.FileHandler(log_fname,
                                                      mode='a')
        log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        date_fmt = '%Y-%m-%d %H:%M:%S'
        f_format = logging.Formatter(log_fmt, date_fmt)
        self.acquisition_logger.setFormatter(f_format)
        logger.addHandler(self.acquisition_logger)

        # This section is for testing with a virtual spectrometer
        #######################################################################
        # if acquisition_mode == 'acquire_darks':
        #     self.spectrometer.fpath = 'Example/dark.txt'
        # if acquisition_mode == 'acquire_cont':
        #     self.spectrometer.fpath = 'Example/spectrum_00366.txt'
        #######################################################################

        # Pull the plotting data from the GUI
        widgetData = {'gas_params':    self.gas_table.getData(),
                      'bgpoly_params': self.bgpoly_table.getData(),
                      'offset_params': self.offset_table.getData(),
                      'shift_params':  self.shift_table.getData()}
        for label in self.widgets:
            widgetData[label] = self.widgets.get(label)

        # Stop the background acquisition
        self.scopeWorker.stop()
        self.scopeThread.quit()
        self.scopeThread.wait()

        # Initialise the acquisition thread and worker
        self.specThread = QThread()
        self.specWorker = AcqSpecWorker(self.spectrometer, self.gps,
                                        widgetData)
        self.specWorker.moveToThread(self.specThread)

        # Assign signals depending on mode
        if mode == 'acquire_darks':
            self.specThread.started.connect(self.specWorker.acquire_dark)
            self.progress.setRange(0, 100)
            self.specWorker.setDark.connect(self.set_dark_spectrum)

        else:
            self.specThread.started.connect(self.specWorker.acquire_spec)
            self.progress.setRange(0, 0)
            self.specWorker.setSpec.connect(self.set_meas_spectrum)

        # Assign signals
        self.specWorker.plotSpec.connect(self.plot_spectrum)
        self.specWorker.progress.connect(self.update_progress)
        self.specWorker.error.connect(self.update_error)
        self.specWorker.finished.connect(self.acquisition_complete)
        self.specWorker.finished.connect(self.specThread.quit)
        self.specWorker.finished.connect(self.specWorker.deleteLater)
        self.specThread.finished.connect(self.specThread.deleteLater)
        self.specThread.start()

        # Disable the start/acquisition buttons and enable the pause/stop
        # buttons
        self.rt_start_btn.setEnabled(False)
        self.rt_pause_btn.setEnabled(True)
        self.rt_stop_btn.setEnabled(True)
        self.spec_connect_btn.setEnabled(False)
        self.acquire_darks_btn.setEnabled(False)
        self.rt_flag_btn.setEnabled(False)

        # If running real time, launch the analyser loop
        if mode == 'acquire_cont' and self.rt_fitting_flag:
            self.begin_analysis('rt_analyse')

# =============================================================================
#   Gui Theme
# =============================================================================

    def change_theme(self):
        """Change the theme."""
        if self.theme == 'Light':
            self.changeThemeDark()
            self.theme = 'Dark'
        elif self.theme == 'Dark':
            self.changeThemeLight()
            self.theme = 'Light'

    @pyqtSlot()
    def changeThemeDark(self):
        """Change theme to dark."""
        darkpalette = QPalette()
        darkpalette.setColor(QPalette.Window, QColor(53, 53, 53))
        darkpalette.setColor(QPalette.WindowText, Qt.white)
        darkpalette.setColor(QPalette.Base, QColor(25, 25, 25))
        darkpalette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
        darkpalette.setColor(QPalette.ToolTipBase, Qt.black)
        darkpalette.setColor(QPalette.ToolTipText, Qt.white)
        darkpalette.setColor(QPalette.Text, Qt.white)
        darkpalette.setColor(QPalette.Button, QColor(53, 53, 53))
        darkpalette.setColor(QPalette.Active, QPalette.Button,
                             QColor(53, 53, 53))
        darkpalette.setColor(QPalette.ButtonText, Qt.white)
        darkpalette.setColor(QPalette.BrightText, Qt.red)
        darkpalette.setColor(QPalette.Link, QColor(42, 130, 218))
        darkpalette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        darkpalette.setColor(QPalette.HighlightedText, Qt.black)
        darkpalette.setColor(QPalette.Disabled, QPalette.ButtonText,
                             Qt.darkGray)
        QApplication.instance().setPalette(darkpalette)

        # Update graphs
        self.graphwin.setBackground('k')
        self.scopewin.setBackground('k')
        self.mapwin.setBackground('k')
        pen = pg.mkPen('w', width=1)

        axes = self.plot_axes + [self.scope_ax, self.map_ax]

        for ax in axes:
            ax.getAxis('left').setPen(pen)
            ax.getAxis('right').setPen(pen)
            ax.getAxis('top').setPen(pen)
            ax.getAxis('bottom').setPen(pen)
            ax.getAxis('left').setTextPen(pen)
            ax.getAxis('bottom').setTextPen(pen)
        # self.scope_ax.getAxis('left').setPen(pen)
        # self.scope_ax.getAxis('right').setPen(pen)
        # self.scope_ax.getAxis('top').setPen(pen)
        # self.scope_ax.getAxis('bottom').setPen(pen)

    @pyqtSlot()
    def changeThemeLight(self):
        """Change theme to light."""
        QApplication.instance().setPalette(self.style().standardPalette())
        self.graphwin.setBackground('w')
        self.scopewin.setBackground('w')
        self.mapwin.setBackground('w')
        pen = pg.mkPen('k', width=1)

        axes = self.plot_axes + [self.scope_ax, self.map_ax]

        for ax in axes:
            ax.getAxis('left').setPen(pen)
            ax.getAxis('right').setPen(pen)
            ax.getAxis('top').setPen(pen)
            ax.getAxis('bottom').setPen(pen)
            ax.getAxis('left').setTextPen(pen)
            ax.getAxis('bottom').setTextPen(pen)
        # self.scope_ax.getAxis('left').setPen(pen)
        # self.scope_ax.getAxis('right').setPen(pen)
        # self.scope_ax.getAxis('top').setPen(pen)
        # self.scope_ax.getAxis('bottom').setPen(pen)


# Cliet Code
def main():
    """Run main function."""
    # Create an instance of QApplication
    app = QApplication(sys.argv)

    app.setStyle("Fusion")

    # Show the GUI
    view = MainWindow()
    view.show()

    # Execute the main loop
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
