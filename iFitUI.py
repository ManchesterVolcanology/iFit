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
from PyQt5.QtCore import Qt, QThreadPool
from PyQt5.QtWidgets import (QMainWindow, QWidget, QApplication, QGridLayout,
                             QMessageBox, QLabel, QComboBox, QTextEdit,
                             QLineEdit, QPushButton, QProgressBar, QFrame,
                             QSplitter, QCheckBox, QSizePolicy, QSpacerItem,
                             QTabWidget, QAction, QFileDialog, QScrollArea,
                             QToolBar)

from ifit.gui_functions import (analysis_loop, acquire_spectra, Widgets,
                                SpinBox, DSpinBox, Table, Worker,
                                QTextEditLogger, connect_spectrometer)
from ifit.gui_tools import ILSWindow, FLATWindow, CalcFlux, LDFWindow

__version__ = '3.3'
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

        # Generate the threadpool for launching background processes
        self.threadpool = QThreadPool()

        # Setup widget stylesheets
        QTabWidget().setStyleSheet('QTabWidget { font-size: 18pt; }')

        # Create an empty dictionary to hold the GUI widgets
        self.widgets = Widgets()

        # Build the GUI
        self._createApp()

        # Update widgets from loaded config file
        self.config = {}
        self.config_fname = None
        if os.path.isfile('bin/.config'):
            with open('bin/.config', 'r') as r:
                self.config_fname = r.readline().strip()
            self.load_config(fname=self.config_fname)

    def _createApp(self):
        """Build the main GUI."""
        # Add file menubar
        saveAct = QAction(QIcon('bin/icons/save.png'), '&Save', self)
        saveAct.setShortcut('Ctrl+S')
        saveAct.triggered.connect(partial(self.save_config, False))
        saveasAct = QAction(QIcon('bin/icons/saveas.png'), '&Save As',
                            self)
        saveasAct.setShortcut('Ctrl+Shift+S')
        saveasAct.triggered.connect(partial(self.save_config, True))
        loadAct = QAction(QIcon('bin/icons/open.png'), '&Load', self)
        loadAct.triggered.connect(partial(self.load_config, None))

        # Add tools menubar
        ilsAct = QAction(QIcon('bin/icons/ils.png'), '&Measure ILS', self)
        ilsAct.triggered.connect(self.open_ils_window)
        flatAct = QAction(QIcon('bin/icons/flat.png'), '&Measure\nFlat Field',
                          self)
        flatAct.triggered.connect(self.open_flat_window)
        fluxAct = QAction(QIcon('bin/icons/flux.png'), '&Calculate flux', self)
        fluxAct.triggered.connect(self.open_flux_window)
        ldfAct = QAction(QIcon('bin/icons/ldf.png'),
                         '&Light Dilution\nAnalysis', self)
        ldfAct.triggered.connect(self.open_ldf_window)

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

        # Create a toolbar
        toolbar = QToolBar("Main toolbar")
        self.addToolBar(toolbar)
        toolbar.addAction(saveAct)
        toolbar.addAction(saveasAct)
        toolbar.addAction(loadAct)
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

        # Set a label for the spectrometer ID
        self.connected_flag = False
        layout.addWidget(QLabel('Spectrometer:'), 0, 0)
        self.spec_id = QLabel('Not connected')
        self.spec_id.setToolTip('Spectrometer Serial Number')
        layout.addWidget(self.spec_id, 0, 1)

        # Create a button to connect to a spectrometer
        self.connect_btn = QPushButton('Connect')
        self.connect_btn.setToolTip('Connect or disconnect the spectrometer')
        self.connect_btn.clicked.connect(partial(connect_spectrometer, self))
        layout.addWidget(self.connect_btn, 0, 2)

        # Create a button to acquire a test spectrum
        self.acquire_test_btn = QPushButton('Test Spectrum')
        self.acquire_test_btn.setToolTip('Acquire test spectrum displayed in '
                                         + 'Scope')
        self.acquire_test_btn.clicked.connect(partial(self.begin_acquisition,
                                              'acquire_single'))
        self.acquire_test_btn.setEnabled(False)
        layout.addWidget(self.acquire_test_btn, 0, 3)

        # Create a control for the spectrometer integration time
        layout.addWidget(QLabel('Integration\nTime (ms):'), 1, 0)
        self.widgets['int_time'] = SpinBox(100, [10, 1000000])
        self.widgets['int_time'].setToolTip('Spectrometer integration time')
        layout.addWidget(self.widgets['int_time'], 1, 1)

        # Create a button to update the integration time
        self.update_inttime_btn = QPushButton('Update')
        self.update_inttime_btn.setToolTip('Update Spectrometer Integration'
                                           + ' Time')
        self.update_inttime_btn.clicked.connect(self.update_int_time)
        self.update_inttime_btn.setEnabled(False)
        layout.addWidget(self.update_inttime_btn, 1, 2)

        # Create a button to toggle real-time analysis
        self.rt_fitting_flag = False
        self.rt_flag_btn = QPushButton('Fitting OFF')
        self.rt_flag_btn.setToolTip('Toggle Real Time Fitting')
        self.rt_flag_btn.clicked.connect(self.toggle_fitting)
        self.rt_flag_btn.setEnabled(False)
        self.rt_flag_btn.setStyleSheet("background-color: red")
        layout.addWidget(self.rt_flag_btn, 0, 4)

        # Create a control for the spectrometer coadds
        layout.addWidget(QLabel('Coadds:'), 2, 0)
        self.widgets['coadds'] = SpinBox(10, [1, 1000000])
        self.widgets['coadds'].setToolTip('No. spectra to average')
        layout.addWidget(self.widgets['coadds'], 2, 1)

        # Create a button to update the coadds
        self.update_coadds_btn = QPushButton('Update')
        self.update_coadds_btn.setToolTip('Update Spectrometer Coadds')
        self.update_coadds_btn.clicked.connect(self.update_coadds)
        self.update_coadds_btn.setEnabled(False)
        layout.addWidget(self.update_coadds_btn, 2, 2)

        # Create a control for the number of dark spectra
        layout.addWidget(QLabel('No. Dark\nSpectra:'), 3, 0)
        self.widgets['ndarks'] = SpinBox(10, [1, 1000000])
        self.widgets['ndarks'].setToolTip('Set number of dark spectra to '
                                          + 'measure')
        layout.addWidget(self.widgets['ndarks'], 3, 1)

        # Create a button to acquire the dark spectra
        self.acquire_darks_btn = QPushButton('Acquire')
        self.acquire_darks_btn.setToolTip('Measure Dark Spectra')
        self.acquire_darks_btn.clicked.connect(partial(self.begin_acquisition,
                                               'acquire_darks'))
        self.acquire_darks_btn.setEnabled(False)
        layout.addWidget(self.acquire_darks_btn, 3, 2)

        # Add stereo button for non-liniarity correction
        self.widgets['nonlin_flag'] = QCheckBox('Correct\nNon-Linearity?')
        self.widgets['nonlin_flag'].setToolTip('Turn on correction for non-'
                                               + 'linear intensity response'
                                               + ' correction')
        layout.addWidget(self.widgets['nonlin_flag'], 1, 3)

        # Add stereo button for non-liniarity correction
        self.widgets['eldark_flag'] = QCheckBox('Correct\nElectronic dark?')
        self.widgets['eldark_flag'].setToolTip('Turn on electronic dark '
                                               + 'correction')
        layout.addWidget(self.widgets['eldark_flag'], 2, 3)

        # Add an input for the save selection
        layout.addWidget(QLabel('Output\nFolder:'), 4, 0)
        self.widgets['rt_save_path'] = QLineEdit()
        self.widgets['rt_save_path'].setToolTip('Folder to hold real time '
                                                + 'results')
        layout.addWidget(self.widgets['rt_save_path'], 4, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.setFixedSize(70, 25)
        btn.clicked.connect(partial(browse, self, self.widgets['rt_save_path'],
                                    'folder', None))
        layout.addWidget(btn, 4, 4)

        # Add button to begin analysis
        self.rt_start_btn = QPushButton('Begin!')
        self.rt_start_btn.setToolTip('Begin Spectra Acquisition')
        self.rt_start_btn.clicked.connect(partial(self.begin_acquisition,
                                                  'acquire_cont'))
        self.rt_start_btn.setFixedSize(90, 25)
        self.rt_start_btn.setEnabled(False)
        layout.addWidget(self.rt_start_btn, 5, 1)

        # Add button to pause analysis
        self.rt_pause_btn = QPushButton('Pause')
        self.rt_pause_btn.setToolTip('Pause/Play Spectra Acquisition')
        self.rt_pause_btn.clicked.connect(partial(self.pause))
        self.rt_pause_btn.setFixedSize(90, 25)
        self.rt_pause_btn.setEnabled(False)
        layout.addWidget(self.rt_pause_btn, 5, 2)

        # Add button to stop analysis
        self.rt_stop_btn = QPushButton('Stop')
        self.rt_stop_btn.setToolTip('Stop Spectra Acquisition')
        self.rt_stop_btn.clicked.connect(partial(self.stop))
        self.rt_stop_btn.setFixedSize(90, 25)
        self.rt_stop_btn.setEnabled(False)
        layout.addWidget(self.rt_stop_btn, 5, 3)

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
        self.logBox.setFormatter(logging.Formatter('%(message)s'))
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

# =============================================================================
#       Set up the analysis graphs
# =============================================================================

        self.graphwin = pg.GraphicsLayoutWidget(show=True)
        pg.setConfigOptions(antialias=True)
        # pg.setConfigOptions(useOpenGL=True)

        glayout = QGridLayout(tab1)

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
        p0 = pg.mkPen(color='#1f77b4', width=1.0)
        p1 = pg.mkPen(color='#ff7f0e', width=1.0)
        l0 = ax0.plot(pen=p0, name='Spectrum')
        l1 = ax0.plot(pen=p1, name='Fit')
        l2 = ax1.plot(pen=p0)
        l3 = ax2.plot(pen=p0)
        l4 = ax3.plot(pen=p0)
        l5 = ax3.plot(pen=p1)
        l6 = ax4.plot(pen=p0)

        ax0.addLegend()

        self.plot_lines = [l0, l1, l2, l3, l4, l5, l6]

        # Add the graphs to the layout
        glayout.addWidget(self.graphwin, 0, 0, 1, 7)

# =============================================================================
#      Graph settings
# =============================================================================

        # Create a checkbox to turn plotting on or off
        self.widgets['graph_flag'] = QCheckBox('Show Graphs?')
        self.widgets['graph_flag'].setToolTip('Display graph plots (can slow '
                                              + 'analysis)')
        glayout.addWidget(self.widgets['graph_flag'], 1, 0)

        # Add combo box for the graph parameter
        glayout.addWidget(QLabel('Parameter to graph:'), 1, 1)
        self.widgets['graph_param'] = QComboBox()
        self.widgets['graph_param'].addItems([''])
        # self.widgets['graph_param'].setFixedSize(70, 20)
        glayout.addWidget(self.widgets['graph_param'], 1, 2)

        # Create a checkbox to turn scrolling on or off
        self.widgets['scroll_flag'] = QCheckBox('Scroll Graphs?')
        self.widgets['scroll_flag'].setToolTip('Allow graphs to scroll\n'
                                               + '(limits no. spectra '
                                               + 'displayed)')
        glayout.addWidget(self.widgets['scroll_flag'], 1, 3)

        # Add spinbox for the graph scroll amount
        glayout.addWidget(QLabel('No. Spectra\nTo Display:'), 1, 4)
        self.widgets['scroll_amt'] = SpinBox(100, [1, 10000])
        # self.widgets['scroll_amt'].setFixedSize(70, 20)
        glayout.addWidget(self.widgets['scroll_amt'], 1, 5)

        # Create a checkbox to only display good fits
        self.widgets['good_fit_flag'] = QCheckBox('Only Show\nGood Fits?')
        self.widgets['good_fit_flag'].setToolTip('Only display results for '
                                                 + 'fits that pass the '
                                                 + 'quality checks')
        glayout.addWidget(self.widgets['good_fit_flag'], 2, 0)

        # Add combo box for the graphbackground color
        glayout.addWidget(QLabel('Graph Background:'), 2, 1)
        self.widgets['graph_bg'] = QComboBox()
        self.widgets['graph_bg'].setToolTip('Select light/dark graph '
                                            + 'background')
        self.widgets['graph_bg'].addItems(['Dark', 'Light'])
        self.widgets['graph_bg'].currentTextChanged.connect(self.alt_graph_bg)
        glayout.addWidget(self.widgets['graph_bg'], 2, 2)

        vspacer = QSpacerItem(QSizePolicy.Minimum, QSizePolicy.Expanding)
        glayout.addItem(vspacer, 1, 6, 1, -1)

# =============================================================================
#      Set up the scope plot
# =============================================================================

        self.scopewin = pg.GraphicsLayoutWidget(show=True)
        glayout = QGridLayout(tab3)

        # Make the graph
        self.scope_ax = self.scopewin.addPlot(row=0, col=0)
        self.scope_ax.setDownsampling(mode='peak')
        self.scope_ax.setClipToView(True)
        self.scope_ax.showGrid(x=True, y=True)

        self.scope_line = self.scope_ax.plot([], [], pen=p0)

        # Add the graphs to the layout
        glayout.addWidget(self.scopewin, 0, 0)

# =============================================================================
#       Create settings
# =============================================================================

        # Create tabs for settings
        slayout = QGridLayout(tab2)

        stab1 = QWidget()
        stab2 = QWidget()
        stab3 = QWidget()

        tabwidget = QTabWidget()
        tabwidget.addTab(stab1, 'Model')
        tabwidget.addTab(stab2, 'Spectrometer')
        tabwidget.addTab(stab3, 'Parameters')
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
        nrow += 1

        # Add sterio button for flat correction
        self.widgets['flat_flag'] = QCheckBox('Correct Flat\nSpectrum?')
        self.widgets['flat_flag'].setToolTip('Remove flat-field spectrum '
                                             + 'before fitting')
        layout.addWidget(self.widgets['flat_flag'], nrow, ncol, 1, 2)
        nrow += 1

        # Add spinboxs for the stray light window
        layout.addWidget(QLabel('Stray Light\nWindow: (nm)'), nrow, ncol, 2, 1)
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
        layout.addWidget(QLabel('ILS Parameters:'), nrow, ncol)
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
        self.gas_table = Table(ptab1, 'param', 620, 400)
        self.bgpoly_table = Table(ptab2, 'poly', 550, 400, 'bg_poly')
        self.offset_table = Table(ptab3, 'poly', 550, 400, 'offset')
        self.shift_table = Table(ptab4, 'poly', 550, 400, 'shift')

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

    def alt_graph_bg(self):
        """Change the graph background color."""
        var = self.widgets.get('graph_bg')

        if var == 'Light':
            bgcolor = 'w'
            fgcolor = 'k'
        if var == 'Dark':
            bgcolor = 'k'
            fgcolor = 'w'

        self.graphwin.setBackground(bgcolor)
        self.scopewin.setBackground(bgcolor)
        pen = pg.mkPen(fgcolor, width=1)

        for ax in self.plot_axes:
            ax.getAxis('left').setPen(pen)
            ax.getAxis('right').setPen(pen)
            ax.getAxis('top').setPen(pen)
            ax.getAxis('bottom').setPen(pen)
        self.scope_ax.getAxis('left').setPen(pen)
        self.scope_ax.getAxis('right').setPen(pen)
        self.scope_ax.getAxis('top').setPen(pen)
        self.scope_ax.getAxis('bottom').setPen(pen)

    def closeEvent(self, event):
        """Handle GUI closure."""
        # Pull widget values
        config = {'gas_params':    self.gas_table.getData(),
                  'bgpoly_params': self.bgpoly_table.getData(),
                  'offset_params': self.offset_table.getData(),
                  'shift_params':  self.shift_table.getData()}
        for label in self.widgets:
            config[label] = self.widgets.get(label)

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
        config = {'gas_params':    self.gas_table.getData(),
                  'bgpoly_params': self.bgpoly_table.getData(),
                  'offset_params': self.offset_table.getData(),
                  'shift_params':  self.shift_table.getData()}

        for label in self.widgets:
            config[label] = self.widgets.get(label)

        if asksavepath or self.config_fname is None:
            filter = 'YAML (*.yml *.yaml);;All Files (*)'
            fname, s = QFileDialog.getSaveFileName(self, 'Save Config', '',
                                                   filter)
            if fname != '':
                self.config_fname = fname

        with open(self.config_fname, 'w') as outfile:
            yaml.dump(config, outfile)

        logger.info(f'Config file saved to {self.config_fname}')

        with open('bin/.config', 'w') as w:
            w.write(self.config_fname)

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
        self.statusBar().showMessage('Ready')

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

        self.update_graph_flag = True

        self.update_plots()

    def update_plots(self):
        """Update the plots."""
        # See if the graph data has been updated
        if self.update_graph_flag:

            # Plot the data
            if self.widgets.get('graph_flag') and not self.worker.is_paused:

                # Get the time series data
                plotx = np.array(self.df['Number'].to_numpy(), dtype=float)
                ploty = np.array(self.df[self.key].to_numpy(), dtype=float)

                # Remove failed fits
                fail_idx = np.where(~np.isnan(ploty))
                plotx = plotx[fail_idx]
                ploty = ploty[fail_idx]

                # Remove bad points if desired
                if self.widgets.get('good_fit_flag'):
                    fit_quality = self.df['fit_quality'].to_numpy()[fail_idx]
                    qual_idx = np.where(fit_quality == 1)[0]
                    plotx = plotx[qual_idx]
                    ploty = ploty[qual_idx]

                # Scroll the graphs if desired
                if self.widgets.get('scroll_flag'):
                    npts = self.widgets.get('scroll_amt')
                    if len(plotx) > npts:
                        plotx = plotx[-npts:]
                        ploty = ploty[-npts:]

                # Plot the data
                self.plot_lines[0].setData(self.fit_result.grid,
                                           self.fit_result.spec)
                self.plot_lines[1].setData(self.fit_result.grid,
                                           self.fit_result.fit)
                self.plot_lines[2].setData(*self.spectrum)
                self.plot_lines[3].setData(self.fit_result.grid,
                                           self.fit_result.resid)
                self.plot_lines[4].setData(self.fit_result.grid,
                                           self.fit_result.meas_od[self.key])
                self.plot_lines[5].setData(self.fit_result.grid,
                                           self.fit_result.synth_od[self.key])
                self.plot_lines[6].setData(plotx, ploty)

                # Turn off autoscaling the x axis after the first plot
                if self.autoscale_flag:
                    self.autoscale_flag = False
                    for i in [0, 1, 2, 3]:
                        self.plot_axes[i].enableAutoRange('x', False)

                self.update_graph_flag = False

    def begin_analysis(self, analysis_mode):
        """Set up and start the analysis worker."""
        # Pull the plotting data from the GUI
        widgetData = {'gas_params':    self.gas_table.getData(),
                      'bgpoly_params': self.bgpoly_table.getData(),
                      'offset_params': self.offset_table.getData(),
                      'shift_params':  self.shift_table.getData()}

        for label in self.widgets:
            widgetData[label] = self.widgets.get(label)

        # Initialise the analysis worker
        self.worker = Worker(analysis_loop, analysis_mode, widgetData)
        self.worker.signals.finished.connect(self.analysis_complete)
        self.worker.signals.progress.connect(self.update_progress)
        self.worker.signals.status.connect(self.update_status)
        self.worker.signals.plotter.connect(self.get_plot_data)
        self.worker.signals.error.connect(self.update_error)
        self.threadpool.start(self.worker)

        # Disable the start button and enable the pause/stop buttons
        self.start_btn.setEnabled(False)
        self.pause_btn.setEnabled(True)
        self.stop_btn.setEnabled(True)

        # Set plot x limits where known
        self.autoscale_flag = True

        # Initialise the plotting timer
        self.update_graph_flag = False

# =============================================================================
#   Acquisition Loop Setup
# =============================================================================

    def acquisition_complete(self):
        """Slot to run once the acquisition worker is finished."""
        # Renable the start button
        self.rt_start_btn.setEnabled(True)
        self.rt_pause_btn.setEnabled(False)
        self.rt_stop_btn.setEnabled(False)
        self.connect_btn.setEnabled(True)
        self.acquire_test_btn.setEnabled(True)
        self.acquire_darks_btn.setEnabled(True)
        self.update_inttime_btn.setEnabled(True)
        self.update_coadds_btn.setEnabled(True)
        self.rt_flag_btn.setEnabled(True)
        self.rt_pause_btn.setText('Pause')

        # Reset the range on the progress bar
        self.progress.setRange(0, 100)

        # Turn off the logger
        logger.removeHandler(self.acquisition_logger)

        # Set the status bar
        self.statusBar().showMessage('Ready')

    def catch_spectrum(self, spec_data):
        """Slot to catch the spectra acquired by the acquisition loop."""
        self.last_spectrum, self.last_info, plot_flag = spec_data

        self.scope_line.setData(*self.last_spectrum)

        try:
            self.worker.set_spectrum(self.last_info['fname'])
        except AttributeError:
            pass

    def begin_acquisition(self, acquisition_mode):
        """Set up and start the acquisition worker."""
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
        # if acquisition_mode == 'acquire_single':
        #     self.spectrometer.fpath = 'Example/spectrum_00000.txt'
        # if acquisition_mode == 'acquire_darks':
        #     self.spectrometer.fpath = 'Example/dark.txt'
        # if acquisition_mode == 'acquire_cont':
        #     self.spectrometer.fpath = 'Example/spectrum_00366.txt'
        #######################################################################

        # Set the progress bar to busy
        if acquisition_mode != 'acquire_darks':
            self.progress.setRange(0, 0)

        # Pull the plotting data from the GUI
        widgetData = {'gas_params':    self.gas_table.getData(),
                      'bgpoly_params': self.bgpoly_table.getData(),
                      'offset_params': self.offset_table.getData(),
                      'shift_params':  self.shift_table.getData()}

        for label in self.widgets:
            widgetData[label] = self.widgets.get(label)

        # Initialise the acquisition worker
        self.acq_worker = Worker(acquire_spectra, acquisition_mode, widgetData,
                                 self.spectrometer)
        self.acq_worker.signals.finished.connect(self.acquisition_complete)
        self.acq_worker.signals.spectrum.connect(self.catch_spectrum)
        self.acq_worker.signals.progress.connect(self.update_progress)
        self.acq_worker.signals.status.connect(self.update_status)
        self.acq_worker.signals.error.connect(self.update_error)
        self.threadpool.start(self.acq_worker)

        # Disable the start/acquisition buttons and enable the pause/stop
        # buttons
        self.rt_start_btn.setEnabled(False)
        self.rt_pause_btn.setEnabled(True)
        self.rt_stop_btn.setEnabled(True)
        self.connect_btn.setEnabled(False)
        self.acquire_test_btn.setEnabled(False)
        self.acquire_darks_btn.setEnabled(False)
        self.update_inttime_btn.setEnabled(False)
        self.update_coadds_btn.setEnabled(False)
        self.rt_flag_btn.setEnabled(False)

        # If running real time, launch the analyser loop
        if acquisition_mode == 'acquire_cont' and self.rt_fitting_flag:
            self.begin_analysis('rt_analyse')

        # Set plot x limits where known
        self.autoscale_flag = True

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
            self.worker.pause()
        except AttributeError:
            pass
        try:
            self.acq_worker.pause()
        except AttributeError:
            pass

    def stop(self):
        """Kill the worker loop."""
        try:
            self.worker.kill()
            logger.info('Analysis stopped')
        except AttributeError:
            pass
        try:
            self.acq_worker.kill()
            logger.info('Acquisition stopped')
        except AttributeError:
            pass


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


class QHLine(QFrame):
    """Horizontal line widget."""

    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)


class QVLine(QFrame):
    """Vertical line widget."""

    def __init__(self):
        super(QVLine, self).__init__()
        self.setFrameShape(QFrame.VLine)
        self.setFrameShadow(QFrame.Sunken)


# Cliet Code
def main():
    """Run main function."""
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
    view = MainWindow()
    view.show()

    # Execute the main loop
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
