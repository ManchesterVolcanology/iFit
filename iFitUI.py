import os
import sys
import logging
import numpy as np
import pyqtgraph as pg
from functools import partial
from PyQt5.QtGui import QIcon, QPalette, QColor
from PyQt5.QtCore import Qt, QThreadPool
from PyQt5.QtWidgets import (QMainWindow, QWidget, QApplication, QGridLayout,
                             QMessageBox, QLabel, QComboBox, QTextEdit,
                             QLineEdit, QPushButton, QProgressBar, QFrame,
                             QSplitter, QCheckBox, QSizePolicy, QSpacerItem,
                             QTabWidget, QAction)

from ifit.gui_functions import (control_loop, Widgets, SpinBox, Table, browse,
                                load_config, save_config, Worker,
                                QTextEditLogger)

__version__ = '3.3'
__author__ = 'Ben Esse'


class MainWindow(QMainWindow):
    """View for the iFit GUI"""

    def __init__(self):
        """View initialiser"""
        super().__init__()

        # Set the window properties
        self.setWindowTitle(f'iFit {__version__}')
        self.statusBar().showMessage('Ready')
        self.setGeometry(40, 40, 400, 140)
        self.setWindowIcon(QIcon('bin/icon.ico'))

        # Set the window layout
        self.generalLayout = QGridLayout()
        self._centralWidget = QWidget(self)
        self.setCentralWidget(self._centralWidget)
        self._centralWidget.setLayout(self.generalLayout)

        # Generate the threadpool
        self.threadpool = QThreadPool()

        # Setup widget stylesheets
        QTabWidget().setStyleSheet('QTabWidget { font-size: 18pt; }')

        # Create an empty dictionary to hold the GUI widgets
        self.widgets = Widgets()

        # Build the GUI
        self._createApp()

        # Update widgets from loaded config file
        load_config(self, fname='bin/config.yaml')

    def _createApp(self):
        """Handles building the main GUI"""

        # Add file menubar
        saveAct = QAction('&Save', self)
        saveAct.setShortcut('Ctrl+S')
        saveAct.triggered.connect(partial(save_config, self, False))
        saveasAct = QAction('&Save As', self)
        saveasAct.setShortcut('Ctrl+Shift+S')
        saveasAct.triggered.connect(partial(save_config, self, True))
        loadAct = QAction('&Load', self)
        loadAct.triggered.connect(partial(load_config, None))

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(saveAct)
        fileMenu.addAction(saveasAct)
        fileMenu.addAction(loadAct)

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
        """Builds the main GUI controls"""
        layout = QGridLayout(self.controlFrame)
        layout.setAlignment(Qt.AlignTop)

        # Create an option menu for the spectra format
        layout.addWidget(QLabel('Format:'), 0, 0)
        self.widgets.add('spec_type', QComboBox())
        self.widgets['spec_type'].addItems(['iFit',
                                            'Master.Scope',
                                            'Spectrasuite',
                                            'OpenSO2',
                                            'FLAME',
                                            'Basic'])
        self.widgets['spec_type'].setFixedSize(100, 20)
        layout.addWidget(self.widgets['spec_type'], 0, 1)

        # Add an input for the spectra selection
        layout.addWidget(QLabel('Spectra:'), 1, 0)
        self.widgets['spec_fnames'] = QTextEdit()
        layout.addWidget(self.widgets['spec_fnames'], 1, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.clicked.connect(partial(browse, self.widgets['spec_fnames'],
                                    'multi'))
        layout.addWidget(btn, 1, 4)

        # Add an input for the dark selection
        layout.addWidget(QLabel('Darks:'), 2, 0)
        self.widgets['dark_fnames'] = QTextEdit()
        layout.addWidget(self.widgets['dark_fnames'], 2, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.clicked.connect(partial(browse, self.widgets['dark_fnames'],
                                    'multi'))
        layout.addWidget(btn, 2, 4)

        # Add an input for the save selection
        layout.addWidget(QLabel('Save:'), 3, 0)
        self.widgets['save_path'] = QLineEdit()
        layout.addWidget(self.widgets['save_path'], 3, 1, 1, 3)
        btn = QPushButton('Browse')
        btn.clicked.connect(partial(browse, self.widgets['save_path'],
                                    'save'))
        layout.addWidget(btn, 3, 4)

        # Add button to begin analysis
        self.start_btn = QPushButton('Begin!')
        self.start_btn.clicked.connect(partial(self.begin))
        self.start_btn.setFixedSize(90, 25)
        layout.addWidget(self.start_btn, 4, 1)

        # Add button to pause analysis
        self.pause_btn = QPushButton('Pause')
        self.pause_btn.clicked.connect(partial(self.pause))
        self.pause_btn.setFixedSize(90, 25)
        layout.addWidget(self.pause_btn, 4, 2)

        # Add button to stop analysis
        self.stop_btn = QPushButton('Stop')
        self.stop_btn.clicked.connect(partial(self.stop))
        self.stop_btn.setFixedSize(90, 25)
        layout.addWidget(self.stop_btn, 4, 3)

# =============================================================================
#   Generate the program outputs
# =============================================================================

    def _createOuput(self):
        """Builds the main GUI visual ouputs"""
        layout = QGridLayout(self.outputFrame)
        layout.setAlignment(Qt.AlignTop)

        # Add a progress bar
        self.progress = QProgressBar(self)
        self.progress.setFixedSize(400, 25)
        layout.addWidget(self.progress, 0, 0, 1, 6)

        # Add numerical outputs
        layout.addWidget(QLabel('Last amt:'), 1, 0)
        self.last_amt = QLabel('-')
        layout.addWidget(self.last_amt, 1, 1)
        layout.addWidget(QLabel('+/-'), 1, 2)
        self.last_err = QLabel('-')
        layout.addWidget(self.last_err, 1, 3)

        self.logBox = QTextEditLogger(self)

        # log to text box
        self.logBox.setFormatter(logging.Formatter('%(message)s'))
        logging.getLogger().addHandler(self.logBox)
        logging.getLogger().setLevel(logging.INFO)

        # self.logBox.widget.setFixedSize(400, 150)
        layout.addWidget(self.logBox.widget, 2, 0, 1, 6)

        # log to file
        if not os.path.isdir('bin/'):
            os.makedirs('bin/')
        fh = logging.FileHandler('bin/iFit.log')
        fh.setLevel(logging.INFO)
        fmt = '%(asctime)s %(levelname)s %(module)s %(funcName)s %(message)s'
        fh.setFormatter(logging.Formatter(fmt))
        logging.getLogger().addHandler(fh)

        msg = 'Welcome to iFit! Written by Ben Esse'
        self.logBox.widget.appendPlainText(msg)

# =============================================================================
#   Set up graphs and settings
# =============================================================================

    def _createGraphs(self):
        """Build the graphical display and program settings"""
        layout = QGridLayout(self.graphFrame)

        # Generate tabs for the gaphs and settings
        tab1 = QWidget()
        tab2 = QWidget()

        # Form the tab widget
        tabwidget = QTabWidget()
        tabwidget.addTab(tab1, 'Graphs')
        tabwidget.addTab(tab2, 'Settings')
        tabwidget.setStyleSheet('QTabWidget { font-size: 18pt; }')
        layout.addWidget(tabwidget, 0, 0)

# =============================================================================
#       Set up the graphs
# =============================================================================

        graphwin = pg.GraphicsLayoutWidget(show=True)
        pg.setConfigOptions(antialias=True)
        graphwin.resize(1000, 600)

        glayout = QGridLayout(tab1)

        # Make the graphs
        ax0 = graphwin.addPlot(row=0, col=0)
        ax1 = graphwin.addPlot(row=0, col=1)
        ax2 = graphwin.addPlot(row=1, col=0)
        ax3 = graphwin.addPlot(row=1, col=1)
        ax4 = graphwin.addPlot(row=2, col=0, colspan=2)

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
        p0 = pg.mkPen(color='#1f77b4', width=2)
        p1 = pg.mkPen(color='#ff7f0e', width=2)
        l0 = ax0.plot([], [], pen=p0, symbol='+', symbolSize=5,
                      symbolBrush=('#1f77b4'), name='Spectrum')
        l1 = ax0.plot([], [], pen=p1, name='Fit')
        l2 = ax1.plot([], [], pen=p0)
        l3 = ax2.plot([], [], pen=p0, symbol='+', symbolSize=5,
                      symbolBrush=('#1f77b4'))
        l4 = ax3.plot([], [], pen=p0, symbol='+', symbolSize=5,
                      symbolBrush=('#1f77b4'))
        l5 = ax3.plot([], [], pen=p1)
        l6 = ax4.plot([], [], pen=p0, symbol='+', symbolSize=5,
                      symbolBrush=('#1f77b4'))

        ax0.addLegend()

        self.plot_lines = [l0, l1, l2, l3, l4, l5, l6]
        self.plot_axes = [ax0, ax0, ax1, ax2, ax3, ax3, ax4]

        for ax in self.plot_axes:
            ax.showGrid(x=True, y=True)

        # Add the graphs to the layout
        glayout.addWidget(graphwin, 0, 0, 1, 7)

# =============================================================================
#      Graph settings
# =============================================================================

        # Create a checkbox to turn plotting on or off
        self.widgets['graph_flag'] = QCheckBox('Show Graphs?')
        glayout.addWidget(self.widgets['graph_flag'], 1, 0)

        # Add combo box for the graph parameter
        glayout.addWidget(QLabel('Parameter to graph:'), 1, 1)
        self.widgets['graph_param'] = QComboBox()
        self.widgets['graph_param'].addItems([''])
        # self.widgets['graph_param'].setFixedSize(70, 20)
        glayout.addWidget(self.widgets['graph_param'], 1, 2)

        # Create a checkbox to turn scrolling on or off
        self.widgets['scroll_flag'] = QCheckBox('Scroll Graphs?')
        glayout.addWidget(self.widgets['scroll_flag'], 1, 3)

        # Add spinbox for the graph scroll amount
        glayout.addWidget(QLabel('No. Spectra\nTo Display:'), 1, 4)
        self.widgets['scroll_amt'] = SpinBox(100, [1, 10000])
        # self.widgets['scroll_amt'].setFixedSize(70, 20)
        glayout.addWidget(self.widgets['scroll_amt'], 1, 5)

        vspacer = QSpacerItem(QSizePolicy.Minimum, QSizePolicy.Expanding)
        glayout.addItem(vspacer, 1, 6, 1, -1)

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
        nrow = 0
        ncol = 1

        # Add spinboxs for the fit window
        layout.addWidget(QLabel('Fit Window:\n    (nm)'), nrow, ncol, 2, 1)
        self.widgets['fit_lo'] = SpinBox(310, [0, 10000])
        self.widgets['fit_lo'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['fit_lo'], nrow, ncol+1)
        nrow += 1
        self.widgets['fit_hi'] = SpinBox(320, [0, 10000])
        self.widgets['fit_hi'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['fit_hi'], nrow, ncol+1)
        nrow += 1

        # Add spinbox for the model grid padding
        layout.addWidget(QLabel('Model Grid\nPadding (nm):'), nrow, ncol)
        self.widgets['model_padding'] = SpinBox(1.0, [0, 10000])
        self.widgets['model_padding'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['model_padding'], nrow, ncol+1)
        nrow += 1

        # Add spinbox for the model grid spacing
        layout.addWidget(QLabel('Model Grid\nSpacing (nm):'), nrow, ncol)
        self.widgets['model_spacing'] = SpinBox(0.01, [0, 10])
        self.widgets['model_spacing'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['model_spacing'], nrow, ncol+1)
        nrow += 1

        # Add combo box for interpolation method
        layout.addWidget(QLabel('Interpolation\nMethod:'), nrow, ncol)
        self.widgets['interp_method'] = QComboBox()
        self.widgets['interp_method'].addItems(['cubic', 'linear'])
        self.widgets['interp_method'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['interp_method'], nrow, ncol+1)

        # New column
        layout.addWidget(QVLine(), 0, ncol+2, 10, 1)
        nrow = 0
        ncol += 3

        # Add sterio button for dark correction
        self.widgets['dark_flag'] = QCheckBox('Correct Dark\nSpectrum?')
        layout.addWidget(self.widgets['dark_flag'], nrow, ncol, 1, 2)
        nrow += 1

        # Add sterio button for flat correction
        self.widgets['flat_flag'] = QCheckBox('Correct Flat\nSpectrum?')
        layout.addWidget(self.widgets['flat_flag'], nrow, ncol, 1, 2)
        nrow += 1

        # Add spinboxs for the stray light window
        layout.addWidget(QLabel('Stray Light\nWindow: (nm)'), nrow, ncol, 2, 1)
        self.widgets['stray_lo'] = SpinBox(280, [0, 10000])
        self.widgets['stray_lo'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['stray_lo'], nrow, ncol+1)
        nrow += 1
        self.widgets['stray_hi'] = SpinBox(290, [0, 10000])
        self.widgets['stray_hi'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['stray_hi'], nrow, ncol+1)
        nrow += 1

        # Add sterio button for stray light correction
        self.widgets['stray_flag'] = QCheckBox('Remove?')
        layout.addWidget(self.widgets['stray_flag'], nrow-2, ncol+2, 2, 1)

        # Add spinbox to control spike removal
        layout.addWidget(QLabel('Spike Limit\n(counts):'), nrow, ncol)
        self.widgets['spike_limit'] = SpinBox(1000, [0, 10000000])
        self.widgets['spike_limit'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['spike_limit'], nrow, ncol+1)
        self.widgets['despike_flag'] = QCheckBox('Remove?')
        layout.addWidget(self.widgets['despike_flag'], nrow, ncol+2)
        nrow += 1

        # New column
        layout.addWidget(QVLine(), 0, ncol+3, 10, 1)
        nrow = 0
        ncol += 4

        # Add combo box for residual display
        layout.addWidget(QLabel('Residual Display:'), nrow, ncol)
        self.widgets['resid_type'] = QComboBox()
        self.widgets['resid_type'].addItems(['Percentage', 'Absolute'])
        self.widgets['resid_type'].setFixedSize(100, 20)
        layout.addWidget(self.widgets['resid_type'], nrow, ncol+1)
        nrow += 1

        # Add sterio button for auto-update of fit params
        self.widgets['update_flag'] = QCheckBox('Auto-update\nFit Parameters?')
        layout.addWidget(self.widgets['update_flag'], nrow, ncol, 1, 2)
        nrow += 1

        # Add spinbox for the residual limit
        layout.addWidget(QLabel('Residual Limit:'), nrow, ncol)
        self.widgets['resid_limit'] = SpinBox(1.0, [0, 10000])
        self.widgets['resid_limit'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['resid_limit'], nrow, ncol+1)
        nrow += 1

        # Add spinboxs for the intensity limits
        layout.addWidget(QLabel('Intensity Limits:'), nrow, ncol, 2, 1)
        self.widgets['lo_int_limit'] = SpinBox(0, [0, 100000])
        self.widgets['lo_int_limit'].setFixedSize(70, 20)
        layout.addWidget(self.widgets['lo_int_limit'], nrow, ncol+1)
        nrow += 1
        self.widgets['hi_int_limit'] = SpinBox(70000, [0, 100000])
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
        self.widgets['ils_mode'].addItems(['Params', 'File', 'Manual'])
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
        btn.clicked.connect(partial(browse, self.widgets['ils_path'],
                                    'single'))
        layout.addWidget(btn, nrow, ncol+3)
        nrow += 1

        # Add an input for the flat spectrum
        layout.addWidget(QLabel('Flat Spectrum:'), nrow, ncol)
        self.widgets['flat_path'] = QLineEdit()
        self.widgets['flat_path'].setFixedSize(300, 25)
        layout.addWidget(self.widgets['flat_path'], nrow, ncol+1, 1, 2)
        btn = QPushButton('Browse')
        btn.setFixedSize(100, 25)
        btn.clicked.connect(partial(browse, self.widgets['flat_path'],
                                    'single'))
        layout.addWidget(btn, nrow, ncol+3)
        nrow += 1

        # Add an input for the wavelength calibration
        layout.addWidget(QLabel('Wavelength\nCalibration:'), nrow, ncol)
        self.widgets['wl_calib'] = QLineEdit()
        self.widgets['wl_calib'].setFixedSize(300, 25)
        layout.addWidget(self.widgets['wl_calib'], nrow, ncol+1, 1, 2)
        btn = QPushButton('Browse')
        btn.setFixedSize(100, 25)
        btn.clicked.connect(partial(browse, self.widgets['wl_calib'],
                                    'single'))
        layout.addWidget(btn, nrow, ncol+3)
        nrow += 1

        layout.addWidget(QHLine(), nrow, 0, 1, 10)
        nrow += 1

        # Add inputs for the manual ILS parameters
        layout.addWidget(QLabel('FWEM:'), nrow, ncol)
        self.widgets['fwem'] = QLineEdit()
        self.widgets['fwem'].setFixedSize(100, 20)
        layout.addWidget(self.widgets['fwem'], nrow, ncol+1)
        self.widgets['fwem_fit'] = QCheckBox('Fit?')
        layout.addWidget(self.widgets['fwem_fit'], nrow, ncol+2)
        nrow += 1

        layout.addWidget(QLabel('k:'), nrow, ncol)
        self.widgets['k'] = QLineEdit()
        self.widgets['k'].setFixedSize(100, 20)
        layout.addWidget(self.widgets['k'], nrow, ncol+1)
        self.widgets['k_fit'] = QCheckBox('Fit?')
        layout.addWidget(self.widgets['k_fit'], nrow, ncol+2)
        nrow += 1

        layout.addWidget(QLabel('a_w:'), nrow, ncol)
        self.widgets['a_w'] = QLineEdit()
        self.widgets['a_w'].setFixedSize(100, 20)
        layout.addWidget(self.widgets['a_w'], nrow, ncol+1)
        self.widgets['a_w_fit'] = QCheckBox('Fit?')
        layout.addWidget(self.widgets['a_w_fit'], nrow, ncol+2)
        nrow += 1

        layout.addWidget(QLabel('a_k:'), nrow, ncol)
        self.widgets['a_k'] = QLineEdit()
        self.widgets['a_k'].setFixedSize(100, 20)
        layout.addWidget(self.widgets['a_k'], nrow, ncol+1)
        self.widgets['a_k_fit'] = QCheckBox('Fit?')
        layout.addWidget(self.widgets['a_k_fit'], nrow, ncol+2)
        nrow += 1

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
        btn.clicked.connect(partial(browse, self.widgets['frs_path'],
                                    'single'))
        layout.addWidget(btn, nrow, ncol+2)
        nrow += 1

        # Create the tabs for the parameters
        ptab1 = QWidget()
        ptab2 = QWidget()
        ptab3 = QWidget()
        ptab4 = QWidget()

        tabwidget = QTabWidget()
        tabwidget.setFixedWidth(600)
        tabwidget.addTab(ptab1, 'Absorbers')
        tabwidget.addTab(ptab2, 'Polynomial')
        tabwidget.addTab(ptab3, 'Offset')
        tabwidget.addTab(ptab4, 'Shift')
        layout.addWidget(tabwidget, 1, 0, 1, 6)

        vspacer = QSpacerItem(QSizePolicy.Minimum, QSizePolicy.Expanding)
        layout.addItem(vspacer, nrow, 0, 1, -1)

        hspacer = QSpacerItem(QSizePolicy.Expanding, QSizePolicy.Minimum)
        layout.addItem(hspacer, 0, 5, -1, 1)

        # Create the absorber and polynomial tables
        self.gas_table = Table(ptab1, 'param', 550)
        self.bgpoly_table = Table(ptab2, 'poly', 250, 'bg_poly')
        self.offset_table = Table(ptab3, 'poly', 250, 'offset')
        self.shift_table = Table(ptab4, 'poly', 250, 'shift')

        # Link the parameter table to the plot parameter combobox
        self.gas_table.cellChanged.connect(self.update_plot_params)

    def update_plot_params(self):
        """Updates plot parameter options"""
        rows = self.gas_table.getData()
        params = [r[0] for r in rows]
        self.widgets['graph_param'].clear()
        self.widgets['graph_param'].addItems(params)

    def closeEvent(self, event):
        """Handle GUI closure"""
        reply = QMessageBox.question(self, 'Message',
                                     "Are you sure to quit?",
                                     QMessageBox.Yes | QMessageBox.No,
                                     QMessageBox.No)

        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

# =============================================================================
#   Analysis Loop Setup
# =============================================================================

    def thread_complete(self):
        """Slot to run once the worker is finished"""
        # Renable the start button
        self.start_btn.setEnabled(True)

    def update_progress(self, prog):
        """Slot to update the progress bar"""
        self.progress.setValue(prog)

    def update_plots(self, plot_info):
        """Slot to update the plots"""
        # Unpack the data
        fit_result, spectrum, df = plot_info
        key = self.widgets.get('graph_param')

        # Update the numerical output
        amt = fit_result.params[key].fit_val
        err = fit_result.params[key].fit_err
        self.last_amt.setText(f'{amt:.03g}')
        self.last_err.setText(f'{err:.03g}')

        if self.widgets.get('graph_flag'):

            plotx = df['Number'].dropna().to_numpy()
            ploty = df[key].dropna().to_numpy()

            # Check the length of the time series parameters
            if self.widgets.get('scroll_flag'):
                nscroll = self.widgets.get('scroll_amt')
                npts = len(ploty)
                if npts > nscroll:
                    plotx = plotx[int(npts-nscroll):]
                    ploty = ploty[int(npts-nscroll):]

            # Pack the data to plot
            data = [[fit_result.grid, fit_result.spec],
                    [fit_result.grid, fit_result.fit],
                    spectrum,
                    [fit_result.grid, fit_result.resid],
                    [fit_result.grid, fit_result.meas_od[key]],
                    [fit_result.grid, fit_result.synth_od[key]],
                    [plotx, ploty]]

            # And plot!
            for i, l in enumerate(self.plot_lines):
                x, y = data[i]
                if np.nanmax(y) > 1e6:
                    order = int(np.ceil(np.log10(np.nanmax(y)))) - 1
                    y = y / 10**order
                    self.plot_axes[i].setLabel('left',
                                               f'Fit value (1e{order})')

                l.setData(x, y)

    def begin(self):
        """Function to set up and start the analysis worker"""
        self.worker = Worker(control_loop, gui=self)
        self.worker.signals.finished.connect(self.thread_complete)
        self.worker.signals.progress.connect(self.update_progress)
        self.worker.signals.plotter.connect(self.update_plots)
        self.threadpool.start(self.worker)

    def pause(self):
        """Pauses the worker loop"""
        self.worker.pause()

    def stop(self):
        """Kills the worker loop"""
        self.worker.kill()
        logging.info('Analysis aborted')


class QHLine(QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)


class QVLine(QFrame):
    def __init__(self):
        super(QVLine, self).__init__()
        self.setFrameShape(QFrame.VLine)
        self.setFrameShadow(QFrame.Sunken)


# Cliet Code
def main():
    """Main function"""
    # Create an instance of QApplication
    app = QApplication(sys.argv)

    app.setStyle("Fusion")

    # Now use a palette to switch to dark colors:
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(53, 53, 53))
    palette.setColor(QPalette.WindowText, Qt.white)
    palette.setColor(QPalette.Base, QColor(25, 25, 25))
    palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    palette.setColor(QPalette.ToolTipBase, Qt.black)
    palette.setColor(QPalette.ToolTipText, Qt.white)
    palette.setColor(QPalette.Text, Qt.white)
    palette.setColor(QPalette.Button, QColor(53, 53, 53))
    palette.setColor(QPalette.ButtonText, Qt.white)
    palette.setColor(QPalette.BrightText, Qt.red)
    palette.setColor(QPalette.Link, QColor(42, 130, 218))
    palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    palette.setColor(QPalette.HighlightedText, Qt.black)
    app.setPalette(palette)

    # Show the GUI
    view = MainWindow()
    view.show()

    # Execute the main loop
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
