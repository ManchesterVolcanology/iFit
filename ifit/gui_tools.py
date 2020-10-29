import sys
import numpy as np
import pyqtgraph as pg
from functools import partial
from scipy.optimize import curve_fit
from PyQt5.QtGui import QIcon, QPalette, QColor
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QMainWindow, QWidget, QApplication, QGridLayout,
                             QLabel, QTextEdit, QLineEdit, QPushButton,
                             QFileDialog, QScrollArea)

try:
    from .make_ils import super_gaussian
except ImportError:
    from make_ils import super_gaussian


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
    view = FLATWindow()
    view.show()

    # Execute the main loop
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
