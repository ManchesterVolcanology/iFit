# iFit

iFit is a program designed to retrieve SO2 column densities from scattered UV sunlight spectra, with the focus of measuring the SO2 flux from volcanoes. It also allows control of Ocean Optics USB series spectrometers, with simultaneous fitting.

## Install

### Python Script

iFit requires Python 3.6+ with `numpy` and `scipy` for basic operation. `matplotlib` and `pandas` are also useful for plotting and handelling data, but not strictly required. To use the GUI also requires `PySide6`, `pyqtgraph`, `pyqtdarktheme`, `pyyaml`, `utm` and `pandas`, as well as [`python-seabreeze`](https://github.com/ap--/python-seabreeze) for acquiring spectra with Ocean Insight spectrometers, and `pyserial` for communicating with a USB GPS unit.

Most libraries are available in Anaconda, which offers easy package management.  But first, to get the iFit code, clone the repository to your computer from the command line using

```
git clone https://github.com/ManchesterVolcanology/iFit.git
```

If you are working on Windows and do not have git, checkout git for windows.

### Anaconda Install

It is suggested to install Python and the various dependancies using Anaconda or Miniconda from https://anaconda.org/.

Once installed, create the environment for the iFit GUI from the terminal:

```
conda create -n ifit -c conda-forge numpy scipy pandas pyyaml pyserial utm pyside6 pyqtgraph seabreeze qt_theme::pyqtdarktheme
```

Note here we are using `ifit` as the environment name, but this can be changed to whatever you like. The current version of iFit is tested with certain versions of Python libraries, so you can also use the provided `environemnt.yml` to avoid any potential compatability issues with future versions of Python or the libraries used.

### Python Seabreeze

To install the relevant drivers for the python-seabreeze library to control Ocean Optics spectrometers run the following:

```
seabreeze_os_setup
```

Note that you may have to run this in a new terminal with root/administrator privilages. See the [`python-seabreeze docs`](https://python-seabreeze.readthedocs.io/en/latest/) for more details. Massive thanks to Andreas Poehlmann for creating and maintaining python-seabreeze, without which I would have spent an inordinate amount of time trying (and likely failing) to talk to spectrometers.

### Running the iFIt GUI

iFit should now be ready to run. To open the GUI, run the following within the iFit directory:

```
python iFitUI.py
```

### Executable

Executables have been compiled for Windows (32- and 64-bit versions) and are available on the [`releases`](https://github.com/benjaminesse/iFit/releases) page. These are currently in beta for further testing, please let me know if you have any issues with these. Currently there are no executables available for Mac or Linux, if you would like these please get in touch.

### Building the executable

The iFit executable is built with pyinstaller. If you want to build it yourself use the following command:

``pyinstaller -F -n iFitUI_win64 -w -i bin\icons\main.ico iFitUI.py``

This builds the executable into a single file, tells the program not to open a command window when launching and adds the icon file.
