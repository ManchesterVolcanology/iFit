# iFit

iFit is a program designed to retrieve SO2 column densities from scattered UV sunlight spectra, with the focus of measuring the SO2 flux from volcanoes. It also allows control of Ocean Optics USB series spectrometers, with simultaneous fitting.

iFit's documentation can be found here: https://ifit.readthedocs.io/en/latest/

## Install

### Python Script

iFit requires Python 3.6+ with `numpy` and `scipy` for basic operation. `matplotlib` and `pandas` are also useful for plotting and handelling data, but not strictly required. To use the GUI also requires `PySide6`, `pyqtgraph`, `pyyaml` and `pandas`, as well as [`python-seabreeze`](https://github.com/ap--/python-seabreeze) for acquiring spectra with Ocean Insight spectrometers.

Most libraries are available in Anaconda, which offers easy package management.  But first, to get the iFit code, clone the repository to your computer from the command line using

``git clone https://github.com/benjaminesse/iFit.git``

If you are working on Windows and do not have git, checkout git for windows.

### Anaconda Install

Install Anaconda or Miniconda from https://anaconda.org/.

Once installed, use conda to create the environment for iFit. Not here we will use iFit as the environment name, but feel free to name it whatever you like.

```
conda create -n iFit
conda activate iFit
```

This create and activates the environment. Now install the required libraries:

```
conda install numpy scipy tqdm pandas pyyaml pyserial
conda install -c conda-forge utm pyside2 pyqtgraph seabreeze
pip install pyqtdarktheme
```

Note that pyqtdarktheme is not available on conda so must be installed using pip.

### Standard Python install

iFit can also be installed using the standard Python distribution. Firstly download the latest Python release from https://www.python.org/.

From within the iFit directory (or another location you choose) create the virtual environment:

```
python -m venv venv
```

This will create a local directory holding the virtual environment. To activate it run:

```
venv\Scripts\activate
```

on windows or

```
source venv/bin/activate
```

on Unix or MacOS. To test the install type `python` and hit enter. This should open an interactive command prompt for python. To exit this type `exit()` and hit enter.

Next we need to install the libraries required for iFit. For basic functionality you only need:

```
pip install numpy scipy
```

For full functionality, including the GUI, install the following:

```
pip install numpy scipy matplotlib tqdm pandas pyyaml pyserial utm PySide2 pyqtgraph seabreeze pyqtdarktheme
```

### Python Seabreeze

To install the relevant drivers for the python-seabreeze library to control Ocean Insight (was Ocean Optics) spectrometers run the following:

```
seabreeze_os_setup
```

Note that you may have to run this in a new terminal with root/administrator privilages. See the [`python-seabreeze docs`](https://python-seabreeze.readthedocs.io/en/latest/) for more details. Massive thanks to Andreas Poehlmann for creating and maintaining python-seabreeze, without which I would have spent an inordinate amount of time trying (and likely failing) to talk to spectrometers.

### Running the iFIt GUI

iFit should now be ready to run. To open the GUI, run:

```
python iFitUI.py
```

### Executable

Executables have been compiled for Windows (32- and 64-bit versions) and are available on the [`releases`](https://github.com/benjaminesse/iFit/releases) page. These are currently in beta for further testing, please let me know if you have any issues with these. Currently there are no executables available for Mac or Linux, if you would like these please get in touch.

### Building the executable

The iFit executable is build with pyinstaller. If you want to build it yourself use the following command:

``pyinstaller -F -n iFitUI_win64 -w -i bin\icons\main.ico iFitUI.py``

This builds the executable into a single file, tells the program not to open a command window when launching and adds the icon file.
