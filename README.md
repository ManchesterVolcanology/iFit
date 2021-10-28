# iFit

iFit is a program designed to retrieve SO<sub>2</sub> column densities from scattered UV sunlight spectra, with the focus of measuring the SO<sub>2</sub> flux from volcanoes. It also allows control of Ocean Optics USB series spectrometers, with simultaneous fitting.

iFit's documentation can be found here: https://ifit.readthedocs.io/en/latest/

## Install

### Executable
Executables have been compiled for Windows (32- and 64-bit versions) and are available on the [`releases`](https://github.com/benjaminesse/iFit/releases) page. These are currently in beta for further testing, please let me know if you have any issues with these. Currently there are no executables available for Mac or Linux, if you would like these please get in touch.

### Building the executable
The iFit executable is build with pyinstaller. If you want to build it yourself use the following command:

```pyinstaller -F -n iFitUI_win64 -w -i bin\icons\main.ico iFitUI.py```

This builds the executable into a single file, tells the program not to open a command window when launching and adds the icon file.

### Python Script
iFit requires Python 3.6+ with `numpy` and `scipy` for basic operation. `matplotlib` and `pandas` are also useful for plotting and handelling data, but not strictly required. To use the GUI also requires `PyQT5`, `pyqtgraph`, `pyyaml` and `pandas`, as well as [`python-seabreeze`](https://github.com/ap--/python-seabreeze) for acquiring spectra with Ocean Insight spectrometers.

To get started with Anaconda create a new environment:

```conda create -n myenv python=3.6 numpy scipy matplotlib tqdm pandas pyqt pyqtgraph pyyaml```

where `myenv` is the name of the environment. The environment can then be activated with:

```
conda activate myenv
```

Note that currently (checked October 2021) the codnda version of pyqtgraph is outdated. If you have problems with pyqtgraph then uninstall from conda and add with pip:
```
conda uninstall pyqtgraph
pip install pyqtgraph
```

This should allow the use of the basic code as well as the GUI. If wanting to acquire spectra then also run:

```
conda install -c conda-forge seabreeze
```

To install `python-seabreeze`, the library for cotrolly Ocean Optics (or Ocean Insight) spectrometers. To install the relevant drivers shut the terminal, reopen a new one and run

```
seabreeze_os_setup
```

after activating the environment. See the [`python-seabreeze docs`](https://python-seabreeze.readthedocs.io/en/latest/) for more details. Massive thanks to Andreas Poehlmann for creating and maintaining python-seabreeze, without which I would have spent an inordinate amount of time trying (and likely failing) to talk to spectrometers.
