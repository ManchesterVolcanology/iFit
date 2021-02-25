# iFit

iFit is a program designed to retrieve SO<sub>2</sub> column densities from scattered UV sunlight spectra, with the focus of measuring the SO<sub>2</sub> flux from volcanoes. It also allows control of Ocean Optics USB series spectrometers, with simultaneous fitting.

iFit's documentation can be found here: https://ifit.readthedocs.io/en/latest/

## Install

### Executable
Executables have been compiled for Windows (32- and 64-bit versions) and are available on the [`releases`](https://github.com/benjaminesse/iFit/releases) page. These are currently in beta for further testing, please let me know if you have any issues with these. Currently there are no executables available for Mac or Linux, if you would like these please get in touch.

### Python Script
iFit requires Python 3.6+ with `numpy` and `scipy` for basic operation. `matplotlib` and `pandas` are also useful for plotting and handelling data, but not strictly required. To use the GUI also requires `PyQT5`, `pyqtgraph`, `pyyaml` and `pandas`, as well as [`python-seabreeze`](https://github.com/ap--/python-seabreeze) for acquiring spectra with Ocean Insight spectrometers.

To get started with Anaconda create a new environment:

```conda create -n myenv python=3.6 numpy scipy matplotlib tqdm pandas pyqt pyqtgraph pyyaml```

where `myenv` is the name of the environment. The environment can then be activated with:

```
conda activate myenv
```

This should allow the use of the basic code as well as the GUI. If wanting to acquire spectra then also run:

```
conda install -c conda-forge seabreeze
seabreeze_os_setup
```

To install `python-seabreeze` and install the relevant drivers. See the [`python-seabreeze docs`](https://python-seabreeze.readthedocs.io/en/latest/) for more details. Massive thanks to Andreas Poehlmann for creating and maintaining python-seabreeze, without which I would have spent an inordinate amount of time trying (and likely failing) to talk to spectrometers.

## Icons

Some icons by [Yusuke Kamiyamane](https://p.yusukekamiyamane.com/). Licensed under a [Creative Commons Attribution 3.0 License](https://creativecommons.org/licenses/by/3.0/).
