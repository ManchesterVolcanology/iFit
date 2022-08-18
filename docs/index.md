## Welcome to iFit Documentation Test

This will be populated with the documentation shortly...

## Introduction

iFit is a spectral analysis tool designed to retrieve volcanic |SO2| slant column densities (SCDs) from passive measurements of scattered UV sunlight. It is a similar method to the traditional Differential Optical Absorption Spectroscopy (DOAS) approach commonly used by volcanologists, but has two main differences:

1) The forward model and all calculations therein are computed in intensity on a high resolution (0.01 nm spaced) wavelength grid, instead of in optical depth. This avoids the need for several corrections required by traditional DOAS, namely the I0-effect and saturation effects.

2) A high resolution literature solar spectrum is used in place of the measured reference spectrum, removing the chance for contamination of the reference spectrum with |SO2|.

More details on iFit can be found in following paper: [Esse et al. (2020)](https://doi.org/10.1016/j.jvolgeores.2020.107000)

## Installation

To install iFit simply clone this repository. iFit is written in Python and requires version 3.6 or higher, with [`numpy`](https://numpy.org/doc/stable/) and [`scipy`](http://docs.scipy.org/doc/scipy/reference). An example script is given in `iFit.py` (note this also uses [`matplotlib](https://matplotlib.org), [`pandas'](http://pandas.pydata.org/pandas-docs/dev) and [`tqdm`](https://tqdm.github.io/)).

There is also a Graphical User Interface (GUI), accessible by running ``iFitUI.py``. This additionally requires [`PyQt5`](https://doc.qt.io/qtforpython/), [`pyqtgraph`](http://www.pyqtgraph.org/), [`pyyaml`](https://github.com/yaml/pyyaml) and [`pandas`](http://pandas.pydata.org/pandas-docs/dev).

The GUI can also handle spectra acquisition using the [`python-seabreeze`](https://github.com/ap--/python-seabreeze) library with [`Ocean Insight`](https://www.oceaninsight.com/home) (previously Ocean Optics) spectrometers.