# iFit

iFit is a program to retrieve SO2 column densities from scattered UV sunlight spectra, with the focus of measuring the SO2 flux from volcanoes. It also allows control of Ocean Optics USB series spectrometers, with simultaneous fitting.

## Install
The easiest way to use iFit is through the latest release, however this is currently only available for Windows 64-bit systems. If you require iFit for another operating system either use the source code (after installing python and the required libraries) or contact me and I will try to help!

iFit uses the Python Seabreeze library to talk to the spectormeters. Details can be found here: https://github.com/ap--/python-seabreeze

Whether using the executable or the source code you will need the drivers for the spectrometers. These drivers and how to install them can be found in the documentation for Python Seabreeze.

## calc_flux.py
A second program is included with iFit called calc_flux. This allows the user to calculate the flux of So2 from a volcano for a traverse measurement. The user inputs the SO2 time series from iFit, a GPS track and the wind speed to generate the flux.
