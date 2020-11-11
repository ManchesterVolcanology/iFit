iFit
#####

What is iFit?
==============

iFit is a spectral analysis tool designed to retrieve volcanic |SO2| slant column densities (SCDs) from passive measurements of scattered UV sunlight. It is a similar method to the traditional Differential Optical Absorption Spectroscopy (DOAS) approach commonly used by volcanologists, but has two main differences:

1) The forward model and all calculations therein are computed in intensity on a high resolution (0.01 nm spaced) wavelength grid, instead of in optical depth. This avoids the need for several corrections required by traditional DOAS, namely the I0-effect and saturation effects.

2) A high resolution literature solar spectrum is used in place of the measured reference spectrum, removing the chance for contamination of the reference spectrum with |SO2|.

More details on iFit can be found in following paper: `Esse et al. (2020) <https://doi.org/10.1016/j.jvolgeores.2020.107000>`_

Installation
=============

To install iFit simply clone this repository. iFit is written in Python and requires version 3.6 or higher, with numpy, scipy and matplotlib. An example script if given in `iFit.py`.

There is also a Graphical User Interface (GUI) written using PyQT. This requires PyQT5, pyqtgraph, pyyaml and pandas.

Alternatively an executable is available. Currently this is only for Windows (64-bit), however if you require an alternative version please get in touch and I will try to help.

.. Substitutions
.. |SO2| replace:: SO\ :sub:`2`
