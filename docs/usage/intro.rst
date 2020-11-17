iFit
#####

Introduction
============

iFit is a spectral analysis tool designed to retrieve volcanic |SO2| slant column densities (SCDs) from passive measurements of scattered UV sunlight. It is a similar method to the traditional Differential Optical Absorption Spectroscopy (DOAS) approach commonly used by volcanologists, but has two main differences:

1) The forward model and all calculations therein are computed in intensity on a high resolution (0.01 nm spaced) wavelength grid, instead of in optical depth. This avoids the need for several corrections required by traditional DOAS, namely the I0-effect and saturation effects.

2) A high resolution literature solar spectrum is used in place of the measured reference spectrum, removing the chance for contamination of the reference spectrum with |SO2|.

More details on iFit can be found in following paper: `Esse et al. (2020) <https://doi.org/10.1016/j.jvolgeores.2020.107000>`_

Installation
============

To install iFit simply clone this repository. iFit is written in Python and requires version 3.6 or higher, with `numpy <https://numpy.org/doc/stable/>`_ and `scipy <http://docs.scipy.org/doc/scipy/reference>`_. An example script is given in ``iFit.py`` (note this also uses `matplotlib <https://matplotlib.org>`_, `pandas <http://pandas.pydata.org/pandas-docs/dev>`_ and `tqdm <https://tqdm.github.io/>`_).

There is also a Graphical User Interface (GUI), accessible by running ``iFitUI.py``. This additionally requires `PyQt5 <https://doc.qt.io/qtforpython/>`_, `pyqtgraph <http://www.pyqtgraph.org/>`_, `pyyaml <https://github.com/yaml/pyyaml>`_ and `pandas <http://pandas.pydata.org/pandas-docs/dev>`_.

The GUI can also handle spectra acquisition using the `python-seabreeze <https://github.com/ap--/python-seabreeze>`_ library with `Ocean Insight <https://www.oceaninsight.com/home>`_ (previously Ocean Optics) spectrometers.

The easiest way to install is using Anaconda. First create a new environment:

.. code-block:: bash

  conda create -n myenv python=3.6 numpy scipy matplotlib pyqt pyqtgraph pyyaml pandas

Setting myenv to whatever name you want to use for the environment. Then activate the environment once it is created:

.. code-block:: bash

  conda activate myenv

Finally, install python-seabreeze if required

.. code-block:: bash

  conda install -c conda-forge python-seabreeze

Alternatively an executable is available. Currently this is only for Windows (64-bit), however if you require an alternative version please get in touch and I will try to help.

.. Substitutions
.. |SO2| replace:: SO\ :sub:`2`
