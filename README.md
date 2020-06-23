# iFit

iFit is a program designed to retrieve SO2 column densities from scattered UV sunlight spectra, with the focus of measuring the SO2 flux from volcanoes. It also allows control of Ocean Optics USB series spectrometers, with simultaneous fitting.

## Install

iFit is written to work with basic python libraries contained in scientific distributions (e.g. Anaconda). Simply clone iFit and run with Python! Work is underway on producing an executable for windows (as well other operating systems if required).

## Using iFit
An example script is given in iFit.py, but more details are given here. iFit retrieves SO2 slant column densities by fitting measured scattered sunlight UV spectra that contain absorption features after passing through a volcanic plume. The basic principal works on the Beer-Lambert law, using the following equation:

<img src="https://www.codecogs.com/eqnedit.php?latex=I(\lambda)&space;=&space;G(x)&space;\otimes&space;\left(&space;I_0(\lambda)&space;\cdot&space;P'(\lambda)&space;\cdot&space;\exp&space;\left(&space;\Sigma_i&space;\left[&space;-\sigma_i(\lambda)&space;\cdot&space;a_i&space;\right&space;]&space;\right&space;)\right&space;)&space;&plus;&space;I_{offset}(\lambda)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?I(\lambda)&space;=&space;G(x)&space;\otimes&space;\left(&space;I_0(\lambda)&space;\cdot&space;P'(\lambda)&space;\cdot&space;\exp&space;\left(&space;\Sigma_i&space;\left[&space;-\sigma_i(\lambda)&space;\cdot&space;a_i&space;\right&space;]&space;\right&space;)\right&space;)&space;&plus;&space;I_{offset}(\lambda)" title="I(\lambda) = G(x) \otimes \left( I_0(\lambda) \cdot P'(\lambda) \cdot \exp \left( \Sigma_i \left[ -\sigma_i(\lambda) \cdot a_i \right ] \right )\right ) + I_{offset}(\lambda)" />

## Code function
iFit works on two main python objects: `Parameters` and `Analyser`.

### Parameters

A `Parameters` object contains all the information on the fit parameters used to analyse a spectrum. Each `Parameter` in the `Parameters` object has the following values:
- `name` identifies the `Parameter` (and so must be unique)
- `value` gives the initial guess for that `Parameter` in the fit
- `vary` controls whether that `Parameter` is allowed to be varied by the model

Additionally for gas cross-sections the `xpath` value is used to set the file path to the cross-section file.

So a `Parameters` object could be generated like this:

```python
from ifit.parameters import Parameters

# Create parameter dictionary
params = Parameters()

# Add the gases
params.add('SO2',  value=1.0e16, vary=True, xpath='Ref/SO2.txt')
params.add('O3',   value=1.0e19, vary=True, xpath='Ref/O3.txt')
params.add('Ring', value=0.1,    vary=True, xpath='Ref/Ring.txt')
    
# Add background polynomial parameters
params.add('bg_poly0', value=0.0, vary=True)
params.add('bg_poly1', value=0.0, vary=True)
params.add('bg_poly2', value=0.0, vary=True)
params.add('bg_poly3', value=1.0, vary=True)

# Add intensity offset parameters
params.add('offset0', value=0.0, vary=True)

# Add wavelength shift parameters
params.add('shift0', value=0.0, vary=True)
params.add('shift1', value=0.1, vary=True)
```

This defines three gas `Parameter`s for SO2, O3 and Ring, as well as the polynomial coefficients for the background polynomial, intensity offset and wavelength shift. Once the `Parameters` is defined the `Analyser` can be generated.

### Analyser
The `Analyser` handles the actual analysis of the spectra. It must be generated first, defining certain settings for the analysis, as well as the `Parameters` already defined:

```python

from ifit.spectral_analysis import Analyser

# Generate the analyser
analyser = Analyser(params       = params,
                    fit_window   = [310, 320],
                    frs_path     = 'Ref/sao2010.txt',
                    stray_flag   = True,
                    stray_window = [280, 290])
```

This will generate an analyser that will fit the emasured spectra between 310 - 320 nm, performing a stray light correction using the measured intensities between 280 - 290 nm.

Measured spectra can then be analysed by using `analyser.fit_spectrum`:

```python
fit = analyser.fit_spectrum([x,y])
```

In this case x and y are the measured spectum wavelengths and intensities respectively. This returns a `FitResult` object which holds the fit data and useful information.

The `FitResult` object contains a copy of the `Parameters` object that was passed to the `Analyser` with the `fit_val` and `fit_err` values populated with the optimised value and associated error for each `Parameter`. It also contains:
- `grid` the wavelength grid of the fit window
- `spec` the measured spectrum (after pre-processing) in the fit window
- `fit` the optimised model spectrum
- `resid` the residual between the measurement and the model

## In Progress
This is version 3.0 and is still in progress. Currently only post analysis is coded, but real time analysis with Ocean Optics spectrometer control is being added.

This update required a significant rewrite of the iFit software but should make it more user friendly and flexible for future applications. 
