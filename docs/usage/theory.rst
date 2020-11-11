Theory
#######

iFit retrieves |SO2| slant column densities by fitting measured scattered sunlight UV spectra that contain absorption features after passing through a volcanic plume. The basic principal works on the Beer-Lambert law:

:math:`I(\lambda) = I_0(\lambda) \cdot \exp\left(\sum_{i}\left[-\sigma_i (\lambda) \cdot a_i \right] \right)`

Where :math:`I(\lambda)` is the measured intensity of light, :math:`I_0(\lambda)` is the initial intensity and :math:`\sigma_i (\lambda)` and :math:`a_i` are the absorption cross-section and column density of species :math:`i` respectively.

This is for an ideal case, but in reality the forward model must be odified to take into account additional signals and instrument effects. Firstly, two polynomial terms are included: :math:`P(\lambda)` to take account of broadband aerosol scattering and overal changes in intesity, and :math:`I_{offset}(\lambda)` to take into account any offset intensity due to stray light in the detector or uncorrected instrument offset.

The spectrometers typically used to measure volcanic |SO2| are not able to resolve the sharp absorption lines in either the solar spectrum or the |SO2| cross-section, so the forward model must be smoothed in order to compare with a measured spectrum. This is achieved by convoloving the modelled spectrum with the instrument line shape (ILS), :math:`G(x)`.

Adding these factors gives the iFit forward model:

:math:`I(\lambda) = G(x) \otimes \left( I_0^*(\lambda) \cdot P(\lambda) \cdot \exp \left(\sum_{i}\left[-\sigma_i (\lambda) \cdot a_i \right] \right)  \right) + I_{offset}(\lambda)`

Finally, thermal fluctuations can cause expansion or contraction of components within the spectrometer, which causes changes to the spectrometer wavelength calibration (however the wavelength of the emasured light changes across the detector pixels). This shift is taken into account by a shift and squeeze to the modelled wavelength grid.

To retrieve the |SO2| SCD the iFit forward model is fitted to each measured spectrum using a non-linear least-squares minimization in which the fitted parameters are varied in order to minimize the residual between the model and measurement.

.. Substitutions
.. |SO2| replace:: SO\ :sub:`2`
