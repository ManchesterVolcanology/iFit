"""Contains functions to calculate light dilution curves."""
import logging
import numpy as np

logger = logging.getLogger(__name__)


def generate_ld_curves(analyser, spectrum, wb1=[306, 316], wb2=[312, 322],
                       so2_lims=[0, 1e19], so2_step=5e17, ldf_lims=[0.0, 0.9],
                       ldf_step=0.1):
    """Generate light dilution curves.

    Parameters
    ----------
    analyser : ifit Analyser object
        Analyser used to fit the data. This must be initialised to cover both
        wavebands used and should not have LDF as a parameters.
    spectrum : 2D numpy array
        The spectrum used to make the light dilution curves as [x, y].
    wb1 : tuple, optional
        The wavelength limits of wave band 1 (the shorter band). Default is
        [306, 316].
    wb2: tuple, optional
        The wavelength limits of wave band 2 (the longer band). Default is
        [312, 322].
    so2_lims : tuple, optional
        The SO2 SCD limits to produce the light dilution curves over, measured
        in molecules/cm2. Default is [0, 1e19].
    so2_step : float, optional
        The step size for the SO2 SCD grid. Default is 5e17
    ldf_lims : tuple, optional
        The LDF grid limits to produce the curves. Default is [0.0, 0.9].
    ldf_step : float, optional
        The step size for the LDF grid. Default is 0.1.

    Returns
    -------
    something : float
        Desc.
    """
    # Build the SO2 and LDF grids
    so2_grid = np.arange(so2_lims[0], so2_lims[1]+so2_step, so2_step)
    ldf_grid = np.arange(ldf_lims[0], ldf_lims[1]+ldf_step, ldf_step)

    # Fit the initial spectrum across the whole window
    logger.info('Running initial fit...')
    init_fit = analyser.fit_spectrum(spectrum)

    # Update the parameters with the fit results
    analyser.params.update_values(init_fit.params.popt_list())

    # Add the LDF to the analyser params
    analyser.params.add('LDF', value=0.0, vary=False)

    # Initialise an empty array to hold the synthetic spectra
    shape = [len(so2_grid), len(init_fit.spec), len(ldf_grid)]
    spectra = np.zeros(shape)

    # Generate the syntheic spectra
    logger.info('Generating synthetic spectra...')
    for j, ldf in enumerate(ldf_grid):
        for i, so2 in enumerate(so2_grid):

            # Update SO2 and LDF parameters
            analyser.params['SO2'].set(value=so2)
            analyser.params['LDF'].set(value=ldf)

            # Extract the parameter list
            fit_params = analyser.params.fittedvalueslist()

            # Create the synthetic spectrum
            spectra[i, ..., j] = analyser.fwd_model(init_fit.grid, *fit_params)

    # Turn off instrument corrections in the analyser
    analyser.stray_flag = False
    analyser.flat_flag = False
    analyser.dark_flag = False
    analyser.despike_flag = False
    analyser.bad_pixels = None

    # Create array to store answers
    ld_results = np.zeros([len(so2_grid)*len(ldf_grid), 6])

    # Reset the SO2 and LDF parameters
    analyser.params['SO2'].set(value=1.0e16)
    analyser.params['LDF'].set(value=0.0)

    # Loop through each synthetic spectrum
    logger.info('Analysing spectra. This can take some time...')
    for j, ldf in enumerate(ldf_grid):
        logger.info(f'- Running LDF = {ldf:.02f}')
        for i, so2 in enumerate(so2_grid):

            # Extract the syntheteic spectrum
            spectrum = [init_fit.grid, spectra[i, ..., j]]

            # Analyse spectrum in waveband 1
            fit1 = analyser.fit_spectrum(spectrum, fit_window=wb1,
                                         update_params=True)

            # Analyse spectrum in waveband 2
            fit2 = analyser.fit_spectrum(spectrum, fit_window=wb2,
                                         update_params=True)

            # Add results to the output array
            row_n = j*len(so2_grid) + i
            ld_results[row_n] = [ldf,
                                 so2,
                                 fit1.params['SO2'].fit_val,
                                 fit1.params['SO2'].fit_err,
                                 fit2.params['SO2'].fit_val,
                                 fit2.params['SO2'].fit_err]

    logger.info('Light dilution calculations complete!')

    return ld_results
