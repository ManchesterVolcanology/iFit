import numpy as np
from scipy.interpolate import griddata
from scipy.optimize import curve_fit

from ifit_lib.make_ils import make_ils
from ifit_lib.smooth import smooth
from ifit_lib.make_poly import make_poly

#========================================================================================
#=========================================fit_spec=======================================
#========================================================================================

def fit_spec(common, y, grid, q = None):
    
    '''
    Function to fit measured spectrum using a full forward model including a solar 
    spectrum background polynomial, ring effect, shift, stretch, and gas amts for so2, 
    no2, o3, bro
    
    INPUTS:
    -------
    common: common dictionary of parameters and variables passed from the main program to 
              subroutines
    y:      intensity data from the measured spectrum 
    grid:   measurement wavelength grid over which the fit occurs
    q:      Queue to which to add the output if threaded (default = None)
    
    OUTPUTS:
    --------
    results:     resulting optimised fit parameters
    cov:         covarience matrix
    y:           processed spectral intensity data (after extracting the desired window 
                   and correction for dark, flat and stray light)
    fitted_flag: boolian variable to tell the main program if the fit was achieved

    '''  
    
    # Unpack the inital fit parameters
    params = np.ones(len(common['params']))
    i = 0
    
    for key, val in common['params'].items():
        params[i] = val
        i += 1
    '''
    # Unpack the solar, ring spectra and gas cross sections and make global for forward
    global sol
    sol = common['sol']
    global ring
    ring = common['ring']
    global so2_xsec
    so2_xsec = common['so2_xsec']
    global no2_xsec
    no2_xsec = common['no2_xsec']
    global o3_xsec
    o3_xsec = common['o3_xsec']
    global bro_xsec
    bro_xsec = common['bro_xsec']
    global model_grid
    model_grid = common['model_grid']
    global poly_n
    poly_n = common['poly_n']
    global ils_flag
    ils_flag = common['ils_flag']
    global set_ils_width
    set_ils_width = common['ils_width']
    global ils_gauss_weight
    ils_gauss_weight = common['ils_gauss_weight']
    global ldf_flag
    ldf_flag = common['ldf_flag']
    global set_ldf
    set_ldf = float(common['ldf'])
    '''
    global com
    com = common

    # Remove the dark spectrum from the measured spectrum
    if common['dark_flag'] == True:
        y = np.subtract(y, common['dark'])
    
    # Remove stray light
    if common['stray_flag'] == True:
        y = np.subtract(y, np.average(y[common['stray_i1']:common['stray_i2']]))
    
    # Find weighting from noise in spectrum
    sigma = np.divide(y, smooth(y, 3))
    
    # Cut desired wavelength window
    y = y[common['ind1']:common['ind2']]
    sigma = sigma[common['ind1']:common['ind2']]

    # Divide by flat spectrum
    if common['flat_flag'] == True:
        y = np.divide(y, common['flat'])

    # Appempt to fit!
    try:
        # Fit
        results, cov = curve_fit(ifit_fwd, grid, y, p0 = params, sigma = sigma)

        # Form so2 transmittance spectrum
        trans_so2 = y / com['F_no_so2']
                        
        # Fit successful
        fitted_flag = True
    
    # If fit fails, report and carry on
    except RuntimeError:
        results = np.zeros(len(params))
        cov = np.zeros((len(params),len(params)))
        trans_so2 = np.zeros(len(y))
        fitted_flag = False

    # Generate a dictionary of errors
    # NOTE this is just variation in the fitting params, this does not include 
    #  systematic errors!
    err_dict = {}

    for m, l in enumerate(common['params'].keys()):
        err_dict[l] = np.sqrt(np.diag(cov))[m]
    
    if q == None:                     
        return results, err_dict, y, trans_so2, com['so2_spec'], fitted_flag
    
    else:
        output = results, err_dict, y, trans_so2, com['so2_spec'], fitted_flag
        q.put(('fit', output))


#========================================================================================
#=========================================ifit_fwd=======================================
#========================================================================================

def ifit_fwd(grid,*args):
    
    '''
    INPUTS:
    -------
    grid:  measurement wavelength grid
    *args: parameters used for fitting. These are n polynomial parameters, wavelength 
             shift and stretch, ring amount and gas amounts for SO2, NO2, O3 and Bro.
             Also optionally the ils width and ldf

    OUTPUTS:
    --------
    F: Fitted spectrum interpolated onto the  spectrometer wavelength grid             
    '''

     # Unpack polynomialparameters
    p = np.zeros(com['poly_n'])
    for i in range(com['poly_n']):
        p[i] = (args[i])

    # Unpack rest
    i += 1
    shift = args[i]
    i += 1
    stretch = args[i]
    i += 1
    ring_amt = args[i]
    i += 1
    so2_amt = args[i]
    i += 1
    no2_amt = args[i]
    i += 1
    o3_amt = args[i]
    i += 1
    bro_amt = args[i]

    # Unpack optional parameters
    if com['ils_flag'] == True:
        i += 1
        ils_width = args[i]
        
    else:
        ils_width = com['ils_width']#set_ils_width
    
    if com['ldf_flag'] == True:
        i += 1
        ldf = args[i]
        i += 1
        #ldf_grad = args[i]
        
    else:
        ldf = com['ldf']#set_ldf

    # Construct background polynomial
    bg_poly = make_poly(com['model_grid'], p)
    
    # Build gas transmittance spectra
    so2_T = np.exp(-(np.multiply(com['so2_xsec'], so2_amt)))
    no2_T = np.exp(-(np.multiply(com['no2_xsec'], no2_amt)))
    o3_T  = np.exp(-(np.multiply(com['o3_xsec'], o3_amt)))
    bro_T  = np.exp(-(np.multiply(com['bro_xsec'], bro_amt)))
    
    # Calculate ring effect
    ring_T = np.multiply(com['sol'], np.multiply(com['ring'], ring_amt))
    sol_T = np.add(com['sol'], ring_T)
    
    # Generate raw, high res spectrum with and without SO2
    raw_F = np.multiply(bg_poly, sol_T)
    raw_F = np.multiply(raw_F, no2_T)
    raw_F = np.multiply(raw_F, o3_T)
    raw_F = np.multiply(raw_F, bro_T)
    raw_F_no_so2 = raw_F              # Create forward model without so2 contribution
    raw_F = np.multiply(raw_F, so2_T)
    
    raw_F_no_o3 = np.multiply(bg_poly, sol_T)
    raw_F_no_o3 = np.multiply(raw_F_no_o3, so2_T)
    raw_F_no_o3 = np.multiply(raw_F_no_o3, no2_T)
    raw_F_no_o3 = np.multiply(raw_F_no_o3, bro_T)
    
    raw_F_no_bro = np.multiply(bg_poly, sol_T)
    raw_F_no_bro = np.multiply(raw_F_no_bro, so2_T)
    raw_F_no_bro = np.multiply(raw_F_no_bro, no2_T)
    raw_F_no_bro = np.multiply(raw_F_no_bro, o3_T)
    
    # Add light dilution effect
    light_d = np.multiply(raw_F_no_so2, ldf)
    raw_F = np.add(raw_F, light_d)
    raw_F_no_so2 = np.add(raw_F_no_so2, light_d)
    raw_F_no_o3  = np.add(raw_F_no_o3, light_d)
    raw_F_no_bro = np.add(raw_F_no_bro, light_d)

    # Avoid unphysical ils widths
    if ils_width < 0.0:
        ils_width = 0.1
    if ils_width > 2.0:
        ils_width = 2.0    
        
    # Convolve high res raw_F with ILS
    ils = make_ils(ils_width, (com['model_grid'][1] - com['model_grid'][0]),
                   com['ils_gauss_weight'])
    F_conv = np.convolve(raw_F, ils, 'same')
    F_conv_no_so2 = np.convolve(raw_F_no_so2, ils, 'same')
    F_conv_no_o3  = np.convolve(raw_F_no_o3 , ils, 'same')
    F_conv_no_bro = np.convolve(raw_F_no_bro, ils, 'same')
    
    # Apply shift and stretch to the model_grid
    shift_model_grid = np.add(com['model_grid'], shift)
    line = np.subtract(com['model_grid'], min(com['model_grid']))
    line = np.divide(line, max(line))
    shift_model_grid = np.add(shift_model_grid, np.multiply(line, stretch))
    
    # Interpolate onto measurement wavelength grid
    F = griddata(shift_model_grid, F_conv, grid)
    
    com['F_no_so2'] = griddata(shift_model_grid, F_conv_no_so2, grid)
    com['so2_spec'] = griddata(shift_model_grid, np.convolve(so2_T, ils, 'same'), grid)

    com['F_no_o3'] = griddata(shift_model_grid, F_conv_no_o3, grid)
    com['o3_spec'] = griddata(shift_model_grid, np.convolve(o3_T, ils, 'same'), grid)

    com['F_no_bro'] = griddata(shift_model_grid, F_conv_no_bro, grid)
    com['bro_spec'] = griddata(shift_model_grid, np.convolve(bro_T, ils, 'same'), grid)
    
    # Remove solar residual
    if com['solar_resid_flag'] == 'Remove':
        F = np.multiply(F, com['solar_resid'])
        com['F_no_so2'] = np.multiply(com['F_no_so2'], com['solar_resid'])
        com['F_no_o3'] = np.multiply(com['F_no_o3'], com['solar_resid'])
        com['F_no_bro'] = np.multiply(com['F_no_bro'], com['solar_resid'])
    
    return F
