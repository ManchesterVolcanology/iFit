import numpy as np
from scipy.interpolate import griddata
from scipy.optimize import curve_fit

from ifit_lib.make_ils import make_ils
from ifit_lib.make_poly import make_poly
from ifit_lib.get_shift import get_shift
from ifit_lib.smooth import smooth

#========================================================================================
#======================================== fit_spec ======================================
#========================================================================================

def fit_spec(common, spectrum, grid, q = None):
    
    '''
    Function to fit measured spectrum using a full forward model including a solar 
    spectrum background polynomial, ring effect, wavelength shift and stretch, and gas
    amounts for so2, no2, o3, bro
    
    INPUTS:
    -------
    common: dictionary
        Common dictionary of parameters and variables passed from the main program
        to subroutines
        
    spectrum: 2D array
        Intensity data from the measured spectrum 
        
    grid: 1D array
        Measurement wavelength grid over which the fit occurs
        
    q: queue (optional)
        Queue to which to add the output if threaded (default = None)
    
    OUTPUTS:
    --------
    fit_dict: dictionary
        Dictionary of optimised parameters
        
    err_dict: dictionary
        Dictionary of error in the optimised parameters
        
    y: array 
        Measured spectrum, corrected for dark, bias and flat response, in the fitting 
        window
        
    fit: array
        Fitted spectrum
        
    gas_T: dictionary
        Dictionary of synthetic and measured absorbance spectra
        
    fitted_flag: bool
        Flag showing if the fit was successful or not

    ''' 
            
    # Cretae a copy of common for forward model to access
    global com
    com = common

    # Unpack spectrum
    x, y = spectrum
    
    # Pre-calculate shift
    if common['Fit shift'] == 'Pre-calc':
        
        # Find first guess for the shift
        shift, err = get_shift(x, y, common) 
        
        # If sensible, update shift value
        if err == 0:
            common['params']['shift'][0] = shift 

    # Remove the dark spectrum from the measured spectrum
    if common['dark_flag'] == True:
        y = np.subtract(y, common['dark'])
    
    # Remove stray light
    if common['stray_flag'] == True:
        y = np.subtract(y, np.average(y[common['stray_idx']]))
    
    # Cut desired wavelength window
    y = y[common['fit_idx']]

    # Divide by flat spectrum
    if common['flat_flag'] == True:
        y = np.divide(y, common['flat'])
        
    # Set weighting for fitting function
    if common['fit_weight'] == 'SO2 xsec':
        so2_poly = np.polyval(np.polyfit(com['model_grid'], com['so2_xsec'], 3),
                              com['model_grid'])
        sigma = np.divide(1, np.abs(np.subtract(common['so2_xsec'], so2_poly)))
        sigma = griddata(com['model_grid'], sigma, grid, 'cubic')
        
    if common['fit_weight'] == 'Noise':
        noise = np.abs(np.divide(y, smooth(y, 3)))
        sigma = np.divide(1, noise)
        
        sigma = None
        
    if common['fit_weight'] == 'None':
        sigma = None
    
    # Unpack the inital fit parameters
    fit_params = []
    
    for key, val in common['params'].items():
       
        if val[1] == 'Fit':
            fit_params.append(val[0])
    
    # Create dictionary to hold gas trans data
    gas_T = {}
    
    # Appempt to fit!
    try:
        # Fit
        popt, pcov = curve_fit(ifit_fwd_model, grid, y, p0 = fit_params, sigma = sigma)

        # Calculate transmittance spectra
        fit = ifit_fwd_model(grid, *popt, calc_trans_flag = True)
        
        # Form absorbance spectra
        for g in ['so2', 'no2', 'o3', 'bro', 'ring']:
            gas_T['meas_abs_' + g] = -np.log(np.divide(y, com['F_no_' + g]))
            gas_T['synth_abs_' + g]  = -np.log(com[g + '_spec'])
        
        #################################################################################
        '''
        np.savetxt('Dump/spec.txt', np.column_stack((grid, y)))
        np.savetxt('Dump/fit.txt', np.column_stack((grid, fit)))
        
        for key, val in gas_T.items():
            np.savetxt('Dump/' + key + '.txt', np.column_stack((grid, gas_T[key])))
        '''
        #################################################################################
                        
        # Fit successful
        fitted_flag = True
    
    # If fit fails, report and carry on
    except RuntimeError:
        # Fill returned arrays with zeros
        popt = np.zeros(len(fit_params))
        pcov = np.zeros((len(fit_params),len(fit_params)))
        
        fit = np.zeros(len(grid))
        
        # Form absorbance spectra
        for g in ['so2', 'no2', 'o3', 'bro', 'ring']:
            gas_T['meas_abs_' + g] = np.zeros(len(grid))
            gas_T['synth_abs_' + g]  = np.zeros(len(grid))
            
        # Turn off fitted flag
        fitted_flag = False

    # Unpack fit results
    fit_dict  = {}
    m = 0
    for key, val in common['params'].items():
        if val[1] == 'Fit':
            fit_dict[key] = popt[m]
            m+=1

    # Generate a dictionary of errors
    # NOTE this is covarience in fitting params and does not include systematic errors!
    err_dict = {}
    m = 0
    for key, val in common['params'].items():
        if val[1] == 'Fit':
            err_dict[key] = np.sqrt(np.diag(pcov))[m]
            m+=1
    
    # Return results, either to a queue if threaded, or as an array if not
    if q == None:                     
        return fit_dict, err_dict, y, fit, gas_T, fitted_flag
    
    else:
        output = fit_dict, err_dict, y, fit, gas_T, fitted_flag
        q.put(('fit', output))


#========================================================================================
#======================================== ifit_fwd ======================================
#========================================================================================

def ifit_fwd_model(grid, *fit_params, calc_trans_flag = False):
    
    '''
    iFit forward model to fit measured UV sky spectra
    
    INPUTS:
    -------
    grid, array
        Measurement wavelength grid
        
    *args, list
        Forward model state vector 
        
    calc_trans_flag, bool (optional)
        Flag whether or not to calculate transmission spectra. Not required for fitting,
        but helpful for analysis of fit quality

    OUTPUTS:
    --------
    F, array
        Fitted spectrum interpolated onto the spectrometer wavelength grid             
    '''

    # Unpack params
    p = {}
    i = 0

    for key, val in com['params'].items():
           
        if val[1] == 'Fit':
            p[key] = fit_params[i]
            i += 1
            
        if val[1] in ['Fix', 'Pre-calc', 'File']:
            p[key] = val[0]
            
        if val[1] == 'N/A':
            p[key] = 0
       
    # Unpack polynomial parameters
    poly_coefs = np.zeros(com['poly_n'])
    for i in range(com['poly_n']):
        poly_coefs[i] = (fit_params[i])

    # Construct background polynomial
    bg_poly = make_poly(com['model_grid'], poly_coefs)
    
    # Build gas transmittance spectra
    so2_T = np.exp(-(np.multiply(com['so2_xsec'], p['so2_amt'])))
    no2_T = np.exp(-(np.multiply(com['no2_xsec'], p['no2_amt'])))
    o3_T  = np.exp(-(np.multiply(com['o3_xsec'],  p['o3_amt'])))
    bro_T = np.exp(-(np.multiply(com['bro_xsec'], p['bro_amt'])))

    
    # Background sky spectrum
    sol_spec = np.multiply(bg_poly, com['sol'])
    
    
    # Ring effect
    ring_T = np.exp(np.multiply(com['ring'], p['ring_amt']))
    sol_T = np.multiply(sol_spec, ring_T)
    
    
    # Background sky spectrum
    bg_spec = np.multiply(np.multiply(sol_T, no2_T), o3_T)
       
    
    # Include gasses and areosol
    raw_F = np.multiply(np.multiply(bg_spec, bro_T), so2_T)
    plume_F = raw_F
    
    # Calc light dilution effect
    light_d = np.multiply(bg_spec, p['ldf'])
    
    # Combine with normal spectrum
    plume_F = np.multiply(plume_F, 1-p['ldf'])
    raw_F = np.add(plume_F, light_d)
    
    # Make sure ILS is sensible, then form and convolve with high res raw_F
    ils = make_ils(p['ils_width'], (com['model_grid'][1] - com['model_grid'][0]),
                   com['gauss_weight'])
    F_conv = np.convolve(raw_F, ils, 'same')

    
    # Apply shift and stretch to the model_grid
    shift_model_grid = np.add(com['model_grid'], p['shift'])
    line = np.linspace(0, 1, num = len(shift_model_grid))
    shift_model_grid = np.add(shift_model_grid, np.multiply(line, p['stretch']))
    
    # Interpolate onto measurement wavelength grid
    F = griddata(shift_model_grid, F_conv, grid, method = 'cubic')
    
    if calc_trans_flag == True:
        
        # Generate raw spectra without selected gas
        raw_F_no_so2 = np.multiply(sol_T, no2_T)
        raw_F_no_so2 = np.multiply(raw_F_no_so2, o3_T)
        raw_F_no_so2 = np.multiply(raw_F_no_so2, bro_T)  
        
        raw_F_no_no2 = np.multiply(sol_T, no2_T)
        raw_F_no_no2 = np.multiply(raw_F_no_no2, o3_T)
        raw_F_no_no2 = np.multiply(raw_F_no_no2, bro_T) 
    
        raw_F_no_o3 = np.multiply(sol_T, so2_T)
        raw_F_no_o3 = np.multiply(raw_F_no_o3, no2_T)
        raw_F_no_o3 = np.multiply(raw_F_no_o3, bro_T)
        
        raw_F_no_bro = np.multiply(sol_T, so2_T)
        raw_F_no_bro = np.multiply(raw_F_no_bro, no2_T)
        raw_F_no_bro = np.multiply(raw_F_no_bro, o3_T)
        
        # Generate raw spectrum without ring
        raw_F_no_ring = np.multiply(sol_spec, so2_T)
        raw_F_no_ring = np.multiply(raw_F_no_ring, no2_T)
        raw_F_no_ring = np.multiply(raw_F_no_ring, o3_T)
        raw_F_no_ring = np.multiply(raw_F_no_ring, bro_T)
        
        # Add light dilution
        raw_F_no_so2 = np.add(np.multiply(raw_F_no_so2, 1 - p['ldf']), light_d)
        raw_F_no_no2 = np.add(np.multiply(raw_F_no_no2, 1 - p['ldf']), light_d)
        raw_F_no_o3  = np.add(np.multiply(raw_F_no_o3,  1 - p['ldf']), light_d)
        raw_F_no_bro = np.add(np.multiply(raw_F_no_bro, 1 - p['ldf']), light_d)
        raw_F_no_ring= np.add(np.multiply(raw_F_no_ring,1 - p['ldf']), light_d)
        
        # Convolve with the ils
        F_conv_no_so2 = np.convolve(raw_F_no_so2, ils, 'same')
        F_conv_no_no2 = np.convolve(raw_F_no_no2, ils, 'same')
        F_conv_no_o3  = np.convolve(raw_F_no_o3,  ils, 'same')
        F_conv_no_bro = np.convolve(raw_F_no_bro, ils, 'same')   
        F_conv_no_ring= np.convolve(raw_F_no_ring,ils, 'same')     
        
        # Place onto spectrometer grid
        com['F_no_so2'] = griddata(shift_model_grid, F_conv_no_so2, grid)
        com['so2_spec'] = griddata(shift_model_grid, np.convolve(so2_T,ils,'same'), grid, 
                                   method = 'cubic')
    
        com['F_no_no2'] = griddata(shift_model_grid, F_conv_no_no2, grid)
        com['no2_spec'] = griddata(shift_model_grid, np.convolve(no2_T,ils,'same'), grid, 
                                   method = 'cubic')
    
        com['F_no_o3'] = griddata(shift_model_grid, F_conv_no_o3, grid)
        com['o3_spec'] = griddata(shift_model_grid, np.convolve(o3_T, ils, 'same'), grid, 
                                  method = 'cubic')
    
        com['F_no_bro'] = griddata(shift_model_grid, F_conv_no_bro, grid)
        com['bro_spec'] = griddata(shift_model_grid, np.convolve(bro_T,ils,'same'), grid, 
                                   method = 'cubic')
    
        com['F_no_bro'] = griddata(shift_model_grid, F_conv_no_bro, grid)
        com['bro_spec'] = griddata(shift_model_grid, np.convolve(bro_T,ils,'same'), grid, 
                                   method = 'cubic')
        
        com['F_no_ring'] = griddata(shift_model_grid, F_conv_no_ring, grid)
        com['ring_spec'] = griddata(shift_model_grid, np.convolve(ring_T,ils,'same'), 
                                    grid, method = 'cubic')

    # Remove solar residual
    if com['solar_resid_flag'] == 'Remove':
        F = np.multiply(F, com['solar_resid'])
        
        if calc_trans_flag == True:
            com['F_no_so2'] = np.multiply(com['F_no_so2'], com['solar_resid'])
            com['F_no_no2'] = np.multiply(com['F_no_no2'], com['solar_resid'])
            com['F_no_o3']  = np.multiply(com['F_no_o3'],  com['solar_resid'])
            com['F_no_bro'] = np.multiply(com['F_no_bro'], com['solar_resid'])
            com['F_no_ring']= np.multiply(com['F_no_ring'],com['solar_resid'])
        
    return F
