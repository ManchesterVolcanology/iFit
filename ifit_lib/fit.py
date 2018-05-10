import numpy as np
from scipy.interpolate import griddata
from scipy.optimize import curve_fit

from ifit_lib.make_ils import make_ils
from ifit_lib.smooth import smooth

#========================================================================================
#=========================================fit_spec=======================================
#========================================================================================

# Function to fit measured spectrum using a full forward model including a solar spectrum
#  background polynomial, ring effect, shift, stretch, and gas amts for so2, no2, o3, bro

# INPUTS: common; common dictionary of parameters and variables passed from the main 
#                  program to subroutines
#         y; intensity data from the measured spectrum 
#         grid; measurement wavelength grid over which the fit occurs

# OUTPUTS: results; resulting optimised fit parameters
#          cov; covarience matrix
#          y; processed spectral intensity data (after extracting the desired window and 
#              correction for dark, flat and stray light)
#          fitted_flag; boolian variable to tell the main program if the fit was achieved

def fit_spec(common, y, grid, fwd_model):
    
    # Unpack the inital fit parameters
    params = common['params']
    
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
    global ils
    ils = common['ils']
    global ils_gauss_weight
    ils_gauss_weight = common['ils_gauss_weight']
    #global solar_resid
    #solar_resid = common['solar_resid']
    global ldf
    ldf = float(common['ldf'])

    # If solar flag is false, set resid spectrum to ones
    #if settings['solar_resid_flag'] != 'Remove':
    #    solar_resid = 1
    
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

    # Unpack initial parameter names and values
    names = []
    init_vals = []

    for i in params:
        names.append(i[0])
        init_vals.append(i[1])
     
    # Appempt to fit!
    try:
        results, cov = curve_fit(fwd_model, grid, y, p0 = init_vals, sigma = sigma)
                                 
        # Fit successful
        fitted_flag = True
    
    # If fit fails, report and carry on
    except RuntimeError:
        results = np.zeros(len(init_vals))
        cov = np.zeros((len(init_vals),len(init_vals)))
        fitted_flag = False

    # Generate a dictionary of errors
    # NOTE this is just variation in the fitting params, this does not include 
    #  systematic errors!
    err_dict = {}
    m = 0
    for l in common['params']:
        err_dict[l[0]] = np.sqrt(np.diag(cov))[m]
        m += 1
                        
    return results, err_dict, y, fitted_flag
                    
#========================================================================================
#=======================================gen_fit_norm=====================================
#========================================================================================
                    
def gen_fit_norm(grid, fit_params):
    
    # Unpack fit results
    a,b,c,d,e,shift,stretch,ring_amt,so2_amt,no2_amt,o3_amt = fit_params
    
    # Feed into forward model to recreate the fit
    fit = ifit_fwd(grid,a,b,c,d,e,shift,stretch,ring_amt,so2_amt,no2_amt,o3_amt)
    
    return fit

#========================================================================================
#=========================================ifit_fwd=======================================
#========================================================================================

# Function describing the forward model used to fit the measured spectra

# INPUTS: grid; measurement wavelength grid over which the fit occurs
#         other params; fit paramters as decribed in the main program

# OUTPUTS: F; constructed forward model

def ifit_fwd(grid,a,b,c,d,e,shift,stretch,ring_amt,so2_amt,no2_amt,o3_amt):

    # Construct background polynomial
    bg_poly = np.polyval((a,b,c,d), model_grid)
    
    # Build gas transmittance spectra
    so2_T = np.exp(-(np.multiply(so2_xsec, so2_amt)))
    no2_T = np.exp(-(np.multiply(no2_xsec, no2_amt)))
    o3_T  = np.exp(-(np.multiply(o3_xsec, o3_amt)))
    
    # Calculate ring effect
    ring_T = np.multiply(sol, np.multiply(ring, ring_amt))
    sol_T = np.add(sol, ring_T)
    
    # Generate raw, high res spectrum with and without SO2
    raw_F = np.multiply(bg_poly, sol_T)
    raw_F = np.multiply(raw_F, no2_T)
    raw_F = np.multiply(raw_F, o3_T)
    raw_F_no_so2 = raw_F              # Create forward model without so2 contribution
    raw_F = np.multiply(raw_F, so2_T)
    
    # Add light dilution effect
    light_d = np.multiply(raw_F_no_so2, ldf)
    raw_F = np.add(raw_F, light_d)

    # Convolve high res raw_F with ILS
    F_conv = np.convolve(raw_F, ils, 'same')
    
    # Apply shift and stretch to the model_grid
    shift_model_grid = np.add(model_grid, shift)
    line = np.subtract(model_grid, min(model_grid))
    line = np.divide(line, max(line))
    shift_model_grid = np.add(shift_model_grid, np.multiply(line, stretch))
    
    # Interpolate onto measurement wavelength grid
    F = griddata(shift_model_grid, F_conv, grid)
    
    # Remove solar residual
    #F = np.multiply(F, solar_resid)
    
    return F



#========================================================================================
#========================================================================================
#=========================================with ils=======================================
#========================================================================================
#========================================================================================
                    
#========================================================================================
#=======================================gen_fit_ils======================================
#========================================================================================
                    
def gen_fit_ils(grid, fit_params):
    
    # Unpack fit results
    a,b,c,d,e,shift,stretch,ring_amt,so2_amt,no2_amt,o3_amt,ils_width = fit_params
    
    # Feed into forward model to recreate the fit
    fit = ifit_ils_fwd(grid, a, b, c, d, e, shift, stretch, ring_amt, so2_amt, no2_amt, 
                       o3_amt, ils_width)
    
    return fit

#========================================================================================
#=======================================ifit_ils_fwd=====================================
#========================================================================================

# Function describing the forward model used to fit the measured spectra

# INPUTS: grid; measurement wavelength grid over which the fit occurs
#         other params; fit paramters as decribed in the main program

# OUTPUTS: F; constructed forward model

def ifit_ils_fwd(grid,a,b,c,d,e,shift,stretch,ring_amt,so2_amt,no2_amt,o3_amt,ils_width):
    
    # Construct background polynomial
    bg_poly = np.polyval((a,b,c,d,e), model_grid)
    
    # Build gas transmittance spectra
    so2_T = np.exp(-(np.multiply(so2_xsec, so2_amt)))
    no2_T = np.exp(-(np.multiply(no2_xsec, no2_amt)))
    o3_T  = np.exp(-(np.multiply(o3_xsec, o3_amt)))
    
    # Calculate ring effect
    ring_T = np.multiply(sol, np.multiply(ring, ring_amt))
    sol_T = np.add(sol, ring_T)
    
    # Generate raw, high res spectrum with and without SO2
    raw_F = np.multiply(bg_poly, sol_T)
    raw_F = np.multiply(raw_F, o3_T)
    raw_F = np.multiply(raw_F, no2_T)
    raw_F_no_so2 = raw_F              # Create forward model without so2 contribution
    raw_F = np.multiply(raw_F, so2_T)
    
    # Add light dilution effect
    light_d = np.multiply(raw_F_no_so2, ldf)
    raw_F = np.add(raw_F, light_d)
  
    # Convolve high res raw_F with ILS
    ils = make_ils(ils_width, model_grid[1] - model_grid[0], ils_gauss_weight)
    F_conv = np.convolve(raw_F, ils, 'same')
    F_no_so2_conv = np.convolve(raw_F_no_so2, ils, 'same')
    
    # Apply shift and stretch to the model_grid
    shift_model_grid = np.add(model_grid, shift)
    line = np.subtract(model_grid, min(model_grid))
    line = np.divide(line, max(line))
    shift_model_grid = np.add(shift_model_grid, np.multiply(line, stretch))
    
    # Interpolate 
    F = griddata(shift_model_grid, F_conv, grid)
    global F_no_so2
    F_no_so2 = griddata(shift_model_grid, F_no_so2_conv, grid)
    global so2_spec
    so2_spec = griddata(shift_model_grid, so2_T, grid)
    
    # Remove solar residual
    #F = np.multiply(F, solar_resid)
    #F_no_so2 = np.multiply(F_no_so2, solar_resid)
    
    return F

#========================================================================================
#========================================================================================
#=========================================with ldf=======================================
#========================================================================================
#========================================================================================

#========================================================================================
#=======================================gen_fit_ldf======================================
#========================================================================================
                    
def gen_fit_ldf(grid, fit_params):
    
    # Unpack fit results
    a,b,c,d,e,shift,stretch,ring_amt,so2_amt,no2_amt,o3_amt,ldf = fit_params
    
    # Feed into forward model to recreate the fit
    fit = ifit_ldf_fwd(grid, a, b, c, d, e, shift, stretch, ring_amt, so2_amt, no2_amt, 
                       o3_amt, ldf)
    
    return fit


#========================================================================================
#=======================================ifit_ldf_fwd=====================================
#========================================================================================

# Function describing the forward model used to fit the measured spectra

# INPUTS: grid; measurement wavelength grid over which the fit occurs
#         other params; fit paramters as decribed in the main program

# OUTPUTS: F; constructed forward model

def ifit_ldf_fwd(grid,a,b,c,d,e,shift,stretch,ring_amt,so2_amt,no2_amt,o3_amt,ldf):

    # Construct background polynomial
    bg_poly = np.polyval((a,b,c,d,e), model_grid)
    
    # Build gas transmittance spectra
    so2_T = np.exp(-(np.multiply(so2_xsec, so2_amt)))
    no2_T = np.exp(-(np.multiply(no2_xsec, no2_amt)))
    o3_T  = np.exp(-(np.multiply(o3_xsec,  o3_amt)))
    
    # Calculate ring effect
    ring_T = np.multiply(sol, np.multiply(ring, ring_amt))
    sol_T = np.add(sol, ring_T)

    # Generate raw, high res spectrum with and without SO2
    raw_F = np.multiply(bg_poly, sol_T)
    raw_F = np.multiply(raw_F, o3_T)
    raw_F = np.multiply(raw_F, no2_T)
    raw_F_no_so2 = raw_F              # Create forward model without so2 contribution
    raw_F = np.multiply(raw_F, so2_T)
    
    # Add light dilution effect
    #if ldf < 0.0 or ldf > 1.0:
    #    ldf = 0.0
    light_d = np.multiply(raw_F_no_so2, ldf)
    raw_F = np.add(raw_F, light_d)
  
    # Convolve high res raw_F with ILS
    F_conv = np.convolve(raw_F, ils, 'same')
    F_no_so2_conv = np.convolve(raw_F_no_so2, ils, 'same')
    
    # Apply shift and stretch to the model_grid
    shift_model_grid = np.add(model_grid, shift)
    line = np.subtract(model_grid, min(model_grid))
    line = np.divide(line, max(line))
    shift_model_grid = np.add(shift_model_grid, np.multiply(line, stretch))
    
    # Interpolate onto measurement wavelength grid
    F = griddata(shift_model_grid, F_conv, grid)
    global F_no_so2
    F_no_so2 = griddata(shift_model_grid, F_no_so2_conv, grid)
    global so2_spec
    so2_spec = griddata(shift_model_grid, so2_T, grid)
    
    # Remove solar residual
    #F = np.multiply(F, solar_resid)
    #F_no_so2 = np.multiply(F_no_so2, solar_resid)
    
    return F
