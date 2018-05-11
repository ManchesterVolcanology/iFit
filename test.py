# -*- coding: utf-8 -*-
"""
Created on Wed May  9 11:05:30 2018

@author: mqbpwbe2
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

from ifit_lib.read_spectrum import read_spectrum
from ifit_lib.find_nearest import extract_window
from ifit_lib.make_ils import make_ils

def ifit_fwd(grid,a,b,c,d,e,shift,stretch,ring_amt,so2_amt,no2_amt,o3_amt):

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
    raw_F = np.multiply(raw_F, no2_T)
    raw_F = np.multiply(raw_F, o3_T)
    raw_F = np.multiply(raw_F, so2_T)

    # Convolve high res raw_F with ILS
    ils = make_ils(0.52, model_grid[1]-model_grid[0], 1.0)
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


# Create common
common = {}

common['model_resolution'] = 0.02
common['sol_path']         = 'data_bases/sao2010.txt'
common['ring_path']        = 'data_bases/ring.txt'
common['so2_path']         = 'data_bases/SO2_293K.dat'
common['no2_path']         = 'data_bases/No2_223l.dat'
common['o3_path']          = 'data_bases/o3_223l.dat'
common['bro_path']         = 'data_bases/BrO_Cross_298K.txt'

common['wave_start'] = 306
common['wave_stop']  = 319.9 

# Build model grid, a high res grid on which the forward model is build. It extends
#  2 nm beyond the measurement grid and has a spacing controlled by the user
npts = ((common['wave_stop'] + 2) - (common['wave_start'] - 2)) * 50
global model_grid
model_grid = np.linspace(common['wave_start'] - 2, common['wave_stop'] + 2, num = npts + 1)

# Import flat spectrum and extract window of interest
flat_grid, flat = np.loadtxt('data_bases/flat_USB2+H15972.txt', unpack = True)
x, i1, i2 = extract_window(flat_grid, common['wave_start'], common['wave_stop'])
flat = flat[i1:i2]

# Import solar reference spectrum
sol_x, sol_y = np.loadtxt(common['sol_path'], unpack = True)

# Interpolate onto model_grid
sol = griddata(sol_x, sol_y, model_grid)

'''
# Import solar residual spectrum
self.print_output('Importing solar residual spectrum...', add_line = False)
grid, solar_resid = np.loadtxt(settings['resid_path'], unpack = True)
self.print_output('Residual imported', add_line = False)
'''

# Import ring spectrum and interpolate onto the model_grid
ring_x, ring_y = np.loadtxt(common['ring_path'], unpack = True)
ring = griddata(ring_x, ring_y, model_grid)

# Import SO2 data
so2_xsec = np.loadtxt(common['so2_path'], skiprows=1)
so2_xsec = griddata(so2_xsec[:,0], so2_xsec[:,1], model_grid)


# Import NO2 data
no2_xsec = np.loadtxt(common['no2_path'], skiprows=43)
no2_xsec = griddata(no2_xsec[:,0], no2_xsec[:,2], model_grid)


# Import O3 data
o3_xsec = np.loadtxt(common['o3_path'], skiprows=44)
o3_xsec = griddata(o3_xsec[:,0], o3_xsec[:,2], model_grid)


# Import BrO data
bro = np.loadtxt(common['bro_path'], skiprows=12)

# As BrO is in terms of wavenumber(cm^-1), need to convert to nm
bro[:,0] = 10000000/bro[:,0]

# Reverse array so wavelength is accending
bro = bro[::-1]

# Interpolate onto the grid
bro_xsec = griddata(bro[:,0], bro[:,1], model_grid)

# Filepath to iFit output
fpath = 'Results/iFit/2017-05-18/iFit_out.csv'

# Read data file
data = pd.read_csv(fpath)

# Unpack fit parameters
n = 45
spec_fpath = data['File'][n]
a = data['a'][n]
b = data['b'][n]
c = data['c'][n]
d = data['d'][n]
e = data['e'][n]
shift = data['shift'][n]
stretch = data['stretch'][n]
ring_amt = data['ring'][n]
so2_amt = data['so2'][n]
no2_amt = data['no2'][n]
o3_amt = data['o3'][n]

# Read in fitted spectrum
x, y, date, time, spec_no = read_spectrum(spec_fpath, 'IFRiT')
grid, idx0, idx1 = extract_window(x, common['wave_start'], common['wave_stop'])

# Get dark spectrum
darkpath = 'C:/Users/mqbpwbe2/Dropbox (The University of Manchester)/Fieldwork/05 2017 Montserrat/2017-05-18/helitrav1/rt_so2_out (ours)/dark/spectrum_00000.txt'
x, dark, date, time, spec_no = read_spectrum(darkpath, 'IFRiT')

# Remove dark and flat
y = np.subtract(y, dark)
y = y[idx0:idx1]
spec = np.divide(y, flat)



# Form forward model
fit = ifit_fwd(grid,a,b,c,d,e,shift,stretch,ring_amt,so2_amt,no2_amt,o3_amt)

# Plot
plt.plot(grid, spec, label = 'Spectrum')
plt.plot(grid, fit, label = 'Fit')
plt.legend()
plt.show()