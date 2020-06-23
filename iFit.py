# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 19:01:27 2020

@author: mqbpwbe2
"""

import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from ifit.parameters import Parameters
from ifit.spectral_analysis import Analyser
from ifit.load_spectra import read_spectrum, average_spectra

# =============================================================================
# Define analysis files
# =============================================================================

# Master file path to spectra
mpath = 'C:/Users/mqbpwbe2/Documents/UV_spectra/Etna/2007_01_02_Etna/070102_etna.'
        
# Set locations of dark and measurement spectra
dark_path = mpath
meas_path = mpath

# Set the location to save the spectra blocks
save_path = 'iFit_out.csv'

# Set the spectra type
spec_type = 'Master.Scope'

# Set the dark, reference and measurement sectra numbers
dark_num = np.arange(6, 30)
meas_num = np.arange(50, 1000)

# Build the spectra filenames
print('Building file paths...')
dark_fnames = [f'{dark_path}{n:05d}.Master.Scope' for n in dark_num 
               if os.path.isfile(f'{dark_path}{n:05d}.Master.Scope')]
meas_fnames = [f'{meas_path}{n:05d}.Master.Scope' for n in meas_num 
               if os.path.isfile(f'{meas_path}{n:05d}.Master.Scope')]
print('Done!\n')

# Control whether plotting is turned on or off
plotting_flag = True

# =============================================================================
# Parameter Setup
# =============================================================================

# Create parameter dictionary
params = Parameters()

# Add the gases
params.add('SO2',  value=1.0e16, vary=True, xpath='Ref/SO2_295K.txt')
params.add('O3',   value=1.0e19, vary=True, xpath='Ref/O3_243K.txt')
params.add('Ring', value=0.1,    vary=True, xpath='Ref/Ring_DOASIS.txt')
    
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

# Add ILS parameters
params.add('fwem', value=1.4, vary=True)
params.add('k',    value=2.0, vary=True)
params.add('a_w',  value=0.0, vary=False)
params.add('a_k',  value=0.0, vary=False)

# Generate the analyser
analyser = Analyser(params,
                    fit_window=[310, 320],
                    frs_path = 'Ref/sao2010.txt')

print(params.pretty_print(cols='all'))

# Read in the dark spectrum
x, dark = average_spectra(dark_fnames, spec_type)
analyser.dark_spec = dark

print('Done!\n')

# =============================================================================
# Generate the figure
# =============================================================================

if plotting_flag:
    
    print('Making plot canvas...')

    # Make the figure and define the subplot grid
    fig = plt.figure(figsize=[10,6.4])
    gs = GridSpec(2,2)
    
    # Define axes
    ax0 = fig.add_subplot(gs[0,0])
    ax1 = fig.add_subplot(gs[1,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,1])
    
    # Define plot lines
    l0, = ax0.plot([], [], 'C0x-') # Measured spectrum
    l1, = ax0.plot([], [], 'C1-')  # Model fit
    
    l2, = ax1.plot([], [], 'C0x-') # Residual
    
    l3, = ax2.plot([], [], 'C0x-') # Measured OD
    l4, = ax2.plot([], [], 'C1-')  # Fit OD
    
    l5, = ax3.plot([], [], 'C2o-') # Time series
    
    print('Done!\n')

# =============================================================================
# Run analysis
# =============================================================================

print('Beginning analysis...')

# Make a list of column names
cols = ['Number', 'Time']
for par in params:
    cols += [par, f'{par}_err']
cols += ['fit_quality', 'int_lo', 'int_hi', 'int_av']

# Make a dataframe to hold the fit results
df = pd.DataFrame(index=np.arange(len(meas_fnames)), columns=cols)

for n, fname in tqdm(enumerate(meas_fnames)):
    
    x, y, spec_info, read_err = read_spectrum(fname, spec_type)
    
    fit = analyser.fit_spectrum([x,y], update_params=True, 
                                interp_method='linear',
                                calc_od=['SO2'])

    # Add to the results dataframe
    row = [n, spec_info['time']]
    for par in fit.params.values():
        row += [par.fit_val, par.fit_err]
    row += [fit.nerr, fit.int_lo, fit.int_hi,
            fit.int_av]
    df.loc[n] = row
    
    if plotting_flag:
        # Update the plot
        l0.set_data(fit.grid, fit.spec)
        l1.set_data(fit.grid, fit.fit)
        l2.set_data(fit.grid, fit.resid)
        l3.set_data(fit.grid, fit.meas_od['SO2'])
        l4.set_data(fit.grid, fit.synth_od['SO2'])
        l5.set_data(df['Number'], df['SO2'])
        
        for ax in [ax0, ax1, ax2, ax3]:
            ax.relim()
            ax.autoscale_view()
        
        plt.pause(0.01)
        if n == 0:
            plt.tight_layout()
    
if plotting_flag: plt.show()
print('Done!\n')  

print('Saving results...')
df.to_csv(save_path)
print('Done!')