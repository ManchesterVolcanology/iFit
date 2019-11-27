# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 16:36:31 2019

@author: mqbpwbe2
"""

import glob
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation

from ifit.model_setup import model_setup
from ifit.spectral_analysis import pre_process, fit_spectrum
from ifit.parameters import Parameters

#from openso2.analyse_scan import read_scan


mode = 'single'

#=================================== Single ===================================
if mode == 'single':
    spec_path = 'Data/Masaya/spectrum_00385.txt'
    spec_path = 'Data/Masaya/spectrum_00430.txt'
    dark_path = 'Data/Masaya/dark_00000.txt'

    spec_path = 'F:/Writing/2019_iFit/JVGR/Data/Masaya Traverse/USB2000+/spectra/spectrum_00380.txt'
    dark_path = 'F:/Writing/2019_iFit/JVGR/Data/Masaya Traverse/USB2000+/dark/spectrum_00000.txt'

#=================================== Multi ====================================
if mode == 'multi':
    dark_path = 'Data/Masaya/dark_00000.txt'
    spec_fnames = glob.glob('Data/Masaya/spectrum*')

    spec_fnames = glob.glob('F:/Writing/2019_iFit/JVGR/Data/Masaya Traverse/USB2000+/spectra/spectrum*')
    dark_path = 'F:/Writing/2019_iFit/JVGR/Data/Masaya Traverse/USB2000+/dark/spectrum_00000.txt'

    spec_fnames.sort()

#==================================== Scan ====================================
if mode == 'scan':

    fpath = ''



od_gas = ['SO2', 'Ring', 'O3']

#==============================================================================
#=============================== Model Settings ===============================
#==============================================================================

spec_name = 'FLMS02101'
spec_name = 'USB2+H15972'

# Set up the program settings
settings = {'w_lo':          310.0,
            'w_hi':          320.0,
            'model_spacing': 0.01,
            'model_padding': 1.0,
            'flat_flag':     True,
            'stray_flag':    True,
            'dark_flag':     True,
            'ils_type':      'File',
            'frs_path':      'Ref/sao2010.txt',
            'flat_path':     f'Spectrometer/flat_{spec_name}.txt',
            'ils_path':      f'Spectrometer/ils_params_{spec_name}.txt'
            }

# Set the gas data
gas_data = {}
gas_data['SO2']  = ['Ref/SO2_298K_shifted.txt', 1.0e18, True]
gas_data['NO2']  = ['Ref/NO2_223K.txt',         1.0e18, True]
gas_data['O3']   = ['Ref/O3_253K.txt',          1.0e18, True]
gas_data['Ring'] = ['Ref/Ring.txt',             0.1,    True]
#gas_data['O3_293'] =  ['Ref/O3_293K.txt',  1.0e17, True]

# Add to the settings dictionary
settings['gas_data'] = gas_data

# Build the forward model
common = model_setup(settings)

#==============================================================================
#============================== Parameters setup ==============================
#==============================================================================

# Create parameter dictionary
params = Parameters()

# Set other model settings
bg_poly_n = 4
bl_poly_n = 1
wl_poly_n = 2

for i in range(bg_poly_n):
    params.add(f'bg_poly{i}', value = 1.0, vary = True)
for i in range(bl_poly_n):
    params.add(f'offset{i}', value = 1.0, vary = True)
for i in range(wl_poly_n):
    params.add(f'shift{i}', value = -0.2, vary = True)

# Add the gases
for gas in gas_data:
    params.add(gas, value = gas_data[gas][1], vary = gas_data[gas][2])

# Add to the common
common['params'] = params
print(params.pretty_print())

#==============================================================================
#=================================== Single ===================================
#==============================================================================

if mode == 'single':

    # Read in the spectrum and dark
    x, common['dark'] = np.loadtxt(dark_path, unpack = True)
    meas_spectrum = np.loadtxt(spec_path, unpack = True)

    # Find the fit and stray windows
    common['fit_idx'] = np.where(np.logical_and(x > settings['w_lo'],
                                                x < settings['w_hi']))

    common['stray_idx'] = np.where(np.logical_and(x > 280,  x < 290))

    # Pre-process the spectrum before the fit
    spectrum = pre_process(meas_spectrum, common)
    grid, spec = spectrum

    # Fit the spectrum
    fit_result = fit_spectrum(spectrum, common, calc_od=od_gas)
    print(fit_result.print_result())

    # Save the outputs
    np.savetxt('dump/spectrum.txt', np.column_stack(grid, spec))
    np.savetxt('dump/fit.txt',      fit_result.fit)
    np.savetxt('dump/resid.txt',    fit_result.resid)
    np.savetxt('dump/so2.txt',      fit_result.resid)
    np.savetxt('dump/o3.txt',       fit_result.resid)
    np.savetxt('dump/ring.txt',     fit_result.resid)

    # Set up figure
    fig = plt.figure(figsize = [12,8])
    gs = gridspec.GridSpec(2,2)
    ax0 = fig.add_subplot(gs[0,0])
    ax1 = fig.add_subplot(gs[1,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,1])

    ax0.plot(grid, spec, 'C0o-', label = 'Data')
    ax0.plot(grid, fit_result.fit, 'C1-', lw = 2, label = 'Fit')
    ax1.plot(grid, fit_result.resid, 'C0o-')
    ax2.plot(meas_spectrum[0], meas_spectrum[1])
    for gas in od_gas:
        ax3.plot(grid, fit_result.meas_od[gas], 'C0o-')
        ax3.plot(grid, fit_result.synth_od[gas], 'C1-', lw=2)

    ax0.legend(loc = 'upper left')
    ax0.set_ylabel('Instensity (counts)', fontsize = 14)

    ax1.set_xlabel('Wavelength (nm)', fontsize = 14)
    ax1.set_ylabel('Residual (%)', fontsize = 14)
    ax1.set_ylim(-3, 3)

    plt.tight_layout()
    plt.show()

#==============================================================================
#=================================== Multi ====================================
#==============================================================================

if mode == 'multi':

    # Make a list of column names
    cols = ['File', 'Number', 'Time']
    for par in params:
        cols += [par, f'{par}_err']
    cols += ['fit_quality']

    n_spec = np.arange(len(spec_fnames))

    # Make a dataframe to hold the fit results
    results_df = pd.DataFrame(index = n_spec, columns = cols)

    # Read in the spectrum and dark
    x, common['dark'] = np.loadtxt(dark_path, unpack = True)

    # Find the fit and stray windows
    common['fit_idx'] = np.where(np.logical_and(x > settings['w_lo'],
                                                x < settings['w_hi']))
    if common['stray_flag']:
        common['stray_idx'] = np.where(np.logical_and(x > 280,  x < 290))

    spec_num = np.arange(len(spec_fnames))
    so2 = np.full(len(spec_fnames), np.nan)
    err = np.full(len(spec_fnames), np.nan)

    # Set up figure
    fig = plt.figure(figsize = [8,8])
    gs = gridspec.GridSpec(3,1)
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])
    ax2 = fig.add_subplot(gs[2])

    axes = [ax0, ax1, ax2]

    def animate(i):

        # Unpack the spectrum
        fname = spec_fnames[i]
        x, y = np.loadtxt(fname, unpack = True)

        # Pre-process the spectrum before the fit
        spectrum = pre_process([x, y], common)
        grid, spec = spectrum

        # Fit the spectrum
        fit_result = fit_spectrum(spectrum, common, update_params = True)
        #fit_result.print_result()

        # Add the the results dataframe
        row = [fname, i, dt.datetime.now()]
        for par in fit_result.params.values():
            row += [par.fit_val, par.fit_err]
        row += [fit_result.nerr]

        results_df.loc[i] = row

        so2[i] = fit_result.params['SO2'].fit_val
        err[i] = fit_result.params['SO2'].fit_err

        # Clear the axes
        for ax in axes:
            ax.clear()

        # Plot the data
        ax0.plot(grid, spec, 'C0o', label = 'Data')
        ax0.plot(grid, fit_result.fit, 'C1-', lw = 2, label = 'Fit')
        ax1.plot(grid, fit_result.resid, 'C0o-')
        ax2.errorbar(spec_num, so2, yerr = err, color='g', capsize=2,
                     marker = 'o', ms = 3, alpha = 0.5)

        ax0.legend(loc = 'upper left')
        ax0.set_ylabel('Instensity (counts)', fontsize = 14)

        ax1.set_xlabel('Wavelength (nm)', fontsize = 14)
        ax1.set_ylabel('Residual (%)', fontsize = 14)
        ax1.set_ylim(-3, 3)

        ax2.set_xlabel('Spectrum Number', fontsize = 14)
        ax2.set_ylabel('SO$_2$ SCD (ppm.m)', fontsize = 14)
        ax2.set_xlim(-10, len(spec_fnames)+10)

        plt.tight_layout()

        # Rescale the axes
        for ax in axes:
            ax.relim()
            ax.autoscale_view()
            ax.grid()

    anim = animation.FuncAnimation(fig, animate,
                                   interval=1,
                                   repeat=False,
                                   frames=len(spec_fnames)-1)
    plt.show()

    results_df.to_csv('fit_results.csv')

#==============================================================================
#================================== Scanning ==================================
#==============================================================================

elif mode == 'scan':

    # Read in the scan block
    err, x, info, spec_block = read_scan(fpath)

    # Read in the dark
    common['dark'] = spec_block[0]

    # Find the fit and stray windows
    common['fit_idx'] = np.where(np.logical_and(x > settings['w_lo'],
                                                x < settings['w_hi']))

    common['stray_idx'] = np.where(np.logical_and(x > 280,  x < 290))

    spec_num = np.arange(len(spec_block[1:]))
    so2 = np.full(len(spec_block[1:]), np.nan)
    err = np.full(len(spec_block[1:]), np.nan)

    # Set up figure
    fig = plt.figure(figsize = [8,8])
    gs = gridspec.GridSpec(3,1)
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])
    ax2 = fig.add_subplot(gs[2])

    axes = [ax0, ax1, ax2]

    def animate(i):

        # Unpack the spectrum
        y = spec_block[i+1]

        # Pre-process the spectrum before the fit
        spectrum = pre_process([x, y], common)
        grid, spec = spectrum

        # Fit the spectrum
        fit_result = fit_spectrum(spectrum, common, update_params = True)
        fit_result.print_result()

        so2[i] = fit_result.params['SO2_295'].fit_val/2.463e15
        err[i] = fit_result.params['SO2_295'].fit_err/2.463e15

        for ax in axes:
            ax.clear()

        ax0.plot(grid, spec, 'C0o', label = 'Data')
        ax0.plot(grid, fit_result.fit, 'C1-', lw = 2, label = 'Fit')
        ax1.plot(grid, fit_result.resid, 'C0o-')
        ax2.errorbar(spec_num, so2, yerr = err, color='g', marker='o', capsize=2)

        ax0.legend(loc = 'upper left')
        ax0.set_ylabel('Instensity (counts)', fontsize = 14)

        ax1.set_xlabel('Wavelength (nm)', fontsize = 14)
        ax1.set_ylabel('Residual (%)', fontsize = 14)
        ax1.set_ylim(-3, 3)

        ax2.set_xlabel('Spectrum Number', fontsize = 14)
        ax2.set_ylabel('SO$_2$ SCD (ppm.m)', fontsize = 14)
        ax2.set_xlim(-10, 112)
        #ax2.set_ylim(-100, 400)

        plt.tight_layout()

        # Rescale the axes
        for ax in axes:
            ax.relim()
            ax.autoscale_view()
            ax.grid()

    anim = animation.FuncAnimation(fig, animate,
                                   interval=1,
                                   repeat=False,
                                   frames=len(spec_block)-1)
    plt.show()


