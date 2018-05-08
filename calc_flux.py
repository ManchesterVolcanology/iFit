import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
from math import cos, pi
import datetime
from scipy.interpolate import griddata
from pandas import read_csv

from ifit_lib.click_zoom import click_zoom
from ifit_lib.read_gps import read_txt_gps, gps_vector, haversine, read_nmea
from ifit_lib.center_of_grav import cog

#========================================================================================
#====================================Program controls====================================
#========================================================================================

# Define the location of the volcano
# Soufriere Hills:
#volc_lon = -62.175384
#volc_lat = 16.710445
#time_diff = 4

# Etna
#volc_lon = 14.9929
#volc_lat = 37.7509
#time_diff = -1

# Papandayan
#volc_lon = 107.7310
#volc_lat = -7.3193
#time_diff = -7

# Bromo
#volc_lon = 112.942708
#volc_lat = -7.942708
#time_diff = -7

# Arjuno
#volc_lon = 112.58944
#volc_lat = -7.76444
#time_diff = -7

# Semeru
#volc_lon = 112.9224
#volc_lat = -8.1077
#time_diff = -7

# Merapi
#volc_lon = 110.4457
#volc_lat = -7.5407
#time_diff = -7

# Ijen
#volc_lon = 114.242
#volc_lat = -8.058
#time_diff = -7

# Masaya
#volc_lon = -86.169058
#volc_lat = 11.984104
#time_diff = 6

# Piton de la Fournaise
volc_lon = 55.708889
volc_lat = -21.2425
time_diff = 0

# Create output folder
out_folder = 'Results/calc_flux/' + str(datetime.date.today()) + '/'
if not os.path.exists(out_folder):
    os.makedirs(out_folder)
    
# Create output .csv file if it doesn't exist
out_fname = out_folder + 'flux_out.csv'

with open(out_fname, 'w') as writer:
    writer.write('calc_flux.py output file \n')
    
# Set font size for graphics
plt.rcParams.update({'font.size': 12})

#========================================================================================
#=================================Read in Volcano info===================================
#======================================================================================== 

# Read in volcano location and time difference
with open('Flux Inputs/volcano_data.txt', 'r') as r:
    
    # Read file
    data = r.readlines()
               
    # Unpack data
    volc_lon  = float(data[0].split(';')[1])
    volc_lat  = float(data[1].split(';')[1])
    time_diff = float(data[2].split(';')[1])   

#========================================================================================
#===================================Read in GPS file=====================================
#========================================================================================
    
# Read GPS data from file
print('Reading GPS data...')
#gps_time, lat, lon = read_nmea('Flux inputs/180430d.TXT', time_diff)
gps_time, lat, lon = read_txt_gps('Flux inputs/gps.txt', time_diff)

print('Done!')

#========================================================================================
#==================================Read in data file=====================================
#========================================================================================

print('Reading traverse data...')

# Read data file
try:
    data = read_csv('Flux inputs/ifrit_out.csv')
    
except FileNotFoundError:
    data = read_csv('Flux inputs/iFit_out.csv')

# Load required data into arrays
txt_time    = data['Time']
txt_so2_amt = data['so2 (ppm.m)']
txt_so2_err = data['so2 error']

# Convert the time to julian time
time = []

for i in txt_time:
    h = float(i[0:2])
    m = float(i[3:5])
    s = float(i[6:8])
    time.append((h * 3600.0 + m * 60.0 + s) / 86400.0)
    
# Convert others from string to number arrays
so2_amt = []
for i in txt_so2_amt:
    so2_amt.append(float(i))
    
so2_err = []
for i in txt_so2_err:
    so2_err.append(float(i))

print('Done!')

#========================================================================================
#====================================Set wind speed======================================
#========================================================================================

# Request wind speed from the user
while True:
    try:
        wind_speed = float(input('Input wind speed (ms-1): '))
        input_flag = True
        break
    
    except ValueError:
        print('Input not a float, please re-enter')
        
# Request wind speed uncertainty from the user
while True:
    try:
        wind_err = float(input('Input wind speed uncertainty (%): ')) / 100
        input_flag = True
        break
    
    except ValueError:
        print('Input not a float, please re-enter')
       
#========================================================================================
#==================================Subtract background===================================
#========================================================================================
'''  
bg_flag = input('Subtract background? y/n: ')

if bg_flag == 'y':
    #Plot SO2 time series
    fig = plt.figure(1, figsize = (12,8))
    ax = fig.add_subplot(111)
    ax.plot(time, so2_amt, 'g', lw=1.5)
    ax.fill_between(time, np.subtract(so2_amt, so2_err), np.add(so2_amt, so2_err), 
                    color = 'lightgreen')
    ax.set_xlabel('Julian Time (Fraction of Day)')
    ax.set_ylabel('SO2 column amount (ppm.m)')
    
    plt.show()
    
    # User inputs the background
    background = float(input('Background level: '))
    
    # Apply correction
    so2_amt = np.subtract(so2_amt, background)
'''
#========================================================================================
#==================================Begin analysis loop===================================
#========================================================================================

# Create list to hold all so2 fluxes
fluxes = []

loop_flag = 'y'
loop = 1

# Keep a record of origional values
time_old = time
so2_old = so2_amt
err_old = so2_err

while loop_flag == 'y':
    
    # Re-set variables
    time = time_old
    so2_amt = so2_old
    so2_err = err_old
    
#========================================================================================
#================================Select flux time window=================================
#========================================================================================

    # Apply two sets of windowing, one rough to get the right region, then one 
    #  accurate to select the plume
    for i in range(2):
        
        # Plot the SO2 time series data
        fig = plt.figure(1, figsize = (12,8))
        ax = fig.add_subplot(111)
        ax.plot(time, so2_amt, 'g', lw=1.5)
        ax.fill_between(time, np.subtract(so2_amt, so2_err), np.add(so2_amt, so2_err), 
                        color = 'lightgreen')
        ax.set_xlabel('Julian Time (Fraction of Day)')
        ax.set_ylabel('SO2 column amount (ppm.m)')
        ax.set_title('Select window')
        
        # Cut to right region of data
        ind1, ind2 = click_zoom(fig, time)
        
        # Cut the so2 amount timeseries to selected window
        time    = time[ind1:ind2]
        so2_amt = so2_amt[ind1:ind2]
        so2_err = so2_err[ind1:ind2]

#========================================================================================
#============================Define peak of SO2 concentration============================
#========================================================================================

    peak_index = cog(so2_amt)

#========================================================================================
#========================Interpolate GPS data onto SO2 time grid=========================
#========================================================================================

    # Interpolate the GPS data onto the so2 time grid
    modlon_old = griddata(gps_time,lon,time_old)
    modlat_old = griddata(gps_time,lat,time_old)
    modlat = griddata(gps_time,lat,time)
    modlon = griddata(gps_time,lon,time)
    
#========================================================================================
#==============================Perform geometric correction==============================
#========================================================================================

    # Define wind vector as the vector from the volcano to the CoM
    wind_vector = haversine(volc_lon, volc_lat, modlon[peak_index], modlat[peak_index])
    
    # Find the distance-bearing vectors between each reading. Bearing is radians 
    #  anticlockwise from East from 0 to 2pi
    displacement, bearing, dir_corr = gps_vector(modlon, modlat, wind_vector[1])

    # Correction factor to account for the non-orthogonality of the plume
    orthog_correction = np.array(())
    for i in bearing:
        correction = abs(cos((pi / 2) - (wind_vector[1] - i)))
        orthog_correction = np.append(orthog_correction, correction)
        
    # Correct for if going the wrong way
    orthog_correction = np.multiply(orthog_correction, dir_corr)    
    
    # Apply the correction
    corr_dists = np.multiply(displacement, orthog_correction)    

    # Calculate cumulative distance
    cuml_dist = np.zeros(1)
    for i in range(len(corr_dists)):
        cuml_dist = np.append(cuml_dist,cuml_dist[i] + corr_dists[i])
    
#========================================================================================
#=============================Plot selected SO2 and GPS data=============================
#========================================================================================
        
    # Create plot grid
    fig = plt.figure(figsize = (12,8))
    gs  = gridspec.GridSpec(2, 3, width_ratios = [1,0.95,0.05])
    ax1 = plt.subplot(gs[0,0])  
    ax2 = plt.subplot(gs[0,1])
    ax3 = plt.subplot(gs[1,0])
    ax4 = plt.subplot(gs[1,1])   

    # Plot the full SO2 amount time series
    ax1.plot(time, so2_amt,  'g', lw = 1.5, label = 'SO2 time series')
    ax1.fill_between(time, np.subtract(so2_amt, so2_err), np.add(so2_amt, so2_err), 
                     color = 'lightgreen')
    ax1.plot(time[peak_index], so2_amt[peak_index], 'ko', label = 'Plume CoM')
    ax1.set_xlabel('Julian Time (Fraction of Day)')
    ax1.set_ylabel('SO2 column amount (ppm.m)')
    ax1.legend(loc = 0, fontsize = 10)

    # Plot the selected traverse points and SO2 amounts wrt the volcano
    sc = ax2.scatter(modlon, modlat, c = so2_amt, s=50, cmap = 'plasma', lw=0.0,
                     alpha = 0.5)
    cbax = plt.subplot(gs[0,2])
    cb = Colorbar(ax = cbax, mappable = sc)
    cb.set_label('SO2 column amount (ppm.m)')
    
    ax2.scatter(volc_lon, volc_lat, c='darkorange', s=180, label = 'Volcano')
    
    ax2.set_xlabel('Longitude (deg)')
    ax2.set_ylabel('Latitude (deg)')
    ax2.legend(loc = 0, fontsize = 10)

    # Plot the full SO2 amount time series overlaid with the selected section
    ax3.plot(time_old, so2_old, 'g-', label = 'Full so2 time series')
    ax3.fill_between(time_old, np.subtract(so2_old, err_old), np.add(so2_old, err_old), 
                     color = 'lightgreen')
    ax3.plot(time, so2_amt, 'r-', label = 'Traverse')
    ax3.fill_between(time, np.subtract(so2_amt, so2_err), np.add(so2_amt, so2_err), 
                     color = 'coral')
    ax3.set_xlabel('Julian Time (Fraction of Day)')
    ax3.set_ylabel('SO2 column amount (ppm.m)')
    ax3.legend(loc = 0, fontsize = 10)
    
    # Plot the full GPS track and overlay the selected window wrt the volcano location
    ax4.plot(modlon_old, modlat_old, 'g-', label = 'Full Track')
    ax4.plot(modlon, modlat, 'r', lw = 1.5, label = 'Traverse')
    ax4.scatter(modlon[0], modlat[0], c='g', s=100, label = 'Start')
    ax4.scatter(modlon[-1], modlat[-1], c='r', s=100, label = 'Stop')
    ax4.scatter(volc_lon, volc_lat, c='darkorange', s=200, label = 'Volcano')
    ax4.scatter(modlon[peak_index], modlat[peak_index], c='k', s=100)
    ax4.set_xlabel('Longitude (deg)')
    ax4.set_ylabel('Latitude (deg)')
    ax4.legend(loc = 0, fontsize = 10)

    # Display the plots
    plt.tight_layout()
    plt.show()

#========================================================================================
#=====================================Calculate flux=====================================
#========================================================================================

    # Indentify only unique distance values to ignore times standing still
    #unique_cuml_dist, ind = np.unique(cuml_dist,return_index=True)
    #unique_so2_amt = np.take(so2_amt, ind)
    
    # Convert so2 amounts from ppm.m to molecules.cm-2
    so2_amt_molec = np.multiply(so2_amt, 2.463e15)

    # Convert so2 amounts from molecules.cm-2 to molecules.m-2
    molec_m2 = so2_amt_molec[1:] * 1.0e4

    # Convert to molecules.m-1 and sum over the traverse
    molec_m = np.sum(np.multiply(molec_m2, corr_dists))

    # Multiply by windspeed
    molec_per_s = molec_m * wind_speed

    # Convert to moles.s-1
    moles_per_s = molec_per_s / 6.022e23

    # Convert to kg.s-1. Molar mass of so2 is 64.066g/mole
    kg_per_s = moles_per_s * 0.064066

    # Convert to tonnes per day, make sure it's positive
    flux = abs(int((kg_per_s * 86400.0) / (1000.0)))
    
#========================================================================================
#==================================Calculate flux error==================================
#========================================================================================
    
    # Calculate error in so2 amount
    total_err = np.sum(so2_err)
    total_so2 = np.sum(so2_amt)
    
    delta_A = total_err / total_so2
    
    # Combine with wind speed error
    delta_F = ( (delta_A)**2 + (wind_err)**2 )**0.5
    
    flux_err = int(flux * delta_F)

    print('Flux = ' + str(flux) + ' (+/- ' + str(flux_err) + ') tonnes per day of SO2')
    
#========================================================================================
#======================================Output data=======================================
#========================================================================================

    # Make strings of volcano and COG positions
    volc_pos = str(volc_lon)+ ',' +str(volc_lat)
    cog_pos = str(modlon[peak_index])+ ',' +str(modlat[peak_index])
    
    keep_flag = input('Keep data? (y/n): ')
        
    if keep_flag == 'y':
        
        # Append flux and error to list
        fluxes.append([flux, flux_err]) 

        with open(out_fname, 'a') as writer:
            # Write header information to the file
            writer.write('Peak number,' + str(loop)           + '\n' \
                         'Flux (Tonnes/day),' + str(flux)     + '\n' \
                         'Windspeed (m/s),' + str(wind_speed) + '\n' \
                         'Volcano Long/Lat,' + volc_pos       + '\n' \
                         'CoM position,' + cog_pos            + '\n' )
             
            # Write column headers            
            writer.write('Julian Time (GMT),SO2 Amount,Error,Longitude,Latitude,' + \
                         'Bearing,Distance,Correction Factor,Cumulative Distance\n')
            
            # Write information to the file
            for i in range(1,len(modlon)):
                writer.write(str(time[i])                + ',' + \
                             str(so2_amt[i])             + ',' + \
                             str(so2_err[i])             + ',' + \
                             str(modlon[i])              + ',' + \
                             str(modlat[i])              + ',' + \
                             str(bearing[i-1])           + ',' + \
                             str(displacement[i-1])      + ',' + \
                             str(orthog_correction[i-1]) + ',' + \
                             str(cuml_dist[i])           + '\n')

            writer.write('\n')
            
            # Replot the data and save the figure
            
            # Create plot grid
            fig = plt.figure(figsize = (12,8))
            gs  = gridspec.GridSpec(2, 3, width_ratios = [1,0.95,0.05])
            ax1 = plt.subplot(gs[0,0])  
            ax2 = plt.subplot(gs[0,1])
            ax3 = plt.subplot(gs[1,0])
            ax4 = plt.subplot(gs[1,1])   
        
            # Plot the full SO2 amount time series
            ax1.plot(time, so2_amt,  'g', lw = 1.5, label = 'SO2 time series')
            ax1.fill_between(time, np.subtract(so2_amt, so2_err), np.add(so2_amt, so2_err), 
                             color = 'lightgreen')
            ax1.plot(time[peak_index], so2_amt[peak_index], 'ko', label = 'Plume CoM')
            ax1.set_xlabel('Julian Time (Fraction of Day)')
            ax1.set_ylabel('SO2 column amount (ppm.m)')
            ax1.legend(loc = 0, fontsize = 10)
        
            # Plot the selected traverse points and SO2 amounts wrt the volcano
            sc = ax2.scatter(modlon, modlat, c = so2_amt, s=50, cmap = 'plasma', lw=0.0,
                             alpha = 0.5)
            cbax = plt.subplot(gs[0,2])
            cb = Colorbar(ax = cbax, mappable = sc)
            cb.set_label('SO2 column amount (ppm.m)')
            
            ax2.scatter(volc_lon, volc_lat, c='darkorange', s=180, label = 'Volcano')
            
            ax2.set_xlabel('Longitude (deg)')
            ax2.set_ylabel('Latitude (deg)')
            ax2.legend(loc = 0, fontsize = 10)
        
            # Plot the full SO2 amount time series overlaid with the selected section
            ax3.plot(time_old, so2_old, 'g-', label = 'Full so2 time series')
            ax3.fill_between(time_old, np.subtract(so2_old, err_old), np.add(so2_old, err_old), 
                             color = 'lightgreen')
            ax3.plot(time, so2_amt, 'r-', label = 'Traverse')
            ax3.fill_between(time, np.subtract(so2_amt, so2_err), np.add(so2_amt, so2_err), 
                             color = 'coral')
            ax3.set_xlabel('Julian Time (Fraction of Day)')
            ax3.set_ylabel('SO2 column amount (ppm.m)')
            ax3.legend(loc = 0, fontsize = 10)
            
            # Plot the full GPS track and overlay the selected window wrt the volcano location
            ax4.plot(modlon_old, modlat_old, 'g-', label = 'Full Track')
            ax4.plot(modlon, modlat, 'r', lw = 1.5, label = 'Traverse')
            ax4.scatter(modlon[0], modlat[0], c='g', s=100, label = 'Start')
            ax4.scatter(modlon[-1], modlat[-1], c='r', s=100, label = 'Stop')
            ax4.scatter(volc_lon, volc_lat, c='darkorange', s=200, label = 'Volcano')
            ax4.scatter(modlon[peak_index], modlat[peak_index], c='k', s=100)
            ax4.set_xlabel('Longitude (deg)')
            ax4.set_ylabel('Latitude (deg)')
            ax4.legend(loc = 0, fontsize = 10)
            
            # Save the figure
            plt.tight_layout()
            fig.savefig(out_folder + 'traverse' + str(loop) + '.png')
            plt.clf()
            
            '''
            # Plot all the traverses saved so far
            plt.scatter(npy_lon, npy_lat, c = npy_so2, s = 50, alpha = 0.5, cmap = 'plasma')
            plt.scatter(volc_lon, volc_lat, c = 'darkorange', s = 200)
            plt.xlabel('Longitude')
            plt.ylabel('Latitude')
            plt.show()
            '''
            
            loop += 1
    else:
        print('Data discarded')
            
    loop_flag = input('Analyse another traverse? (y/n): ')

# Calculate the average flux and uncertainty
av_flux = 0
av_err  = 0
n = 0 

for i in fluxes:
    av_flux += i[0]
    av_err  += (i[1]/i[0])**2
    n += 1
    
av_flux = av_flux / n
av_err  = (av_err / n)**0.5 * av_flux

# Save the fluxes in a text file for ease of access
header = 'Results from calc_flux.py\n' + \
         'NOTE errors are from SO2 fitting and wind speed only\n' + \
         'Flux (tonnes/day),     Error (+/- t/d)'

np.savetxt(out_folder + 'flux_results.txt', fluxes, header = header)

with open(out_folder + 'flux_results.txt', 'a') as r:
    r.write('Average flux = ' + str(int(av_flux)) + ' +/- ' + str(int(av_err)) + \
            ' (t/day)')

'''

# Ask if the user wants togenerate a map of the traverses
map_flag = input('Generate map of traverses? (y/n): ')

if map_flag == 'y':

    # Plot traverses on a map
    print('Generating map...')
    # Set map coordinates from volcano location
    llcrnrlon = volc_lon - 0.08
    llcrnrlat = volc_lat - 0.07
    urcrnrlon = volc_lon + 0.08
    urcrnrlat = volc_lat + 0.12
    
    m = Basemap(projection='merc', resolution = 'h', area_thresh = 0.1,
           llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
           urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
    
    #lon, lat, so2 = np.loadtxt(out_folder + 'trav_data.txt', unpack = True)
    
    m.drawcoastlines()
    m.fillcontinents (color='lightgreen', zorder = 1)
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='blue',zorder=0)
    
    m.drawcounties()
    
    
    x, y = m(modlon, modlat)
    volcx, volcy = m(volc_lon, volc_lat)
    sc = m.scatter(x, y, c=so2_amt, s=50, alpha=0.5, cmap='plasma', zorder=2)
    m.scatter(volcx, volcy, c = 'r', zorder = 2)
    
    cb = m.colorbar(sc, label = 'SO2 path amount (ppm.m)')           
    
    plt.savefig(out_folder + 'traverse_map.png') 
    
    plt.show()
'''