from math import radians, cos, sin, asin, atan2, sqrt, pi
import numpy as np
import pynmea2

from ifit_lib.julian_time import hms_to_julian

#========================================================================================
#======================================= read_gps =======================================
#========================================================================================

def read_gps(fpath, datatype = 'text'):
    
    '''
    Function to read in GPS data
    
    INPUTS
    ------
    fpath, string
        File path to the GPS file
        
    datatype, string (optional)
        Style of GPS data. Must be one of 'text' or 'NMEA' (default is 'text')
    '''
    
    # Load the file
    with open (fpath, 'r') as r:
        
        # Read in the data
        lines = r.readlines()
    
        # Create empty arrays to hold the data
        time = []
        lat = []
        lon = []
        alt = []
        
        if datatype == 'NMEA':
        
            for n, line in enumerate(lines):
                
                if line[1:6] == 'GPGGA':
                    
                    try:
                        msg = pynmea2.parse(line)
                        
                        time.append(msg.timestamp)
                        lat.append(msg.latitude)
                        lon.append(msg.longitude)
                        alt.append(msg.altitude)
            
                    except AttributeError:
                        pass
                    
        if datatype == 'text':
            
            for n, line in enumerate(lines[1:]):
                
                # Split data using tab delimiter
                info = line.split('\t')
                
                # Split date and time with space delimiter
                date_text, time_text = info[1].split(' ')
                
                # Append values to arrays
                time.append(time_text)
                lat.append(float(info[2]))
                lon.append(float(info[3]))
                alt.append(float(info[4]))
                
        
    # Convert time to julian 
    try:           
        time = hms_to_julian(time, str_format='%H:%M:%S', out_format='decimal hours')
    except ValueError:
        time = hms_to_julian(time, str_format='%H:%M:%S.%f', out_format='decimal hours')        
                
    return time, lat, lon, alt


#========================================================================================
#===================================== read_txt_gps =====================================
#========================================================================================

def read_txt_gps(gps_fname):

    '''
    Function to read in gps .txt file 
    (e.g. converted from .nmea by http://www.gpsvisualizer.com)
    
    INPUTS
    ------
    gps_fname, str
        File path to the gps file to load
    
    OUTPUTS
    -------
    time, array
        Time values in decimal hours
        
    lat, array
        Point latitudes
        
    lon, array
        Point longitudes
    '''
    
    # Create empty arrays to hold the outputs
    time = []
    lat  = np.array(())
    lon  = np.array(())
    alt  = np.array(())
    
    # Counts number of lines in the text file
    n_lines = sum(1 for line in open(gps_fname))
    
    # Iterate over the lines, stroing the needed data in arrays
    with open(gps_fname, 'r') as reader:
        
        # Read first title line
        reader.readline()
        
        for i in range(1,n_lines):
            
            # Read the line
            line = reader.readline()
            
            # Split data using tab delimiter
            info = line.split('\t')
            
            # Split date and time with space delimiter
            date_text, time_text = info[1].split(' ')
            
            # Append values to arrays
            lat = np.append(lat, float(info[2]))
            lon = np.append(lon, float(info[3]))
            alt = np.append(alt, float(info[4]))
            time.append(time_text)
            
    # Convert time to julian 
    try:           
        time = hms_to_julian(time, str_format = '%H:%M:%S', out_format = 'decimal hours')
    except ValueError:
        time = hms_to_julian(time, str_format = '%H:%M:%S.%f',out_format='decimal hours')
            
    return time, lat, lon, alt

    
#========================================================================================
#====================================== gps_vector ======================================
#========================================================================================



def haversine(lon1, lat1, lon2, lat2):
    
    '''
    Function to calculate the displacement and bearing between two GPS corrdinates
    
    INPUTS
    ------
    lon1, lat1, floats
        Longitude and latitude of first point
        
    lon2, lat2, floats
        Longitude and latitude of second point
    
    OUTPUTS
    -------
    dist, float
        Distance between two points in meters
        
    bearing, float
        Bearing between points (0 - 2pi clockwise from North)
    '''

    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1    
    
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    
    # Earth radius in meters
    r = 6371000
    
    # Calculate distance
    dist = c * r
    
    # Calculate the bearing
    bearing = atan2(sin(dlon) * cos(lat2), 
                    cos(lat1) * (sin(lat2) - sin(lat1)) * cos(lat2) * cos(dlat))
                    
    # Convert barings to the range (0, 2pi) instead of (-pi, pi)
    if bearing < 0:
        bearing = 2 * pi + bearing
    
    return dist, bearing
    

def gps_vector(lon_arr, lat_arr, wind_bearing):

    # Create zero arrays to hold displacement-bearing vectors
    dist = np.zeros(len(lon_arr) - 1)
    bearing = np.zeros(len(lon_arr) - 1)
    
    # Loop over the input arrays (leaving the last as it requires a difference)
    for i in range(len(lon_arr) - 1):
        
        # Extract the lat/long values from the arrays
        lon1, lat1 = lon_arr[i], lat_arr[i]
        lon2, lat2 = lon_arr[i+1], lat_arr[i+1]
        
        # Input into the haversine function
        dist[i], bearing[i] = haversine(lon1, lat1, lon2, lat2)
        
    # Find relative bearing to wind vector
    rel_bearing = np.subtract(bearing, wind_bearing)
    
    # Convert negative relative bearings
    for i in np.where(rel_bearing < 0):
        rel_bearing[i] = rel_bearing[i] + (2 * pi)
        
    # Create array of modifications to correct for inflections in the traverse path
    dir_corr = np.ones(len(rel_bearing))
    for i in np.where(rel_bearing > pi):
        dir_corr[i] = -1
    
    return dist, bearing, dir_corr