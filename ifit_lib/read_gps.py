from math import radians, cos, sin, asin, atan2, sqrt, pi
import numpy as np

#========================================================================================
#======================================read_txt_gps======================================
#========================================================================================

# Function to read a GPS text file and return the time and coordinates of the readings

# INPUTS:  gps_fname; filepath to the gps text file

# OUTPUTS: date, time, latitude, longitude and altitude of readings

def read_txt_gps(gps_fname, time_diff):

    # Create empty arrays to hold the outputs
    date = []
    time = np.array(())
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
            date.append(date_text)
            lat = np.append(lat,info[2])
            lon = np.append(lon,info[3])
            alt = np.append(alt,info[4])
            
            # Extract time and convert to julian time
            h = float(time_text[0:2])
            m = float(time_text[3:5])
            s = float(time_text[6:8])
            
            t = (h * 3600.0 + m * 60.0 + s) / 86400.0
            
            # Correct for time difference
            t = t - time_diff / 24
            
            if t > 1:
                t = t - 1
            
            # Append to array            
            time = np.append(time,t)
            
    return date, time, lat, lon, alt

    
#========================================================================================
#=======================================gps_vector=======================================
#========================================================================================



def haversine(lon1, lat1, lon2, lat2):
    '''
    Function to calculate the displacement and bearing between two GPS corrdinates
    
    INPUTS
    ------
    lon1, lat1: longitude and latitude of first point
    lon2, lat2: longitude and latitude of second point
    
    OUTPUTS
    -------
    dist: distance between two points in meters
    bearing: bearing between points (0 - 2pi clockwise from North)
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