#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 09:00:20 2019

@author: mqbpwbe2
"""

import numpy as np
from math import sin, cos, atan2

#========================================================================================
#====================================== haversine =======================================
#========================================================================================

def haversine(start_coords, end_coords, radius = 6371000):
    
    '''
    Function to calculate the distance and initial bearing between two points
    
    INPUTS
    ------
    start_coords, tuple
        Start coordinates (lat, lon) in decimal degrees (+ve = north/east)
           
    end_coords, tuple
        End coordinates (lat, lon) in decimal degrees (+ve = north/east)
        
    radius, float
        Radius of the body in meters. Default is set to the Earth radius (6731km)
        
    OUTPUTS
    -------
    distance, float
        The linear distance between the two points in meters
        
    bearing, float
        The initial bearing between the two points (radians)
    '''
    
    # Unpack the coordinates and convert to radians
    lat1, lon1 = np.radians(start_coords)
    lat2, lon2 = np.radians(end_coords)
    
    # Calculate the change in lat and lon
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    # Calculate the square of the half chord length
    a = (sin(dlat/2))**2 + ( cos(lat1) * cos(lat2) * (sin(dlon/2))**2 )
    
    # Calculate the angular distance
    c = 2 * atan2(np.sqrt(a), np.sqrt(1-a))
    
    # Find distance moved
    distance = radius * c

    # Calculate the initial bearing
    bearing = atan2(sin(dlon) * cos(lat2),
                    (cos(lat1) * sin(lat2)) - (sin(lat1) * cos(lat2) * cos(dlon)))
    
    return distance, bearing
    
if __name__ == '__main__':
    dist0, bear0 = haversine([37.7505, 14.9934], [37.7905, 15.1386])
    
    print(dist0, np.degrees(bear0))