"""Contains functions for geometrical calculations on a globe."""
import numpy as np
from math import sin, cos, atan2, asin, pi


# =============================================================================
# haversine
# =============================================================================

def haversine(start_coords, end_coords, radius=6371000):
    """Calculate the distance and initial bearing between two points.

    Parameters
    ----------
    start_coords : tuple
        Start coordinates (lat, lon) in decimal degrees (+ve = north/east)

    end_coords : tuple
        End coordinates (lat, lon) in decimal degrees (+ve = north/east)

    radius : float
        Radius of the body in meters. Default is set to the Earth radius
        (6731km)

    Returns
    -------
    distance : float
        The linear distance between the two points in meters

    bearing : float
        The initial bearing between the two points (radians)
    """
    # Unpack the coordinates and convert to radians
    lat1, lon1 = np.radians(start_coords)
    lat2, lon2 = np.radians(end_coords)

    # Calculate the change in lat and lon
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Calculate the square of the half chord length
    a = (sin(dlat/2))**2 + (cos(lat1) * cos(lat2) * (sin(dlon/2))**2)

    # Calculate the angular distance
    c = 2 * atan2(np.sqrt(a), np.sqrt(1-a))

    # Find distance moved
    distance = radius * c

    # Calculate the initial bearing
    bearing = atan2(sin(dlon) * cos(lat2),
                    (cos(lat1)*sin(lat2)) - (sin(lat1)*cos(lat2)*cos(dlon)))

    # Check the bearing is between 0 and 2pi
    bearing = bearing_check(bearing)

    return distance, bearing


# =============================================================================
# calc_end_point
# =============================================================================

def calc_end_point(start_coords, dist, bearing, radius=6371000):
    """Calculate the final coordinates given a starting position and vector.

    Parameters
    ----------
    start_coords : tuple
        Starting coordinates (lat, lon) in decimal degrees (+ve = north/east)

    dist : float
        The distance moved in meters

    bearing : float
        The bearing of travel in degrees clockwise from north

    radius : float
        Radius of the body in meters. Default is set to the Earth radius
        (6731km)

    Returns
    -------
    end_coords, tuple
        The final coordinates (lat, lon) in decimal degrees (+ve = north/east)
    """
    # Convert the inputs to radians
    lat, lon = np.radians(start_coords)
    theta = np.radians(bearing)

    # Calculate the angular distance moved
    ang_dist = dist / radius

    # Calculate the final latitude
    end_lat = asin(np.add((sin(lat) * cos(ang_dist)),
                          (cos(lat) * sin(ang_dist) * cos(theta))))

    # Calculate the final longitude
    end_lon = lon + atan2(sin(theta) * sin(ang_dist) * cos(lat),
                          cos(ang_dist) - (sin(lat)*sin(end_lat)))

    return np.degrees([end_lat, end_lon])


# =============================================================================
# Add a bearing checker
# =============================================================================

def bearing_check(bearing, radians=True, max_iter=1000):
    """Check for a valid bearing (between 0 -> 360 degrees).

    Parameters
    ----------
    bearing : float
        The bearing value to check.
    radians : bool, optional
        If True then the bearing is treated as radians. Otherwise it is assumed
        to be degrees. Default is True
    max_iter : int
        Maximum number of checks to avoid infinite loops. Default is 1000

    Returns
    -------
    float
        The checked bearing.
    """
    # Check if degrees
    if not radians:
        bearing = np.radians(bearing)

    i = 0

    while bearing < 0 or bearing >= 2*pi:
        if bearing < 0:
            bearing += 2*pi
            continue
        if bearing >= 2*pi:
            bearing -= 2*pi
            continue
        i += 1

        if i >= max_iter:
            msg = f'Max iteration in bearing check reached {max_iter}'
            raise ValueError(msg)

    # Check if degrees
    if not radians:
        bearing = np.degrees(bearing)

    return bearing


if __name__ == '__main__':
    dist, bearing = haversine([37.7505, 14.9934], [37.7905, 15.1386])

    print(dist, np.degrees(bearing))

    lat, lon = calc_end_point([37.7505, 14.9934], 13520, 70.74)

    print(lat, lon)
