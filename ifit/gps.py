"""Provides an interface to a USB GPS device."""
import utm
import time
import serial
import logging
import numpy as np
from threading import Thread
from datetime import datetime
import serial.tools.list_ports


logger = logging.getLogger(__name__)


class GPS():
    """GPS object."""

    def __init__(self, comport=None, filename=None, baudrate=4800,
                 parity=serial.PARITY_NONE, stopbits=serial.STOPBITS_ONE,
                 bytesize=serial.EIGHTBITS):
        """Initialize."""
        if comport is None:
            comport = serial.tools.list_ports.comports()[0].device
        self.serial_port = serial.Serial(comport, baudrate=baudrate,
                                         parity=parity, stopbits=stopbits,
                                         bytesize=bytesize)

        self.filename = filename
        self.timestamp = None
        self.datestamp = None
        self.lat = np.nan
        self.lon = np.nan
        self.alt = np.nan
        self.utm_coords = None
        self.running = True
        self.acquired = False

        self.thread = Thread(target=self._updater, daemon=True)
        self.thread.start()

    def _updater(self):
        while self.running:
            try:
                ser_bytes = self.serial_port.readline()
                decoded_bytes = ser_bytes.decode('utf-8')
                if self.filename is not None and self.filename != '':
                    try:
                        with open(self.filename, 'a') as w:
                            w.write(decoded_bytes.strip() + '\n')
                    except FileNotFoundError:
                        logger.warning(f'Unable to find file {self.filename}'
                                       + ' Disabling GPS file stream.')
                        self.filename = None

                # Exctract location information
                data = decoded_bytes.split(",")

                if 'GGA' in data[0]:
                    self._parse_gpgga(data)

                if 'RMC' in data[0]:
                    self._parse_gprmc(data)

            except UnicodeDecodeError:
                time.sleep(1)

            except serial.SerialException:
                logger.warning('GPS disconnected!')

    def _parse_gpgga(self, data):
        try:
            # Read timestamp
            self.timestamp = datetime.strptime(data[1], '%H%M%S.%f').time()

            # Read lat/lon info
            lat_str = data[2]
            lat_dir = data[3]
            lon_str = data[4]
            lon_dir = data[5]

            # Convert to decimel degrees
            if lat_str != '':
                lat = float(lat_str[:2]) + float(lat_str[2:])/60
                if lat_dir == 'S':
                    lat = -lat
                self.lat = lat
            if lon_str != '':
                lon = float(lon_str[:3]) + float(lon_str[3:])/60
                if lon_dir == 'W':
                    lon = -lon
                self.lon = lon

            # Unpack altitude
            alt_str = data[9]
            if alt_str != '':
                alt = float(alt_str)
                alt_unit = data[10]

                # Convert from feet to meters if required
                if alt_unit == 'F':
                    alt = 0.3048 * alt
                self.alt = alt

            # Convert to UTM units
            try:
                self.utm_coords = utm.from_latlon(lat, lon)
            except UnboundLocalError:
                pass

            self.acquired = True

        except ValueError as e:
            logger.debug(f'Error parsing GPS string\n{e}')
            pass

    def _parse_gprmc(self, data):
        """Parse NMEA GPRMC string."""
        try:
            # Read timestamp
            self.timestamp = datetime.strptime(data[1], '%H%M%S.%f').time()

            # Read date stamp
            self.datestamp = datetime.strptime(data[9], '%d%m%y').date()

            # Read lat/lon info
            lat_str = data[3]
            lat_dir = data[4]
            lon_str = data[5]
            lon_dir = data[6]

            # Convert to decimel degrees
            if lat_str != '':
                lat = float(lat_str[:2]) + float(lat_str[2:])/60
                if lat_dir == 'S':
                    lat = -lat
                self.lat = lat
            if lon_str != '':
                lon = float(lon_str[:3]) + float(lon_str[3:])/60
                if lon_dir == 'W':
                    lon = -lon
                self.lon = lon

        except ValueError as e:
            logger.debug(f'Error parsing GPS string\n{e}')
            pass

    def get_position(self,  time_to_wait=60):
        """Report current time and position."""
        t0 = datetime.now()
        while (datetime.now() - t0).total_seconds() < time_to_wait:

            flags = [~np.isnan(self.lat), ~np.isnan(self.lon),
                     self.timestamp is not None, self.datestamp is not None]

            if np.array(flags).all():
                return [datetime.combine(self.datestamp, self.timestamp),
                        self.lat, self.lon, self.alt]

        logger.warning(f'No GPS fix acquired after {time_to_wait} seconds')

    def close(self):
        """Close the connection."""
        self.running = False
        self.thread.join()
        self.serial_port.close()
        logger.info('GPS serial connection closed')
