"""Provides an interface to a USB GPS device."""
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

    def __init__(self, comport, filename=None, baudrate=4800,
                 parity=serial.PARITY_NONE, stopbits=serial.STOPBITS_ONE,
                 bytesize=serial.EIGHTBITS):
        """Initialize."""
        self.serial_port = serial.Serial(comport, baudrate=baudrate,
                                         parity=parity, stopbits=stopbits,
                                         bytesize=bytesize)

        self.filename = filename
        self.timestamp = None
        self.lat = np.nan
        self.lon = np.nan
        self.alt = np.nan
        self.running = True

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

            # Unpak altitude
            alt_str = data[9]
            if alt_str != '':
                alt = float(alt_str)
                alt_unit = data[10]

                # Convert from feet to meters if required
                if alt_unit == 'F':
                    alt = 0.3048 * alt
                self.alt = alt
        except ValueError:
            logger.warn('Error parsing GPS string.')

    def close(self):
        """Close the connection."""
        self.running = False
        self.thread.join()
        self.serial_port.close()
