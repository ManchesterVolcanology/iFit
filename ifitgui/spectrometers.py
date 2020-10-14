import logging
import numpy as np

try:
    import seatease.spectrometers as st
except ImportError:
    msg = 'seatease module not found, unable to use virtual spectrometer'
    logging.debug(msg)

import seabreeze.spectrometers as sb


class Spectrometer():

    def __init__(self, serial=None, integration_time=100, coadds=10):
        '''
        Wrapper around the python-seabreeze library for controlling Ocean
        Optics spectrometers. For more information see:
        https://github.com/ap--/python-seabreeze

        Parameters
        ----------
        serial : string, optional, default=None
            The serial number of the spectrometer to use. If None connects to
            the first available device. If "virtual" generates a virtual
            spectrometer using the seateeze library

        integration_time : int, optional, default=100
            The integration time of the spectrometer in milliseconds

        coadds : int, optional, default=10
            The number of individual spectra to average for each measurement

        '''

        # Connect to the spectrometer
        if serial == 'virtual':
            self.spec = st.Spectrometer.from_first_available()

        else:
            self.spec = sb.Spectrometer.from_serial_number(serial=serial)

        # Set the initial integration time and coadds
        self.update_coadds(coadds)
        self.update_integration_time(integration_time)

    def update_integration_time(self, integration_time):
        '''Update the spectrometer integrations time (ms)'''

        self.integration_time = integration_time
        self.spec.integration_time_micros(integration_time*1000)
        logging.info(f'Updated integration time to {integration_time} ms')

    def update_coadds(self, coadds):
        '''Update the number of coadds to average each spectrum over'''

        self.coadds = coadds
        logging.info(f'Updated coadds to {coadds}')

    def get_spectrum(self):
        '''Read a spectrum from the spectrometer'''

        # Get the wavelengths
        x = self.spec.wavelengths()

        # Create an empty array to hold the spectra data
        y_arr = np.zeros([self.coadds, len(x)])

        for n in range(self.coadds):
            y_arr[n] = self.spec.intensities()

        y = np.average(y_arr, axis=0)

        return np.row_stack(([x, y]))


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    spec = Spectrometer('virtual')

    x, y0 = spec.get_spectrum()

    spec.update_integration_time(200)
    spec.update_coadds(1)

    x, y1 = spec.get_spectrum()

    plt.plot(x, y0)
    plt.plot(x, y1)
    plt.show()
