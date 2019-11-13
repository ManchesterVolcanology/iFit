# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 09:13:39 2019

@author: mqbpwbe2
"""

import numpy as np
import matplotlib.pyplot as plt

# Import BrO data
bro = np.loadtxt('BrO_Cross_298K.txt', skiprows=12)

# As BrO is in terms of wavenumber(cm^-1), need to convert to nm
bro[:,0] = np.divide(1e7, bro[:,0])

# Reverse array so wavelength is accending
bro = bro[::-1]

plt.plot(bro[:,0], bro[:,1])
plt.show()

h = 'BrO Cross Section at 298K\n' + \
    'Reference: O. C. Fleischmann, J. P. Burrows, University of Bremen, 17.3.2000\n' + \
    'Wavelength (nm), Cross section (cm^2/molecule)'
np.savetxt('BrO_298K.txt', bro)