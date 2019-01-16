# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 15:09:45 2019

@author: mqbpwbe2
"""

import numpy as np
import matplotlib.pyplot as plt

x, y = np.loadtxt('so2_xsec.dat', unpack = True)

plt.plot(x, y)

x, y = np.loadtxt('SO2_293K.dat', unpack = True)

plt.plot(x, y)
plt.show()