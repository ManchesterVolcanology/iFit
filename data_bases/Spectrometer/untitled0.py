# -*- coding: utf-8 -*-
"""
Created on Thu May  2 11:40:47 2019

@author: mqbpwbe2
"""

import numpy as np

x, y = np.loadtxt('flat_FLMS02101.txt', unpack = True)

x += 1.5

np.savetxt('flat_FLMS02101_2.txt', np.column_stack((x, y)))