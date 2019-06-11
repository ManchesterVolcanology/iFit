# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 10:26:06 2019

@author: mqbpwbe2
"""

import numpy as np
import matplotlib.pyplot as plt
'''
wn, xs, se, ne, st, nb = np.loadtxt('SO2_298K_24000-44000.txt', skiprows = 20, unpack = True)

wl = np.divide(1e7, wn)

xs = np.multiply(xs, 1e-18)

wl = np.flip(wl, 0)
xs = np.flip(xs, 0)

np.savetxt('SO2_298K.txt', np.column_stack((wl, xs)))

plt.plot(wl, xs)
plt.show()
'''

wn, xs = np.loadtxt('data.fmt.txt', unpack = True)

wl = np.divide(1e7, wn)

xs = np.multiply(xs, 1e-19)

wl = np.flip(wl, 0)
xs = np.flip(xs, 0)

idx = np.where(wl > 280)
wl = wl[idx]
xs = xs[idx]

dx = np.array([wl[n+1] - wl[n] for n in range(len(wl) - 1)])

wl = np.delete(wl, np.where(dx == 0))
xs = np.delete(xs, np.where(dx == 0))

np.savetxt('SO2_295K.txt', np.column_stack((wl, xs)))

plt.plot(wl, xs)
plt.show()


x, y = np.loadtxt('SO2_295K.txt', unpack = True)

plt.plot(x, y)
plt.show()