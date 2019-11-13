#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 15:06:30 2019

@author: mqbpwbe2
"""

import numpy as np

fname = 'Origional Spectra/O3_xsec.dat'

data = np.loadtxt(fname)

x = data[:,0]

temps = ['293','283','273','263','253','243','233','223','213','203','193']

for n, t in enumerate(temps):
    
    o3 = data[:,n+1]
    
    h = f'O3 absorption cross section at {t}K\n'+\
         'Source: IUP, MolSpec Lab, Serdyuchenko A., Gorshelev V, Weber M.\n'+\
         'Submitted to AMT, June 2013\n' + \
         'Wavelength (nm), Cross-section (cm2/molec)'
         
    np.savetxt(f'O3_{t}K.txt', np.column_stack((x, o3)), header = h)    