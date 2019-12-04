# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 15:29:39 2019

@author: mqbpwbe2
"""

import os
import yaml
import logging
from tkinter import filedialog as fd

#==============================================================================
#================================ load_config =================================
#==============================================================================

def read_config(gui, fpath = None):
    '''Read the config file'''

    if fpath == None:
        fpath = fd.askopenfilename()

    # Open the config file
    try:
        with open(fpath, 'r') as ymlfile:

            config = yaml.load(ymlfile, Loader = yaml.FullLoader)
            
        for label in config:
            
            if label == 'gas_params':
                gui.gastable.set_params(config['gas_params'])
                
            elif label == 'poly_params':
                gui.polytable.set_params(config['poly_params'])
                
            elif label == 'offset_params':
                gui.offsettable.set_params(config['offset_params'])
                
            elif label == 'shift_params':
                gui.polytable.set_params(config['poly_params'])
                
            elif label == 'zspec_fnames':
                gui.spec_fnames = config['zspec_fnames']
                gui.spec_ent.set(f'{len(gui.spec_fnames)} files selected')
                
            elif label == 'zdark_fnames':
                gui.dark_fnames = config['zdark_fnames']
                gui.dark_ent.set(f'{len(gui.dark_fnames)} files selected')
            
            else:
                gui.widgets[label].set(config[label])

        logging.info('Config file loaded')

    except FileNotFoundError:
        logging.warn('Unable to load config file')
        config = {}

    return config

#==============================================================================
#=============================== write_config =================================
#==============================================================================

def write_config(gui, asksavepath=True):
    '''Write the config file'''

    # Pull the configuration of the parameter tables
    gas_params = []
    for row in  gui.gastable._params:
        gas_params.append([par.get() for par in row])
    poly_params = []
    for row in gui.polytable._params:
        poly_params.append([par.get() for par in row])
    offset_params = []
    for row in gui.offsettable._params:
        offset_params.append([par.get() for par in row])
    shift_params = []
    for row in gui.shifttable._params:
        shift_params.append([par.get() for par in row])

    # Populate the config dictionary
    config = {'gas_params':    gas_params,
              'poly_params':   poly_params,
              'offset_params': offset_params,
              'shift_params':  shift_params,
              'zspec_fnames':  gui.spec_fnames,
              'zdark_fnames':  gui.dark_fnames
              }
    
    for label in gui.widgets:
        config[label] = gui.widgets[label].get()

    if asksavepath:
        # Get the desired path
        fpath = fd.asksaveasfilename(defaultextension = '.yaml')

        # Save the config file
        with open(fpath, 'w') as outfile:

            yaml.dump(config, outfile)

            logging.info('Config file saved')

    # Also save the config to the default config location to load on startup
    try:

        if not os.path.exists('bin'):
            os.makedirs('bin')

        with open('bin/config.yaml', 'w') as outfile:
            yaml.dump(config, outfile)

    except FileNotFoundError:
        logging.info('Default config could not be saved')






