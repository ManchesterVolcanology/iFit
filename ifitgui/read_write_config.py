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

        try:
            gui.spec_type.set(config['spec_type'])
            gui.w_lo.set(config['w_lo'])
            gui.w_hi.set(config['w_hi'])
            gui.model_padding.set(config['model_padding'])
            gui.model_spacing.set(config['model_spacing'])
            gui.dark_flag.set(config['dark_flag'])
            gui.flat_flag.set(config['flat_flag'])
            gui.stray_flag.set(config['stray_flag'])
            gui.update_flag.set(config['update_flag'])
            gui.resid_limit.set(config['resid_limit'])
            gui.ils_mode.set(config['ils_mode'])
            gui.ils_path.set(config['ils_path'])
            gui.flat_path.set(config['flat_path'])
            gui.wl_calib.set(config['wl_calib'])
            gui.frs_path.set(config['frs_path'])
            gui.graph_flag.set(config['graph_flag'])
            gui.graph_param.set(config['graph_param'])
            gui.scroll_flag.set(config['scroll_flag'])
            gui.graph_data_n.set(config['graph_data_n'])
            gui.spec_fnames = config['zspec_fnames']
            gui.dark_fnames = config['zdark_fnames']
            gui.spec_ent.set(f'{len(gui.spec_fnames)} spectra selected')
            gui.dark_ent.set(f'{len(gui.dark_fnames)} dark spectra selected')
            gui.save_path.set(config['save_path'])

            gui.gastable.set_params(config['gas_params'])
            gui.polytable.set_params(config['poly_params'])
            gui.offsettable.set_params(config['offset_params'])
            gui.shifttable.set_params(config['shift_params'])

            logging.info('Config file loaded')

        except KeyError:
            logging.warning('Problem loading config file')

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
    config = {'spec_type':     gui.spec_type.get(),
              'zspec_fnames':  gui.spec_fnames,
              'zdark_fnames':  gui.dark_fnames,
              'save_path':     gui.save_path.get(),
              'w_lo':          gui.w_lo.get(),
              'w_hi':          gui.w_hi.get(),
              'model_padding': gui.model_padding.get(),
              'model_spacing': gui.model_spacing.get(),
              'dark_flag':     gui.dark_flag.get(),
              'flat_flag':     gui.flat_flag.get(),
              'stray_flag':    gui.stray_flag.get(),
              'update_flag':   gui.update_flag.get(),
              'resid_limit':   gui.resid_limit.get(),
              'ils_mode':      gui.ils_mode.get(),
              'ils_path':      gui.ils_path.get(),
              'flat_path':     gui.flat_path.get(),
              'wl_calib':      gui.wl_calib.get(),
              'frs_path':      gui.frs_path.get(),
              'gas_params':    gas_params,
              'poly_params':   poly_params,
              'offset_params': offset_params,
              'shift_params':  shift_params,
              'graph_flag':    gui.graph_flag.get(),
              'graph_param':   gui.graph_param.get(),
              'scroll_flag':   gui.scroll_flag.get(),
              'graph_data_n':  gui.graph_data_n.get()
              }

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






