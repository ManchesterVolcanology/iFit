# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 17:47:20 2019

@author: mqbpwbe2
"""
import copy
import numpy as np
from collections import OrderedDict

#==============================================================================
#================================= Parameters =================================
#==============================================================================

class Parameters(OrderedDict):
    '''
    An ordered dictionary of all the Paramter objects that will be included in
    the forward model. Each Parameter has a single entry with a string label,
    value and boolian control on whether it is varied in the fit.

    Based on the Parameters class of lmfit.
    '''

    def __init__(self, *args, **kwargs):

        self.update(*args, **kwargs)

    def add(self, name, value=0, vary=True, xpath=None, lo_bound=-np.inf, 
            hi_bound=np.inf):
        '''
        Method to add a Parameter to the Parameters object.

        **Parameters**

        name : str
            Identifier string for the parameter. Each must be unique

        value : float
            The initial numerical parameter value

        vary : bool, optional
            If True then the parameter is fitted. Otherwise it is fixed to its
            value
        '''

        self.__setitem__(name, Parameter(name     = name,
                                         value    = value,
                                         vary     = vary,
                                         xpath    = xpath,
                                         lo_bound = lo_bound,
                                         hi_bound = hi_bound))

    def add_many(self, param_list):
        '''
        Method to add multiple Parameters to the Parameters object.

        **Parameters**

        param_list : list of Parameter like objects
        '''

        for param in param_list:

            self.__setitem__(param.name, param)


    def update_values(self, new_values):
        '''
        Updates the values of each Parameter in order
        '''
        n = 0
        for name in self:
            if self[name].vary:
                self[name].set(value = new_values[n])
                n += 1


    def valuesdict(self):
        '''Return an ordered dictionary of all parameter values'''
        return OrderedDict((p.name, p.value) for p in self.values())


    def fittedvaluesdict(self):
        '''Return an ordered dictionary of fitted parameter values'''
        return OrderedDict((p.name, p.value) for p in self.values() if p.vary)


    def popt_dict(self):
        ''''Return a dictionary of the optimised parameters'''
        return OrderedDict((p.name, p.fit_val) for p in self.values() if p.vary)


    def valueslist(self):
        '''Return a list of all parameter values'''
        return [(p.value) for p in self.values()]


    def fittedvalueslist(self):
        '''Return a list of the fitted parameter values'''
        return [(p.value) for p in self.values() if p.vary]


    def popt_list(self):
        '''Return a list of the optimised parameters'''
        return [(p.fit_val) for p in self.values() if p.vary]
    
    def bounds(self):
        '''Return a list of the low and high bounds'''
        return [[(p.lo_bound) for p in self.values() if p.vary],
                [(p.hi_bound) for p in self.values() if p.vary]]

    def make_copy(self):
        '''Returns a deep copy of the Parameters object'''
        return copy.deepcopy(self)


    def pretty_print(self, mincolwidth=10, precision=4, cols='basic'):
        '''
        Print the parameters in a nice way

        **Parameters**

        mincolwidth : int, optional, default = 10
            Minimum width of the columns.

        precision : int, optional, default = 4
            Number of significant figures to print to

        cols : str or list
            The columns to be printed. Either "all" for all columns, "basic" 
            for the name, value and if it is fixed or a list
            of the desired column names

        **Returns**

        msg : str
            The formatted message to print
        '''

        # Check colwidth isn't less than 7
        if mincolwidth <= 7:
            mincolwidth = 7

        # Make list of columns
        if cols == 'all': 
            cols = ['name', 'value', 'vary', 'fit_val', 'fit_err']
            
        if cols == 'basic':
            cols = ['name', 'value', 'vary']

        colwidth = [mincolwidth] * (len(cols))

        if 'name' in cols:
            i = cols.index('name')
            colwidth[i] = max([len(name) for name in self]) + 2

        if 'value' in cols:
            i = cols.index('value')
            colwidth[i] = max([len(f'{p.value:.{precision}g}') \
                               for p in self.values()]) + 2

        if 'vary' in cols:
            i = cols.index('vary')
            colwidth[i] = mincolwidth

        if 'fit_val' in cols:
            i = cols.index('fit_val')
            colwidth[i] = max([len(f'{p.fit_val:.{precision}g}') \
                               for p in self.values()]) + 2

        if 'fit_err' in cols:
            i = cols.index('fit_err')
            colwidth[i] = max([len(f'{p.fit_err:.{precision}g}') \
                               for p in self.values()]) + 2

        for n, w in enumerate(colwidth):
            if w < mincolwidth:
                colwidth[n] = mincolwidth

        title = ''
        for n, c in enumerate(cols):
            title += f'|{c:^{colwidth[n]}}'
        title += '|'

        msg = f'\n{"MODEL PARAMETERS":^{len(title)}}\n{title}\n' + \
              f'{"-"*len(title)}\n'

        for name, p in self.items():
            val = f'{p.value:.{precision}g}'
            fval = f'{p.fit_val:.{precision}g}'
            ferr = f'{p.fit_err:.{precision}g}'
            # lo_b = f'{p.lo_bound:.{precision}g}'
            # hi_b = f'{p.hi_bound:.{precision}g}'

            if p.vary: var = 'True'
            else: var = 'False'

            if 'name' in cols:
                msg += f'|{name:^{colwidth[cols.index("name")]}}'
            if 'value' in cols:
                msg += f'|{val:^{colwidth[cols.index("value")]}}'
            if 'vary' in cols:
                msg += f'|{var:^{colwidth[cols.index("vary")]}}'
            # if 'lo_bound' in cols:
            #     msg += f'|{lo_b:^{colwidth[cols.index("lo_bound")]}}'
            # if 'hi_bound' in cols:
            #     msg += f'|{hi_b:^{colwidth[cols.index("hi_bound")]}}'
            if 'fit_val' in cols:
                msg += f'|{fval:^{colwidth[cols.index("fit_val")]}}'
            if 'fit_err' in cols:
                msg += f'|{ferr:^{colwidth[cols.index("fit_err")]}}'

            msg += '|\n'

        return(msg)

#==============================================================================
#================================= Parameter ==================================
#==============================================================================

class Parameter(object):
    '''
    A parameter is a value that can be varied in the fit

    Each parameter has an assosiated name and value and can be set to either
    vary or be fixed in the model

    Based on the Parameter class of lmfit.

    **Parameters**

    name : str
        Identifier string for the parameter. Each must be unique

    value : float
        The initial numerical parameter value

    vary : bool, optional
        If True then the parameter is fitted. Otherwise it is fixed to its
        value
        
    xpath : str, optional
        The file path to the cross-section for this parameter
    '''

    def __init__(self, name, value, vary=True, xpath=None, lo_bound=-np.inf, 
                 hi_bound=np.inf, fit_val=np.nan, fit_err=np.nan):

        self.name = name
        self.value = value
        self.vary = vary
        self.xpath = xpath
        self.lo_bound = lo_bound
        self.hi_bound = hi_bound
        self.fit_val = fit_val
        self.fit_err = fit_err

    def set(self, value=None, vary=None, xpath=None, lo_bound=None, 
            hi_bound=None, fit_val=None, fit_err=None):
        '''Update the properties of a Parameter. All are None by default'''

        if value is not None:
            self.value = value

        if vary is not None:
            self.vary = vary

        if xpath is not None:
            self.xpath = xpath            

        if lo_bound is not None:
            self.lo_bound = lo_bound

        if hi_bound is not None:
            self.hi_bound = hi_bound

        if fit_val is not None:
            self.fit_val = fit_val

        if fit_err is not None:
            self.fit_err = fit_err