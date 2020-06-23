# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 12:00:05 2018

@author: mqbpwbe2
"""

import numpy as np
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


from ifitgui.gui_functions import select_files


# =============================================================================
# Make Input 
# =============================================================================

def make_input(frame, text, var, input_type, row, column, padx = 5, pady = 5,
               command = None, sticky = ['NSEW', None], state=None,
               label_font = ('Verdana', 10), width = None, options = None,
               vals = [0, 1000], increment = 1, rowspan = None,
               columnspan = None):

    '''
    Function to build GUI inputs consisting of a label and an input.

    **Parameters**

    frame : tk.Frame or tk.LabelFrame
        Container in which to place the object
    text : str
        Text to display in the label
    var : tk variable
        Variable to assosiate with the input
    input_type : str
        Type of input to use. Must be one of Entry, Spinbox, OptionMenu,
        Checkbutton or Label
    row : int
        Row number, will be the same for label and input
    column : int
        Column number, the input will be column + 1
    padx : int, optional, default = 5
        X-padding to apply to the label and input.
    pady : int, optional, default = 5
        Y-padding to apply to the label and input.
    command : func, optional, default = 5
        function to run on change to the input value
    sticky : str or tuple of strings, optional, default = None
        Direction to stick the object (compass direction). If given as a tuple
        the first corresponds to the label and the second to the entry.
    state : str, optional, default=None
        Sets the state of the input
    label_font : tuple, optional, default = ('Verdana', 8)
        Font tuple in the form (font, size) for the label.
    width : float, optional, default = None
        Width of the entry. Default is None.
    options : list, optional, default = None
        List of options for an Option Menu. Default is None
    values : tuple or list, optional, default = (0, 10000)
        Sets the range of values for a spinbox. If two values are give it sets
        the limits (from, to)
    increment : int, optional, default = 1
        Value spacing for a spinbox.
    rowspan : int, optional, default = None
        Number of rows the entry will span.
    columnspan : int, optional, default = None
        Number of columns the entry will span.

    **Returns**

    entry : tk object
        Input entry object, type depends on the input_type
    '''

    # Unpack stickyness
    if sticky == None:
        label_sticky = None
        entry_sticky = None
    elif len(sticky) == 2:
        label_sticky = sticky[0]
        entry_sticky = sticky[1]
    else:
        label_sticky = sticky
        entry_sticky = sticky

    # Create the input label
    if text != None:
        label = ttk.Label(frame, text = text, font = label_font)
        label.grid(row=row, column=column, padx=padx, pady=pady,
                   sticky=label_sticky)

    # Check that the entry type is valid
    if input_type not in ['Entry', 'Spinbox', 'OptionMenu', 'Label',
                          'Checkbutton', 'Combobox']:
        raise TypeError('Data entry type "' + input_type + '" not recognised')


    # Normal entry
    if input_type == 'Entry':
        if command == None:
            validate = None
        else:
            validate = "focusout"
        entry = ttk.Entry(frame, textvariable = var, width = width,
                          validate = validate, validatecommand = command)


    # Spinbox
    if input_type == 'Spinbox':

        # Check if range is from:to or a list
        if len(vals) == 2:
            entry = tk.Spinbox(frame, textvariable = var, width = width,
                               from_ = vals[0], to = vals[1],
                               command = command, increment = increment)

        else:
            entry = tk.Spinbox(frame, textvariable = var, width = width,
                               values = vals, command = command,
                               increment = increment)

        # Set first value
        #entry.update(var.get())


    # Option Menu
    if input_type == 'OptionMenu':
        entry = ttk.OptionMenu(frame, var, *options, command = command)
        entry.config(width = width)


    # Label
    if input_type == 'Label':
        entry = ttk.Label(frame, textvariable = var, width = width)


    # Checkbutton
    if input_type == 'Checkbutton':
        entry = ttk.Checkbutton(frame, variable = var)
        
    # Combobox
    if input_type == 'Combobox':
        entry = ttk.Combobox(frame, textvariable=var, width=width,
                             postcommand=command, state=state)
        entry['values'] = options


    # Add entry to the frame
    entry.grid(row=row, column=column+1, padx=padx, pady=pady,
               sticky=entry_sticky, rowspan=rowspan, columnspan=columnspan)

    return entry

# =============================================================================
# Make Table 
# =============================================================================

class ParamTable(tk.Frame):

    '''
    Class to create a parameter table for the gas parameters and xsecs

    **Parameters**

    parent : tk.Tk object
        The gui object

    font : tuple of font and fontsize, optional, default=("Verdana", 10)
        The font for the entries

    header_font : tuple of font and fontsize, optional, default=("Verdana", 14)
        The font for the header
    '''

    def __init__(self, parent, font=('Verdana',10), header_font=('Verdana',14)):

        tk.Frame.__init__(self, parent)

        self.font = font
        self.parent = parent
        self.row = 2

        self.cols = [1,3,5]

        # Make a label
        ttk.Label(self, text='Name', font=header_font
                  ).grid(row=0, column=0, padx=2, pady=2, sticky='NS')
        ttk.Label(self, text='Value', font=header_font
                  ).grid(row=0, column=2, padx=2, pady=2, sticky='NS')
        ttk.Label(self, text='Vary', font=header_font
                  ).grid(row=0, column=4, padx=2, pady=2, sticky='NS')
        ttk.Label(self, text='Xsec Path', font=header_font
                  ).grid(row=0, column=6, padx=2, pady=2, sticky='NS')

        ttk.Separator(self, orient='horizontal'
                      ).grid(column=0, row=1, columnspan=7, sticky='we',
                             padx=10, pady=5)

        for col in self.cols:
            ttk.Separator(self, orient='vertical'
                          ).grid(column=col,
                                 row=0,
                                 sticky='ns',
                                 padx=10,
                                 pady=0,
                                 rowspan=2)

        ttk.Button(self, text="Add Parameter", command=self.add_row
                   ).grid(row=0, column=7, padx=5, pady=5, sticky='W',
                   columnspan=2)

        self._params = []
        self._widgets = []

        self.add_row()



    def clear(self):
        '''Clear the table'''
        for row in self._widgets:
            for w in row:
                w.destroy()

        self._params = []
        self._widgets = []


    def add_row(self, name = '', value = 0, vary = True, xpath = ''):
        '''Adds a row to the table'''

        row_n = tk.IntVar(self, value = self.row-2)
        col_n = 0

        # Make a name entry
        name_h = tk.StringVar(self.parent, value=name)
        name_w = ttk.Entry(self,
                           textvariable = name_h,
                           width=10,
                           font=self.font)
        name_w.grid(row=self.row, column=col_n, padx=2, pady=2)
        col_n =+ 2

        # Make a value entry
        value_h = tk.DoubleVar(self.parent, value=value)
        value_w = ttk.Entry(self,
                            textvariable=value_h,
                            width=10,
                            font=self.font)
        value_w.grid(row=self.row, column=col_n, padx=2, pady=2)
        col_n += 2

        # Make a vary checkbutton
        vary_h = tk.BooleanVar(self.parent, value=vary)
        vary_w = ttk.Checkbutton(self, variable=vary_h)
        vary_w.grid(row=self.row, column=col_n, padx=2, pady=2)
        col_n += 2

        # Make a xpath entry
        xpath_h = tk.StringVar(self.parent, value=xpath)
        xpath_w = ttk.Entry(self,
                            textvariable=xpath_h,
                            width=20,
                            font=self.font)
        xpath_w.grid(row=self.row, column=col_n, padx=2, pady=2)
        col_n += 1

        # Add a browse button
        b_button = ttk.Button(self, text = "Browse", width = 8,
                   command = lambda: select_files(single_file = True,
                                                  holder = xpath_h))
        b_button.grid(row=self.row, column=col_n, padx=5, pady=5,
                      sticky='W')
        col_n += 1

        # Add a remove button
        x_button = ttk.Button(self, text = "X", width = 2,
                              command = lambda: self.remove_row(row_n.get(),
                                                                *[name_w,
                                                                  value_w,
                                                                  vary_w,
                                                                  xpath_w,
                                                                  b_button,
                                                                  x_button]))

        x_button.grid(row=self.row, column=col_n, padx=5, pady=5, sticky='W')

        for col in self.cols:
            ttk.Separator(self, orient='vertical'
                          ).grid(column=col,
                                 row=self.row,
                                 sticky='ns',
                                 padx=10,
                                 pady=0)
        self.row += 1

        self._params.append([name_h, value_h, vary_h, xpath_h])
        self._widgets.append([name_w, value_w, vary_w, xpath_w, b_button,
                              x_button])



    def set_params(self, params):
        '''Update the parameter table'''

        self.clear()
        self.row = 2

        for row in params:
            try:
                self.add_row(row[0], row[1], row[2], row[3])
            except IndexError:
                pass


    def remove_row(self, n, *widgets):
        '''Removes a single row'''
        for widget in widgets:
            widget.destroy()

        del(self._params[n])
        del(self._widgets[n])
        # self.row -= 1

#==============================================================================
#============================= Polynomial Table ===============================
#==============================================================================

class PolyTable(tk.Frame):

    '''
    Class to create a parameter table for the polynomials

    **Parameters**

    parent : tk.Tk object
        The gui object

    font : tuple of font and fontsize, optional, default=("Verdana", 10)
        The font for the entries

    header_font : tuple of font and fontsize, optional, default=("Verdana", 14)
        The font for the header
    '''

    def __init__(self, parent, font=('Verdana',10), header_font=('Verdana',14)):

        tk.Frame.__init__(self, parent)

        self.font = font
        self.parent = parent
        self.row = 2

        self.cols = [1,3]

        # Make a label
        ttk.Label(self, text='Coefficient', font=header_font
                  ).grid(row=0, column=0, padx=2, pady=2, sticky='NS')
        ttk.Label(self, text='Value', font=header_font
                  ).grid(row=0, column=2, padx=2, pady=2, sticky='NS')
        ttk.Label(self, text='Vary', font=header_font
                  ).grid(row=0, column=4, padx=2, pady=2, sticky='NS')

        ttk.Separator(self, orient='horizontal'
                      ).grid(column=0, row=1, columnspan=7, sticky='we',
                             padx=10, pady=5)

        for col in self.cols:
            ttk.Separator(self, orient='vertical'
                          ).grid(column=col,
                                 row=0,
                                 sticky='ns',
                                 padx=10,
                                 pady=0,
                                 rowspan=2)

        ttk.Button(self, text="Add Order", command=self.add_row
                   ).grid(row=0, column=7, padx=5, pady=5, sticky='W',
                   columnspan=2)

        self._params = []
        self._widgets = []

        self.add_row()



    def clear(self):
        '''Clear the table'''
        for row in self._widgets:
            for w in row:
                w.destroy()

        self._params = []
        self._widgets = []


    def add_row(self, value = 0, vary = True):

        row_n = tk.IntVar(self, value = self.row-2)
        col_n = 0

        # Make a name entry
        name_h = tk.StringVar(self.parent, value=f'C{self.row-2}')
        name_w = ttk.Label(self,
                           textvariable = name_h,
                           width=10,
                           font=self.font)
        name_w.grid(row=self.row, column=col_n, padx=2, pady=2, sticky = 'NS')
        col_n =+ 2

        # Make a value entry
        value_h = tk.DoubleVar(self.parent, value=value)
        value_w = ttk.Entry(self,
                            textvariable=value_h,
                            width=10,
                            font=self.font)
        value_w.grid(row=self.row, column=col_n, padx=2, pady=2)
        col_n += 2

        # Make a vary checkbutton
        vary_h = tk.BooleanVar(self.parent, value=vary)
        vary_w = ttk.Checkbutton(self, variable=vary_h)
        vary_w.grid(row=self.row, column=col_n, padx=2, pady=2)
        col_n += 2

        # Add a remove button
        x_button = ttk.Button(self, text = "X", width = 2,
                              command = lambda: self.remove_row(row_n.get(),
                                                                *[name_w,
                                                                  value_w,
                                                                  vary_w,
                                                                  x_button]))

        x_button.grid(row=self.row, column=col_n, padx=5, pady=5, sticky='W')

        for col in self.cols:
            ttk.Separator(self, orient='vertical'
                          ).grid(column=col,
                                 row=self.row,
                                 sticky='ns',
                                 padx=10,
                                 pady=0)
        self.row += 1

        self._params.append([name_h, value_h, vary_h])
        self._widgets.append([name_w, value_w, vary_w, x_button])



    def set_params(self, params):
        '''Update the parameter table'''

        self.clear()
        self.row = 2

        for row in params:

            try:
                self.add_row(row[1], row[2])
            except IndexError:
                pass


    def remove_row(self, n, *widgets):
        '''Removes a single row'''
        for widget in widgets:
            widget.destroy()

        self._params[n] = []
        self._widgets[n] = []
        self.row -= 1


#==============================================================================
#================================ Make Figure =================================
#==============================================================================

class GuiFigure():

    '''
    Class to generate the graphing figure

    **Parameters**

    gridlines : bool, optional, default = True
        If true then the graphs are created with gridlines

    fontsize : int, optional, default = 10
        The fontsize to use on the graph labels
    '''

    def __init__(self,
                 fig_kwargs={},
                 grid=[1,1],
                 axes_info = [],
                 gridlines=True,
                 fontsize=10):

        self.fontsize = fontsize

        # Make the figure and layout
        self.fig = plt.figure(**fig_kwargs)
        gs = gridspec.GridSpec(grid[0], grid[1])

        # Create a list to hold the axes
        self.ax_list = []

        # Create plot axes
        for ax_info in axes_info:
            ax = Axis(self.fig, gs, **ax_info)

            self.ax_list.append(ax)

        # Make it look nice
        plt.tight_layout()

        # Make a flag that will let a rearrange for the first plot, but fix
        #  after to avoid the graphs jumping about
        self.rearange_flag = True


    def set_labels(self, ax, xlabel=None, ylabel=None):
        '''Sets the axis labels'''

        if xlabel != None:
            ax.set_xlabel(xlabel, fontsize = self.fontsize)

        if ylabel != None:
            ax.set_ylabel(ylabel, fontsize = self.fontsize)

    def update_plots(self, data):
        '''Updates the axes with supplied data'''

        # Check that the number of lines matches the data provided
        num_lines = np.sum([len(ax.lines) for ax in self.ax_list])

        if len(data) != num_lines:
            raise ValueError('Number of data must match number of lines')

        # Create a counter
        n = 0

        # Cycle through the axis
        for ax in self.ax_list:

            # Create a flag to control if the axis is rescaled
            relim_flag = False

            for line in ax.lines:

                # For each line unpack and set the data
                x, y = data[n]
                line.set_data(x, y)

                # Check if the data is all nans. This is required to avoid
                #  rescaling an empty plot

                # For a numpy array
                try:
                    if ~np.isnan(y).all():
                        relim_flag = True

                # For a pandas dataframe column
                except TypeError:
                    if not y.isnull().all():
                        relim_flag = True

                # Update the loop counter
                n += 1

            # If any data is plotted rescale the axes
            if relim_flag:
                ax.axis.relim()
                ax.axis.autoscale_view()

        # Make the graphs a nice layout but only for the first plot
        if self.rearange_flag:
            plt.tight_layout()
            self.rearange_flag = False

class Axis():

    def __init__(self, fig, gs,
                 loc=[0,0],
                 lines=[],
                 ax_labels=[None, None],
                 grid=False):

        # Add the axis
        self.axis = fig.add_subplot(gs.new_subplotspec((loc[0], loc[1]),
                                                        rowspan=loc[2],
                                                        colspan=loc[3]))

        # Create an array to hold the lines for the axis
        self.lines = []

        # Create a flag to tell the program whether or not to add a legend
        legend_flag = False

        # Create each line for the plot
        for line_info in lines:
            line, = self.axis.plot([], [], **line_info)

            # If there is a label turn on the legend flag
            if 'label' in line_info: legend_flag = True

            # Add the line to the lines array
            self.lines.append(line)

        # If labels are provided add a legend
        if legend_flag: self.axis.legend()

        # If gidlines are requested add these
        if grid: self.axis.grid(ls='--')

        # Set the axis labels
        self.axis.set_xlabel(ax_labels[0])
        self.axis.set_ylabel(ax_labels[1])