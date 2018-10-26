# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 12:00:05 2018

@author: mqbpwbe2
"""


from tkinter import ttk
import tkinter as tk

def make_input(frame, text, var, input_type, row, column, padx = 5, pady = 5, 
               command = None, sticky = ['NSEW', None], label_font = ('Verdana', 8),
               width = None, options = None, vals = [0, 10], increment = 1):
    
    '''
    Function to build GUI inputs consisting of a label and an input.
    
    INPUTS
    ------
    
    frame, tk.Frame or tk.LabelFrame
        Container in which to place the object
        
    text, str
        Text to display in the label
        
    var, tk variable
        Variable to assosiate with the input
        
    input_type, str
        Type of input to use. Must be one of Entry, Spinbox, OptionMenu, Checkbutton or
        Label
        
    row, int
        Row number, will be the same for label and input
        
    column, int
        Column number, the input will be column + 1
        
    padx, int (optional)
        X-padding to apply to the label and input. Default is 5
        
    pady, int (optional)
        Y-padding to apply to the label and input. Default is 5
        
    command, func
        function to run on change to the input value
        
    sticky, str or tuple of strings (optional)
        Direction to stick the object (compass direction). If given as a tuple the first
        corresponds to the label and the second to the entry. Default is None
        
    label_font, tuple (optional)
        Font tuple in the form (font, size) for the label. Default is ('Verdana', 8)
        
    width, float (optional)
        Width of the entry. Default is None.
        
    options, list (optional)
        List of options for an Option Menu. Default is None
        
    values, tuple or list (optional)
        Sets the range of values for a spinbox. If two values are give it sets the limits
        (from, to)
        
    increment, int
        Value spacing for a spinbox
        
    OUTPUTS
    -------
    label, tk.Label object
        Input label object
    
    entry, tk object
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
    label = ttk.Label(frame, text = text, font = label_font)
    label.grid(row=row, column=column, padx=padx, pady=pady, sticky=label_sticky)
    
    # Check that the entry type is valid
    if input_type not in ['Entry', 'Spinbox', 'OptionMenu', 'Label', 'Checkbutton']:
        raise TypeError('Data entry type "' + input_type + '" not recognised')
        
        
    # Normal entry
    if input_type == 'Entry':
        if command == None:
            validate = None
        else:
            validate = "focusout"
        entry = ttk.Entry(frame, textvariable = var, width = width, validate = validate,
                          validatecommand = command)
        
        
    # Spinbox
    if input_type == 'Spinbox':
        
        # Check if range is from:to or a list
        if len(vals) == 2:
            entry = tk.Spinbox(frame, textvariable = var, width = width, from_ = vals[0],
                               to = vals[1], command = command, increment = increment)
            
        else:
            entry = tk.Spinbox(frame, textvariable = var, width = width, values = vals, 
                               command = command, increment = increment)
            
        # Set first value
        #entry.update(var.get())
            
            
    # Option Menu
    if input_type == 'OptionMenu':
        entry = ttk.OptionMenu(frame, var, *options, command = command)
        entry.config(width = width)
        
    
    # Label
    if input_type == 'Label':
        entry = ttk.Label(frame, textvariable = var)
        
    
    # Checkbutton
    if input_type == 'Checkbutton':
        entry = ttk.Checkbutton(frame, variable = var)
        
        
    # Add entry to the frame
    entry.grid(row=row, column=column+1, padx=padx, pady=pady, sticky=entry_sticky)
    
    return label, entry
    
    
    