# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 14:05:44 2017

@author: mqbpwbe2
"""
import matplotlib.pyplot as plt
from ifit_lib.find_nearest import extract_window


#========================================================================================
#=======================================click_zoom=======================================
#========================================================================================

# Create a function to return mouse clicks
def onclick(event):
    
    # Read data of the mouse click
    click_x = event.xdata
    
    # Append to the global coords variable
    #global coords
    coords.append(click_x)
    
    # Exit after 2 clicks
    if len(coords) == 2:
        plt.close()
        
    return
    
def xclick(event):
    
    # Read data of the mouse click and close
    val = event.xdata
    
    plt.close()
    
    return

def yclick(event):
    
    # Read data of the mouse click
    click_x = event.xdata
    
    # Append to the global coords variable
    #global coords
    coords.append(click_x)
    
    # Exit after 2 clicks
    if len(coords) == 2:
        plt.close()
        
    return
    
#========================================================================================
#=======================================click_zoom=======================================
#========================================================================================

def click_zoom(fig, x_axis):
    
    '''
    Function to find the nearest value from a graph upon a mouse click

    INPUTS
    ------
    fig: figure to interact with
    x_axis: x data of the dataset to be cut down
    
    OUTPUTS
    -------
    ind: tuple of indicies to cut the data set to the desired window
    '''
    
    # Create empty array to store coords
    global coords
    coords = []

    # Connect to plot
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    
    # Disconnect from figure
    fig.canvas.mpl_disconnect(cid)
    
    # Extract the desired window
    x, ind1, ind2 = extract_window(x_axis, coords[0], coords[1])
    
    return ind1, ind2+1

#========================================================================================
#=======================================click_val========================================
#========================================================================================

def click_x(fig):
    global val
    
    # Connect to the plot
    cid = fig.canvas.mpl_connect('button_press_event', xclick)
    plt.show()
    
    # Disconnect from figure
    fig.canvas.mpl_disconnect(cid)
    
    return val
    
def click_y(fig):
    global val
    
    # Connect to the plot
    cid = fig.canvas.mpl_connect('button_press_event', yclick)
    plt.show()
    
    # Disconnect from figure
    fig.canvas.mpl_disconnect(cid)
    
    return val
    
    
    
    
    
    