import numpy as np

#========================================================================================
#======================================update graph======================================
#========================================================================================

# Function to update a figure

# INPUTS: lines; array of data series to update
#         axes; array of axes that correspond to the lines (must be same length and same
#                order as lines)
#         fig; the figure object
#         new_data; array of new data to plot. Has the form [[x1, y1, xlims, ylims], ...]


def update_graph(lines, axes, fig, new_data):

    # Unpack new data
    if len(np.shape(new_data)) > 1:
        xdata = new_data[:,0]
        ydata = new_data[:,1]
        xlims = new_data[:,2]
        ylims = new_data[:,3] 
    
    else:
        xdata = [new_data[0]]
        ydata = [new_data[1]]
        xlims = [new_data[2]]
        ylims = [new_data[3]]
        
        
    # Iterate plotting over each data series
    for i in range(len(lines)):

        # Update data points on the graph
        lines[i].set_xdata(xdata[i])
        lines[i].set_ydata(ydata[i])
        
        # Rescale the axes
        axes[i].relim()
        axes[i].autoscale_view()     
        
        
        # If no limits given, fix as max/min        
        if xlims[i] == False:
            # Manually fix the x axis limits
            axes[i].set_xlim(xdata[i][0], xdata[i][-1])
        else:
            axes[i].set_xlim(xlims[i][0], xlims[i][1])
            
        if ylims[i] == False:
            axes[i].set_ylim(ydata[i][0], ydata[i][-1])
        else:
            axes[i].set_ylim(ylims[i][0], ylims[i][1])            
        
    # Apply changes
    fig.canvas.draw()