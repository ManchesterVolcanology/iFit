import numpy as np

#========================================================================================
#===================================== update graph =====================================
#========================================================================================

# Function to update a figure

# INPUTS: lines; array of data series to update
#         axes; array of axes that correspond to the lines (must be same length and same
#                order as lines)
#         fig; the figure object
#         new_data; array of new data to plot. Has the form [[x1, y1, xlims, ylims], ...]


def update_graph(lines, axes, canvas, new_data):
    
    '''
    Function to update a figure
    
    INPUTS
    ------
    lines, list
        The plots to update
    
    axes, list
        Axes that correspond to the lines (must be same length and order as lines)
    
    canvas, tkagg canvas object
        Canvas that holds the axes
        
    new_data, array
        New data to plot. Has the form [[x1, y1, x1lims, y1lims], 
                                        [x2, y2, x2lims, y2lims],...]
    '''

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

        try:
            # If auto, pad by 10% of range      
            if xlims[i] == 'auto':
                x_min = min(xdata) - abs(max(xdata) - min(xdata)) * 0.1
                x_max = max(xdata) + abs(max(xdata) - min(xdata)) * 0.1
                axes[i].set_xlim(x_min, x_max)
            
            # If false set as limits +/- 1
            elif xlims[i] == False:
                axes[i].set_xlim(min(xdata[i])-1, max(xdata[i])+1)
            
            # If fixed, fix value
            else:
                axes[i].set_xlim(xlims[i][0], xlims[i][1])
            
            # Do same for y axis
            if ylims[i] == 'auto':
                y_min = min(ydata) - abs(max(ydata) - min(ydata)) * 0.1
                y_max = max(ydata) + abs(max(ydata) - min(ydata)) * 0.1
                axes[i].set_ylim(y_min, y_max)
                
            elif ylims[i] == False:
                axes[i].set_ylim(min(ydata[i])-1, max(ydata[i])+1)
            
            else:
                axes[i].set_ylim(ylims[i][0], ylims[i][1]) 
                
        except ValueError:
            pass
        
    # Apply changes
    canvas.draw()