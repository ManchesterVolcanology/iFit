import os
import logging
import platform
import traceback
import tkinter as tk
from tkinter import ttk
import tkinter.messagebox as tkMessageBox
import tkinter.scrolledtext as ScrolledText
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, \
                                              NavigationToolbar2Tk

from ifitgui.build_gui import make_input, GuiFigure, ParamTable, PolyTable
from ifitgui.gui_functions import stop, file_io, analysis_loop, test_spectrum,\
                                  connect_spectrometer, read_darks,\
                                  toggle_fitting
from ifitgui.read_write_config import read_config, write_config

# Define some fonts to use in the program
NORM_FONT = ('TkFixedFont', 10)
LARGE_FONT = ('TkFixedFont', 12, 'bold')


# =============================================================================
# -----------------------------------------------------------------------------
# Set up logging handler
# -----------------------------------------------------------------------------
# =============================================================================

class TextHandler(logging.Handler):

    def __init__(self, text, gui):

        # Run the regular Handler __init__
        logging.Handler.__init__(self)

        # Store a reference to the Text it will log to
        self.text = text
        self.gui = gui

    def emit(self, record):
        msg = self.format(record)

        def append():
            # Add message
            self.text.configure(state='normal')
            self.text.insert(tk.END, msg + '\n')
            self.text.configure(state='disabled')

            # Autoscroll to the bottom
            self.text.yview(tk.END)

        # This is necessary because we can't modify the Text from other threads
        self.text.after(0, append)

        # Force the gui to update
        self.gui.update()


# =============================================================================
# -----------------------------------------------------------------------------
# Main GUI interface
# -----------------------------------------------------------------------------
# =============================================================================

class myGUI(tk.Frame):

    # This class defines the graphical user interface

    def __init__(self, parent, *args, **kwargs):

        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.root = parent

        # Cause exceptions to report in a new window
        tk.Tk.report_callback_exception = self.report_callback_exception

        # Close program on closure of window
        self.root.protocol("WM_DELETE_WINDOW", self.handler)

        # Set the intial window size
        self.root.geometry("1350x800")

        # Make a dictionary to hold the gui widgets
        self.widgets = {}

        # Build the GUI
        self.build_gui()

        # Try loading the default config file
        read_config(self, fpath='bin/config.yaml')

# =============================================================================
# -----------------------------------------------------------------------------
# Build GUI
# -----------------------------------------------------------------------------
# =============================================================================

    def build_gui(self):

        # Build GUI
        self.root.title('iFit v3.2')

        # Set default button Style
        ttk.Style().configure('TButton', width=15, height=20, relief="flat")

        # Set the default style of notebooks
        ttk.Style().configure('TNotebook.Tab', padding=[40, 10],
                              font=LARGE_FONT, anchor="center")
        ttk.Style().configure('Label', anchor='center')

        # Set the window icon
        try:
            self.root.iconbitmap(default='bin/icon.ico')
        except tk.TclError:
            pass

        # Build a menubar to hold options for the user
        menubar = tk.Menu(self)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label='Save',
                             command=lambda: write_config(self, False))
        filemenu.add_command(label='Save as',
                             command=lambda: write_config(self))
        filemenu.add_command(label='Load',
                             command=lambda: read_config(self))
        filemenu.add_separator()
        menubar.add_cascade(label='File', menu=filemenu)
        self.root.config(menu=menubar)

# =============================================================================
#       Add scrollbars for the whole app
# =============================================================================

        # Create a canvas
        canv = tk.Canvas(self.root)
        canv.grid(row=0, column=0, sticky="NSEW")

        # Create the scrollbars
        vsb = ttk.Scrollbar(self.root, orient=tk.VERTICAL,
                            command=canv.yview)
        vsb.grid(row=0, column=1, sticky="NS")

        hsb = ttk.Scrollbar(self.root, orient=tk.HORIZONTAL,
                            command=canv.xview)
        hsb.grid(row=1, column=0, sticky="EW")

        # Configure the canvas to scroll with the scrollbars
        canv.configure(xscrollcommand=hsb.set, yscrollcommand=vsb.set)
        canv.bind('<Configure>',
                  lambda e: canv.configure(scrollregion=canv.bbox("all")))

        # Add bindings to use the mousewheel. Note these are platform dependent
        pltfrm = platform.system()

        # For windows
        if pltfrm == 'Windows':
            canv.bind_all('<MouseWheel>',
                          lambda e: canv.yview_scroll(-1*(e.delta//120),
                                                      "units"))
            canv.bind_all('<Shift-MouseWheel>',
                          lambda e: canv.xview_scroll(-1*(e.delta//120),
                                                      "units"))

        # For Linux
        elif pltfrm == 'Linux':
            canv.bind_all('<Button-4>',
                          lambda e: canv.yview_scroll(-1*(e.delta//120),
                                                      "units"))
            canv.bind_all('<Button-5>',
                          lambda e: canv.yview_scroll(-1*(e.delta//120),
                                                      "units"))
            canv.bind_all('<Shift-Button-4>',
                          lambda e: canv.xview_scroll(-1*(e.delta//120),
                                                      "units"))
            canv.bind_all('<Shift-Button-5>',
                          lambda e: canv.xview_scroll(-1*(e.delta//120),
                                                      "units"))

        # For MacOS
        elif pltfrm == 'Darwin':
            canv.bind_all('<MouseWheel>',
                          lambda e: canv.yview_scroll(-e.delta, "units"))
            canv.bind_all('<Shift-MouseWheel>',
                          lambda e: canv.xview_scroll(-e.delta, "units"))

        # Allow the window to be resized
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

        # Create another frame inside the canvas
        root_frame = tk.Frame(canv)

        # Add the new frame to a window in the canvas
        canv.create_window((0, 0), window=root_frame, anchor="nw")

# =============================================================================
#       Build containers
# =============================================================================

        # Create notebook to flick between real time and post analysis
        nb = ttk.Notebook(root_frame)
        page1 = ttk.Frame(nb)
        page2 = ttk.Frame(nb)

        # Create two frames, one for post analysis and one for real time
        nb.add(page2, text='Acquisition')
        nb.add(page1, text='Post Analysis')
        nb.grid(row=0, column=0, padx=10, pady=10, sticky='NW')

        # Create the frames for real-time analysis
        rt_frame = ttk.Frame(page2)
        rt_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

        # Create frames for post analysis setups
        setup_frame = ttk.Frame(page1)
        setup_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

        # Create a frame to hold the controls
        control_frame = tk.LabelFrame(root_frame, text='Control',
                                      font=LARGE_FONT)
        control_frame.grid(row=1, column=0, padx=10, pady=10, sticky="N")

        graph_settings = ttk.Notebook(root_frame)
        graph_frame = ttk.Frame(graph_settings)
        settings_frame = ttk.Frame(graph_settings)

        graph_settings.add(graph_frame, text='Graphs')
        graph_settings.add(settings_frame, text='Settings')

        graph_settings.grid(row=0, column=1, padx=10, pady=10, rowspan=2,
                            sticky="NW")

# =============================================================================
#       Realtime acquisition controls
# =============================================================================

        # Add control to connect to the spectrometer
        self.spec_id = tk.StringVar(rt_frame, value='Not connected')
        make_input(frame=rt_frame,
                   text='Spectrometer:',
                   row=0, column=0,
                   var=self.spec_id,
                   input_type='Label',
                   width=15,
                   sticky='W')
        self.connect_flag = tk.BooleanVar(rt_frame, value=False)
        self.con_b = ttk.Button(rt_frame, text="Connect", width=15,
                                command=lambda: connect_spectrometer(self))
        self.con_b.grid(row=0, column=2, padx=5, pady=5, sticky='W')

        # Add button to read a test spectrum
        ttk.Button(rt_frame, text='Test Spectrum', width=15,
                   command=lambda: test_spectrum(self)
                   ).grid(row=0, column=3, padx=5, pady=5, sticky='NS')

        # Create a control for the spectrometer integration time
        self.widgets["IntTime"] = tk.IntVar(rt_frame, value=100)
        make_input(frame=rt_frame,
                   text='Integration\nTime (ms):',
                   row=1, column=0,
                   var=self.widgets['IntTime'],
                   input_type='Spinbox',
                   vals=[0, 10000],
                   increment=50,
                   width=15,
                   sticky='W')

        def b_func1():
            self.spec.update_integration_time(self.widgets["IntTime"].get())
        ttk.Button(rt_frame, text='Update', width=15, command=b_func1
                   ).grid(row=1, column=2, padx=5, pady=5, sticky='W')

        # Create a control for the spectrometer coadds
        self.widgets["Coadds"] = tk.IntVar(rt_frame, value=10)
        make_input(frame=rt_frame,
                   text='Coadds:',
                   row=2, column=0,
                   var=self.widgets['Coadds'],
                   input_type='Spinbox',
                   vals=[0, 100],
                   increment=5,
                   width=15,
                   sticky='W')

        def b_func2():
            self.spec.update_coadds(self.widgets["Coadds"].get())
        ttk.Button(rt_frame, text='Update', width=15, command=b_func2
                   ).grid(row=2, column=2, padx=5, pady=5, sticky='W')

        # Add a button to turn fitting on and off
        self.fitting_flag = False
        self.fit_b = tk.Button(rt_frame, text='Fitting OFF', width=10,
                               fg='white', bg='red', font=LARGE_FONT,
                               command=lambda: toggle_fitting(self))
        self.fit_b.grid(row=1, column=3, padx=5, pady=5, rowspan=2,
                        sticky='NSEW')

        # Create a control for the number of darks to collect
        self.widgets["ndarks"] = tk.IntVar(rt_frame, value=10)
        make_input(frame=rt_frame,
                   text='No. Dark\nSpectra:',
                   row=3, column=0,
                   var=self.widgets['ndarks'],
                   input_type='Spinbox',
                   vals=[0, 100],
                   increment=1,
                   width=15,
                   sticky='W')

        # Add buttons to acquire test spectra and dark spectra
        ttk.Button(rt_frame, text='Read Darks', width=15,
                   command=lambda: read_darks(self)
                   ).grid(row=3, column=2, padx=5, pady=5, sticky='W')

        # Set the save path
        self.widgets['rt_save_path'] = tk.StringVar(self, value='')
        ttk.Label(rt_frame, text='Save Path:', font=NORM_FONT
                  ).grid(row=4, column=0, padx=5, pady=5, sticky='W')
        ttk.Entry(rt_frame,
                  font=NORM_FONT,
                  width=30,
                  text=self.widgets['rt_save_path']
                  ).grid(row=4, column=1, padx=5, pady=5, sticky='W',
                         columnspan=2)
        ttk.Button(rt_frame, text="Browse", width=10,
                   command=lambda: file_io(True,
                                           False,
                                           self.widgets['rt_save_path'],
                                           None,
                                           True,
                                           [('Comma Separated', '.csv')])
                   ).grid(row=4, column=3, padx=5, pady=5, sticky='W')

        # Frame to hold the control buttons
        rt_button_frame = ttk.Frame(rt_frame)
        rt_button_frame.grid(row=5, column=0, padx=10, pady=10, columnspan=6,
                             sticky="W")

        # Create button to start
        ttk.Button(rt_button_frame, text='Begin!',
                   command=lambda: analysis_loop(self, rt_flag=True)
                   ).grid(row=0, column=0, padx=25, pady=5)

        # Create button to stop
        self.stop_flag = False
        ttk.Button(rt_button_frame, text='Stop', command=lambda: stop(self)
                   ).grid(row=0, column=1, padx=25, pady=5)

# =============================================================================
#       Set up spectra selection
# =============================================================================

        # Create entry to select spectra type
        spec_options = ['iFit',
                        'iFit',
                        'Master.Scope',
                        'Spectrasuite',
                        'OpenSO2',
                        'FLAME',
                        'Basic']

        self.widgets['spec_type'] = tk.StringVar(setup_frame, value='iFit')
        make_input(frame=setup_frame,
                   text='Format:',
                   row=0, column=0,
                   var=self.widgets['spec_type'],
                   input_type='OptionMenu',
                   options=spec_options,
                   width=12,
                   sticky='W')

        self.spec_fnames = []
        self.dark_fnames = []

        # File dialouge for spectra
        message = 'No spectra selected'
        ttk.Label(setup_frame, text='Spectra:', font=NORM_FONT
                  ).grid(row=1, column=0, padx=5, pady=5, sticky='W')
        self.spec_ent = tk.StringVar(value=message)
        specfp_l = tk.Entry(setup_frame, font=NORM_FONT, width=30,
                            text=self.spec_ent)
        specfp_l.grid(row=1, column=1, padx=5, pady=5, sticky='W',
                      columnspan=2)
        ttk.Button(setup_frame, text="Browse", width=10,
                   command=lambda: file_io(holder=self.spec_fnames,
                                           entry=self.spec_ent)
                   ).grid(row=1, column=3, padx=5, pady=5, sticky='W')

        # File dialouge for darks
        ttk.Label(setup_frame, text='Darks:', font=NORM_FONT
                  ).grid(row=2, column=0, padx=5, pady=5, sticky='W')
        self.dark_ent = tk.StringVar(value=message)
        darkfp_l = tk.Entry(setup_frame, font=NORM_FONT, width=30,
                            text=self.dark_ent)
        darkfp_l.grid(row=2, column=1, padx=5, pady=5, sticky='W',
                      columnspan=2)
        ttk.Button(setup_frame, text="Browse", width=10,
                   command=lambda: file_io(holder=self.dark_fnames,
                                           entry=self.dark_ent)
                   ).grid(row=2, column=3, padx=5, pady=5, sticky='W')

        # Set the save path
        self.widgets['save_path'] = tk.StringVar(self, value='')
        ttk.Label(setup_frame, text='Save Path:', font=NORM_FONT
                  ).grid(row=3, column=0, padx=5, pady=5, sticky='W')
        ttk.Entry(setup_frame,
                  font=NORM_FONT,
                  width=30,
                  text=self.widgets['save_path']
                  ).grid(row=3, column=1, padx=5, pady=5, sticky='W',
                         columnspan=2)
        ttk.Button(setup_frame, text="Browse", width=10,
                   command=lambda: file_io(False,
                                           True,
                                           self.widgets['save_path'],
                                           None,
                                           True,
                                           [('Comma Separated', '.csv')])
                   ).grid(row=3, column=3, padx=5, pady=5, sticky='W')

        # Frame to hold the control buttons
        button_frame = ttk.Frame(setup_frame)
        button_frame.grid(row=4, column=0, padx=10, pady=10, columnspan=6,
                          sticky="W")

        # Create button to start
        ttk.Button(button_frame, text='Begin!',
                   command=lambda: analysis_loop(self, rt_flag=False)
                   ).grid(row=0, column=0, padx=25, pady=5)

        # Create button to stop
        self.stop_flag = False
        ttk.Button(button_frame, text='Stop', command=lambda: stop(self)
                   ).grid(row=0, column=1, padx=25, pady=5)

# =============================================================================
#       Create graphs
# =============================================================================

        # Set up the axes
        axes = [{'loc': [0, 0, 1, 1],
                 'grid': True,
                 'ax_labels': [None, 'Intensity'],
                 'lines': [{'marker': 'o', 'ls': '', 'c': 'C0',
                            'label': 'Data', 'ms': 4},
                           {'lw': 1.5, 'c': 'C1', 'label': 'Fit'}]},
                {'loc': [0, 1, 1, 1],
                 'grid': True,
                 'ax_labels': [None, 'Intensity'],
                 'lines': [{}]},
                {'loc': [1, 0, 1, 1],
                 'grid': True,
                 'ax_labels': ['Wavelength', 'Fit Residual'],
                 'lines': [{'marker': 'o', 'ms': 4}]},
                {'loc': [1, 1, 1, 1],
                 'grid': True,
                 'ax_labels': ['Wavelength', 'Optical Depth'],
                 'lines': [{'marker': 'o', 'c': 'C0', 'ms': 4},
                           {'lw': 1.5, 'c': 'C1'}]},
                {'loc': [2, 0, 1, 2],
                 'grid': True,
                 'ax_labels': ['Spectrum Number', 'Value'],
                 'lines': [{'marker': 'o', 'c': 'C0', 'ms': 4}]}
                ]

        # Build the figure
        self.figure = GuiFigure(grid=[3, 2],
                                fig_kwargs={'figsize': [7, 5]},
                                axes_info=axes)

        # Create the canvas to hold the graph in the GUI
        self.canvas = FigureCanvasTkAgg(self.figure.fig, graph_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=1, column=0, padx=10, pady=10)

        # Add matplotlib toolbar above the plot canvas
        toolbar_frame = tk.Frame(graph_frame, bg='black')
        toolbar_frame.grid(row=0, column=0, sticky='W', padx=5, pady=5)
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        toolbar.update()

# =============================================================================
#       Graph Controls
# =============================================================================

        # Create container for the graph settings
        graph_settings_frame = ttk.Frame(graph_frame)
        graph_settings_frame.grid(row=2, column=0, stick='W', padx=5, pady=5)

        row_n = 0
        col_n = 0

        # Control whether graphs are plotted
        self.widgets['graph_flag'] = tk.BooleanVar(graph_settings_frame,
                                                   value=True)
        make_input(frame=graph_settings_frame,
                   text='Show graphs?',
                   row=row_n, column=col_n,
                   var=self.widgets['graph_flag'],
                   input_type='Checkbutton')
        col_n += 3

        # Update graph param combobox when combobox is clicked
        def update_combobox():

            # Get the gas parameters
            rows = self.gastable._params
            options = []

            # Add the gases to the options
            for row in rows:
                options.append(row[0].get().strip())

            # Set the options
            graph_p = self.widgets['graph_param'].get()
            self.graph_cb['values'] = options
            self.widgets['graph_param'].set(graph_p)

        # Set parameter to graph
        self.widgets['graph_param'] = tk.StringVar(graph_settings_frame)
        self.graph_cb = make_input(frame=graph_settings_frame,
                                   text='Parameter\nto graph:',
                                   row=row_n, column=col_n,
                                   var=self.widgets['graph_param'],
                                   input_type='Combobox',
                                   width=10,
                                   command=update_combobox,
                                   options=[],
                                   state="readonly")
        col_n += 3

        # Control whether the graph scrolls
        self.widgets['scroll_flag'] = tk.BooleanVar(graph_settings_frame,
                                                    value=True)
        make_input(frame=graph_settings_frame,
                   text='Scroll graphs?',
                   row=row_n, column=col_n,
                   var=self.widgets['scroll_flag'],
                   input_type='Checkbutton')
        col_n += 3

        # Set number of measurements to scroll on the graph
        self.widgets['scroll_amt'] = tk.IntVar(graph_settings_frame, value=100)
        make_input(frame=graph_settings_frame,
                   text='No. Spectra\nto display:',
                   row=row_n, column=col_n,
                   var=self.widgets['scroll_amt'],
                   input_type='Spinbox',
                   vals=[0, 1.0e6],
                   increment=100,
                   width=10)
        col_n += 3

        # Create separator
        for col in [2, 5, 8]:
            ttk.Separator(graph_settings_frame, orient='vertical'
                          ).grid(column=col, row=0, rowspan=7, sticky='ns',
                                 padx=10, pady=5)

# =============================================================================
#       Create progress bar
# =============================================================================

        # Make a frame for the progress bar and status
        progress_frame = ttk.Frame(control_frame)
        progress_frame.grid(row=1, column=0, padx=5, pady=5, columnspan=6,
                            sticky='W')

        # Create progress bar
        self.progress = ttk.Progressbar(progress_frame, orient=tk.HORIZONTAL,
                                        length=300, mode='determinate')
        self.progress.grid(row=0, column=0, padx=5, pady=5, sticky='W')

        # Create status indicator
        self.status = tk.StringVar(progress_frame, value='Standby')
        status_e = ttk.Label(progress_frame, textvariable=self.status,
                             width=15)
        status_e.grid(row=0, column=1, padx=5, pady=5, sticky="W")

        # Make a frame for the progress bar and status
        amt_frame = ttk.Frame(control_frame)
        amt_frame.grid(row=2, column=0, padx=5, pady=5, columnspan=6,
                       sticky='W')

        # Create ouput for last so2 amount
        self.last_amt = tk.StringVar(amt_frame, value='-')
        make_input(frame=amt_frame,
                   text='Last amt:',
                   row=0, column=0,
                   var=self.last_amt,
                   input_type='Label',
                   sticky='W',
                   width=10,
                   padx=5)

        # Create ouput for last so2 error
        self.last_err = tk.StringVar(amt_frame, value='-')
        make_input(frame=amt_frame,
                   text='+/-',
                   row=0, column=2,
                   var=self.last_err,
                   input_type='Label',
                   sticky='W',
                   width=10,
                   padx=5)

# =============================================================================
#       Set up text output
# =============================================================================

        # Add text widget to display logging info
        st = ScrolledText.ScrolledText(control_frame, state='disabled',
                                       width=50, height=12)
        st.configure(font='TkFixedFont')
        st.grid(row=4, column=0, padx=10, pady=10, sticky="NW", columnspan=5)

        # Create textLogger
        text_handler = TextHandler(st, self)

        # Logging configuration
        if not os.path.isdir('bin/'):
            os.makedirs('bin/')
        logging.basicConfig(filename='bin/iFit.log',
                            level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S',
                            filemode='w')

        # Add the handler to logger
        logger = logging.getLogger()
        logger.addHandler(text_handler)

# =============================================================================
# =============================================================================
#       Program Settings
# =============================================================================
# =============================================================================

        # Make a notebook to hold the settings pannels
        settings_nb = ttk.Notebook(settings_frame)
        model_frame = ttk.Frame(settings_nb)
        spect_frame = ttk.Frame(settings_nb)
        param_frame = ttk.Frame(settings_nb)
        spec_frame = ttk.Frame(settings_nb)
        settings_nb.add(model_frame, text='Model')
        settings_nb.add(spect_frame, text='Spectrometer')
        settings_nb.add(param_frame, text='Parameters')
        settings_nb.add(spec_frame, text='Spectra')
        settings_nb.grid(row=0, column=0, padx=10, pady=10)

# =============================================================================
#       Model settings
# =============================================================================

        # Create row counters
        row_n = 0
        col_n = 0

        ttk.Label(model_frame, text='Fit Window:\n     (nm)', font=NORM_FONT
                  ).grid(row=row_n, column=col_n, rowspan=2, padx=5, pady=5)

        self.widgets['fit_lo'] = tk.DoubleVar(model_frame, value=310)
        make_input(frame=model_frame,
                   text=None,
                   row=row_n, column=col_n,
                   var=self.widgets['fit_lo'],
                   input_type='Spinbox',
                   width=8,
                   sticky='S')
        row_n += 1

        self.widgets['fit_hi'] = tk.DoubleVar(model_frame, value=320)
        make_input(frame=model_frame,
                   text=None,
                   row=row_n, column=col_n,
                   var=self.widgets['fit_hi'],
                   input_type='Spinbox',
                   width=8,
                   sticky='N')
        row_n += 1

        # Set model grid padding
        self.widgets['model_padding'] = tk.DoubleVar(model_frame, value=1.0)
        make_input(frame=model_frame,
                   text='Model Grid\nPadding (nm):',
                   row=row_n, column=col_n,
                   var=self.widgets['model_padding'],
                   input_type='Spinbox',
                   increment=0.1,
                   width=8,
                   sticky='S')
        row_n += 1

        # Set model grid spacing
        self.widgets['model_spacing'] = tk.DoubleVar(model_frame, value=0.01)
        make_input(frame=model_frame,
                   text='Model Grid\nSpacing (nm):',
                   row=row_n, column=col_n,
                   var=self.widgets['model_spacing'],
                   input_type='Spinbox',
                   increment=0.01,
                   width=8,
                   sticky='S')
        row_n += 1

        # Control the interpolation mode
        interp_options = ['cubic', 'cubic', 'linear']
        self.widgets['interp_method'] = tk.StringVar(model_frame,
                                                     value='cubic')
        make_input(frame=model_frame,
                   text='Interpolation\nMethod:',
                   row=row_n, column=col_n,
                   var=self.widgets['interp_method'],
                   input_type='OptionMenu',
                   options=interp_options,
                   width=8,
                   sticky='W')
        row_n += 1

        # New column
        row_n = 0
        col_n += 3

        # Create separator
        ttk.Separator(model_frame, orient='vertical'
                      ).grid(column=col_n-1, row=0, rowspan=7, sticky='ns',
                             padx=10, pady=5)

        # Control whether or not to remove dark spectra
        self.widgets['dark_flag'] = tk.BooleanVar(model_frame, value=True)
        make_input(frame=model_frame,
                   text='Remove Dark\nSpectrum?',
                   row=row_n, column=col_n,
                   var=self.widgets['dark_flag'],
                   input_type='Checkbutton')
        row_n += 1

        # Control whether or not to remove flat spectra
        self.widgets['flat_flag'] = tk.BooleanVar(model_frame, value=True)
        make_input(frame=model_frame,
                   text='Remove Flat\nSpectrum?',
                   row=row_n, column=col_n,
                   var=self.widgets['flat_flag'],
                   input_type='Checkbutton')
        row_n += 1

        # Control whether to correct for stray light
        self.widgets['stray_flag'] = tk.BooleanVar(model_frame, value=True)
        make_input(frame=model_frame,
                   text='Remove Stray\nLight?',
                   row=row_n, column=col_n,
                   var=self.widgets['stray_flag'],
                   input_type='Checkbutton')
        row_n += 1

        # Set the stray light range
        ttk.Label(model_frame, text='Stray Window:\n       (nm)',
                  font=NORM_FONT).grid(row=row_n, column=col_n, rowspan=2,
                                       padx=5, pady=5)
        self.widgets['stray_lo'] = tk.DoubleVar(model_frame, value=280)
        make_input(frame=model_frame,
                   text=None,
                   row=row_n, column=col_n,
                   var=self.widgets['stray_lo'],
                   input_type='Spinbox',
                   width=8,
                   sticky='S')
        row_n += 1

        self.widgets['stray_hi'] = tk.DoubleVar(model_frame, value=290)
        make_input(frame=model_frame,
                   text=None,
                   row=row_n, column=col_n,
                   var=self.widgets['stray_hi'],
                   input_type='Spinbox',
                   width=8,
                   sticky='N')
        row_n += 1

        # Set the despike limits
        self.widgets['despike_flag'] = tk.BooleanVar(model_frame, value=False)
        make_input(frame=model_frame,
                   text='Remove Spikey\nPixels?',
                   row=row_n, column=col_n,
                   var=self.widgets['despike_flag'],
                   input_type='Checkbutton')
        row_n += 1

        # Set the spike limit
        self.widgets['spike_limit'] = tk.DoubleVar(model_frame, value=1000)
        make_input(frame=model_frame,
                   text='Spike Size\n(counts):',
                   row=row_n, column=col_n,
                   var=self.widgets['spike_limit'],
                   input_type='Spinbox',
                   increment=100,
                   width=8,
                   sticky='W')
        row_n += 1

        # New column
        row_n = 0
        col_n += 4

        # Create separator
        ttk.Separator(model_frame, orient='vertical'
                      ).grid(column=col_n-1, row=0, rowspan=7, sticky='ns',
                             padx=10, pady=5)

        # Create entry to select the residual type
        resid_options = ['Percentage', 'Percentage', 'Absolute']

        self.widgets['resid_type'] = tk.StringVar(model_frame,
                                                  value='Percentage')
        make_input(frame=model_frame,
                   text='Residual Display:',
                   row=row_n, column=col_n,
                   var=self.widgets['resid_type'],
                   input_type='OptionMenu',
                   options=resid_options,
                   width=10,
                   sticky='W')
        row_n += 1

        # Control whether to update fit parameter with the last fit values
        self.widgets['update_flag'] = tk.BooleanVar(model_frame, value=True)
        make_input(frame=model_frame,
                   text='Auto-update\nfit parameters?',
                   row=row_n, column=col_n,
                   var=self.widgets['update_flag'],
                   input_type='Checkbutton')
        row_n += 1

        # Control bound of goodness of fit
        self.widgets['resid_limit'] = tk.DoubleVar(model_frame, value=10)
        make_input(frame=model_frame,
                   text='Residual Limit:',
                   row=row_n, column=col_n,
                   var=self.widgets['resid_limit'],
                   input_type='Spinbox',
                   width=8)
        row_n += 1

        # Control intensity limits
        ttk.Label(model_frame, text='Intensity Limits:',
                  font=NORM_FONT).grid(row=row_n, column=col_n, rowspan=2,
                                       padx=5, pady=5)
        self.widgets['lo_int_limit'] = tk.DoubleVar(model_frame, value=0)
        make_input(frame=model_frame,
                   text=None,
                   row=row_n, column=col_n,
                   var=self.widgets['lo_int_limit'],
                   input_type='Spinbox',
                   vals=[-10000, 100000],
                   sticky='S',
                   increment=1000,
                   width=8)
        row_n += 1

        # Control bound of goodness of fit
        self.widgets['hi_int_limit'] = tk.DoubleVar(model_frame, value=70000)
        make_input(frame=model_frame,
                   text=None,
                   row=row_n, column=col_n,
                   var=self.widgets['hi_int_limit'],
                   input_type='Spinbox',
                   vals=[-10000, 100000],
                   sticky='N',
                   increment=1000,
                   width=8)
        row_n += 1

# =============================================================================
#       Spectrometer Settings
# =============================================================================

        # Create a row counter
        row_n = 0
        col_n = 0

        # Create entries for the ILS parameters
        ils_options = ['Params', 'Params', 'File', 'Manual']
        self.widgets['ils_mode'] = tk.StringVar(spect_frame, value='File')
        make_input(frame=spect_frame,
                   text='Generate ILS:',
                   row=row_n, column=col_n,
                   var=self.widgets['ils_mode'],
                   input_type='OptionMenu',
                   options=ils_options,
                   width=10,
                   sticky='W')

        row_n += 1

        # File path to ILS params
        self.widgets['ils_path'] = tk.StringVar(spect_frame)
        make_input(frame=spect_frame,
                   text='ILS Parameters:',
                   row=row_n, column=col_n,
                   var=self.widgets['ils_path'],
                   input_type='Entry',
                   width=40)
        ttk.Button(spect_frame,
                   text="Browse",
                   command=lambda: file_io(single_file=True,
                                           holder=self.widgets['ils_path'])
                   ).grid(row=row_n, column=2, padx=5, pady=5, sticky='W')

        row_n += 1

        # Flat spectrum
        self.widgets['flat_path'] = tk.StringVar(spect_frame)
        make_input(frame=spect_frame,
                   text='Flat Spectrum:',
                   row=row_n, column=col_n,
                   var=self.widgets['flat_path'],
                   input_type='Entry',
                   width=40)
        ttk.Button(spect_frame,
                   text="Browse",
                   command=lambda: file_io(single_file=True,
                                           holder=self.widgets['flat_path'])
                   ).grid(row=row_n, column=2, padx=5, pady=5, sticky='W')

        row_n += 1

        # Path to the wavelength calibration file
        self.widgets['wl_calib'] = tk.StringVar(spect_frame, value='')
        make_input(frame=spect_frame,
                   text='Wavelength\nCalibration:',
                   row=row_n, column=col_n,
                   var=self.widgets['wl_calib'],
                   input_type='Entry',
                   width=40)
        ttk.Button(spect_frame,
                   text="Browse",
                   command=lambda: file_io(single_file=True,
                                           holder=self.widgets['wl_calib'])
                   ).grid(row=row_n, column=2, padx=5, pady=5, sticky='W')

        row_n += 1

        ttk.Separator(spect_frame, orient='horizontal'
                      ).grid(column=0,
                             row=row_n,
                             sticky='ew',
                             padx=10,
                             pady=0,
                             columnspan=3)

        row_n += 1

        # Make a frame to hold ILS settings
        ils_frame = tk.LabelFrame(spect_frame, text='ILS Parameters',
                                  font=LARGE_FONT)
        ils_frame.grid(row=row_n, column=0, padx=10, pady=10, sticky='NW',
                       columnspan=3)

        # Control the spectrometer ILS width
        self.widgets['fwem'] = tk.DoubleVar(ils_frame)
        make_input(frame=ils_frame,
                   text='FWEM:',
                   row=0, column=0,
                   var=self.widgets['fwem'],
                   input_type='Entry',
                   width=12)
        self.widgets['fwem_fit'] = tk.BooleanVar(ils_frame)
        ttk.Checkbutton(ils_frame, text='Fit?',
                        variable=self.widgets['fwem_fit']
                        ).grid(row=0, column=2)

        # Control the spectrometer ILS k value
        self.widgets['k'] = tk.DoubleVar(ils_frame)
        make_input(frame=ils_frame,
                   text='k:',
                   row=1, column=0,
                   var=self.widgets['k'],
                   input_type='Entry',
                   width=12)
        self.widgets['k_fit'] = tk.BooleanVar(ils_frame)
        ttk.Checkbutton(ils_frame, text='Fit?', variable=self.widgets['k_fit']
                        ).grid(row=1, column=2)

        # Control the spectrometer ILS a_w value
        self.widgets['a_w'] = tk.DoubleVar(ils_frame)
        make_input(frame=ils_frame,
                   text='a_w:',
                   row=2, column=0,
                   var=self.widgets['a_w'],
                   input_type='Entry',
                   width=12)
        self.widgets['a_w_fit'] = tk.BooleanVar(ils_frame)
        ttk.Checkbutton(ils_frame, text='Fit?',
                        variable=self.widgets['a_w_fit']).grid(row=2, column=2)

        # Control the spectrometer ILS a_k value
        self.widgets['a_k'] = tk.DoubleVar(ils_frame)
        make_input(frame=ils_frame,
                   text='a_k:',
                   row=3, column=0,
                   var=self.widgets['a_k'],
                   input_type='Entry',
                   width=12)
        self.widgets['a_k_fit'] = tk.BooleanVar(ils_frame)
        ttk.Checkbutton(ils_frame, text='Fit?',
                        variable=self.widgets['a_k_fit']).grid(row=3, column=2)

# =============================================================================
#       Parameter Settings
# =============================================================================

        # Create a row counter
        row_n = 0
        col_n = 0

        # File path to solar reference file
        self.widgets['frs_path'] = tk.StringVar(param_frame)
        make_input(frame=param_frame,
                   text='Fraunhofer File:',
                   row=row_n, column=col_n,
                   var=self.widgets['frs_path'],
                   input_type='Entry',
                   sticky='E',
                   width=40,
                   columnspan=3)
        ttk.Button(param_frame,
                   text="Browse",
                   command=lambda: file_io(single_file=True,
                                           holder=self.widgets['frs_path'])
                   ).grid(row=row_n, column=4, padx=5, pady=5, sticky='W')

        row_n += 1

        ttk.Separator(param_frame, orient='horizontal'
                      ).grid(column=0, row=row_n, columnspan=6, sticky='we',
                             padx=10, pady=5)
        row_n += 1

        # Make a notebook to hold the different parameter types
        param_nb = ttk.Notebook(param_frame)
        abs_page = ttk.Frame(param_nb)
        poly_page = ttk.Frame(param_nb)
        offset_page = ttk.Frame(param_nb)
        shift_page = ttk.Frame(param_nb)
        param_nb.add(abs_page, text='Absorbers')
        param_nb.add(poly_page, text='Polynomial')
        param_nb.add(offset_page, text='Offset')
        param_nb.add(shift_page, text='Shift')
        param_nb.grid(row=row_n, column=0, padx=10, pady=10, columnspan=10)
        row_n += 1

        # Generate the tables to hold the parameters
        self.gastable = ParamTable(abs_page)
        self.gastable.grid(row=0, column=0, padx=5, pady=5, columnspan=6)

        self.polytable = PolyTable(poly_page)
        self.polytable.grid(row=0, column=0, padx=5, pady=5, columnspan=6,
                            sticky='NSEW')

        self.offsettable = PolyTable(offset_page)
        self.offsettable.grid(row=0, column=0, padx=5, pady=5, columnspan=6)

        self.shifttable = PolyTable(shift_page)
        self.shifttable.grid(row=0, column=0, padx=5, pady=5, columnspan=6)

# =============================================================================
# =============================================================================
#   General program functions
# =============================================================================
# =============================================================================

    # Report exceptions in a new window
    def report_callback_exception(self, *args):

        # Get error and format into a single string
        err = traceback.format_exception(*args)
        err_msg = ''.join(map(str, err))

        # Log error
        logging.exception(err_msg)

        # Report error
        tkMessageBox.showerror('Exception', err_msg)

    # Close program on 'x' button
    def handler(self):

        # Turn on stopping flag
        self.stop_flag = False

        # Open save dialouge
        text = 'Would you like to\nsave the settings?'
        message = tkMessageBox.askquestion('Exit',
                                           message=text,
                                           type='yesnocancel')

        if message == 'yes':
            write_config(self)
            self.quit()

        if message == 'no':
            self.quit()

        if message == 'cancel':
            pass


if __name__ == '__main__':

    root = tk.Tk()
    myGUI(root)
    root.mainloop()
