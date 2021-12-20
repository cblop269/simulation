# from tkinter import *
import tkinter as tk
from tkinter.filedialog import askopenfilename
import numpy as np
from astropy.io import fits
from scipy.constants import c
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from numba import jit
from mpl_toolkits.axes_grid1 import make_axes_locatable
from frontClasses.plotter import Plotter
from backClasses.interferometer import Interferometer


class Window(tk.Tk):

    def __init__(self):
        super().__init__()
        self.config(bg="gray20")
        # self.minsize(1200, 600)
        self.interferometer = Interferometer('/home/seba/Desktop/alma.C34-2.cfg')
        self.filename_antenna = tk.StringVar()
        self.filename_antenna.set(self.get_file_name('/home/seba/Desktop/alma.C34-2.cfg'))

        self.filename_source = tk.StringVar()
        self.create_sub_menu()
        self.create_widgets()
        self.charge_sky_image('/home/seba/Downloads/cameraman(1).fits')
        self.charge_antenna_config('/home/seba/Desktop/alma.C34-2.cfg')
        self.plotter.draw_inputs(self.canvasA, self.interferometer.sky_image, self.interferometer.antenna_pos)

        # self.interferometer.read_image('/home/seba/Downloads/cameraman(1).fits')
        #
        # self.filename_source.set(self.get_file_name('/home/seba/Downloads'))

    def create_widgets(self):
        # Style
        style = tk.ttk.Style()
        style.theme_use('default')
        style.configure('TNotebook', background="gray20")
        style.configure('TNotebook.Tab', background="gray20")
        style.map("TNotebook.Tab", background=[("selected", "gray26")])
        #
        self.notebook = tk.ttk.Notebook(self)
        self.frameA = tk.Frame(self.notebook, bg="gray26")
        self.frameB = tk.Frame(self.notebook, bg="gray26")
        self.frameC = tk.Frame(self.notebook, bg="gray26")
        self.frame_toolbar = tk.Frame(self.frameB, bg="gray26")
        self.notebook.add(self.frameA, text="input")
        self.notebook.add(self.frameB, text="output")
        self.notebook.add(self.frameC, text="export")
        self.notebook.pack(padx=20, pady=20)  # place(relx=0.05, rely=0.05, relwidth=0.9, relheight=0.9)
        self.notebook.select(self.frameA)
        #
        self.plotter = Plotter()
        self.canvasB = FigureCanvasTkAgg(self.plotter.figureB, self.frame_toolbar)
        self.canvasB.get_tk_widget().pack()
        self.frame_toolbar.grid(row=0, column=1, pady=10, padx=0, sticky=tk.NW)
        self.canvasA = FigureCanvasTkAgg(self.plotter.figureA, self.frameA)
        self.canvasA.get_tk_widget().grid(row=0, column=1, rowspan=2, pady=20, padx=20)
        # self.canvasB.get_tk_widget().pack()
        self.toolbar = NavigationToolbar2Tk(self.canvasB, self.frame_toolbar)
        self.toolbar.update()
        # ///////////////////////////////////////////////////////////////////////////////////////////

        #   Frames
        self.frame1 = tk.Frame(self.frameA, bg="gray26")
        self.frame2 = tk.Frame(self.frameA, bg="gray26")
        self.frame3 = tk.Frame(self, bg="gray26")
        self.frame4 = tk.Frame(self.frameC, bg="gray30")
        self.frame5 = tk.Frame(self.frameC, bg="gray30")
        self.frame6 = tk.Frame(self.frameC, bg="gray30")
        self.frame7 = tk.Frame(self.frameC, bg="gray30")

        #   labels
        self.label_settings = tk.Label(self.frame1, text='Simulation parameters', font=("Courier", 16), bg="gray26",
                                            fg="white")
        self.label_latitude = tk.Label(self.frame1, text='observatory latitude', bg="gray26", fg="white")
        self.label_declination = tk.Label(self.frame1, text='source declination', bg="gray26", fg="white")
        self.label_sample_frec = tk.Label(self.frame1, text='Sample frequency', bg="gray26", fg="white")
        self.label_sample_n = tk.Label(self.frame1, text='Sampling interval (s)', bg="gray26", fg="white")
        self.label_hour_angle_S = tk.Label(self.frame1, text='Hour Angle Start', bg="gray26", fg="white")
        self.label_hour_angle_E = tk.Label(self.frame1, text='Hour Angle End', bg="gray26", fg="white")
        self.label_algorythm = tk.Label(self.frame1, text='Algorithm', bg="gray26", fg="white")
        self.label_file_antenna = tk.Label(self.frame2, textvariable=self.filename_antenna, bg="gray26", fg="white")
        self.label_file_source = tk.Label(self.frame2, textvariable=self.filename_source, bg="gray26", fg="white")
        self.label_export_ti = tk.Label(self.frame4, text="Export Transform Image", bg="gray30", fg="white")
        self.label_export_uv = tk.Label(self.frame5, text="Export Visibilities", bg="gray30", fg="white")
        self.label_export_gi = tk.Label(self.frame6, text="Export Grided Image", bg="gray30", fg="white")
        self.label_export_di = tk.Label(self.frame7, text="Export Dirty Image", bg="gray30", fg="white")

        #   Variables
        latitude_str = tk.StringVar(value='-23.0234')
        rad_dec_str = tk.IntVar(value=-60) #-40)
        hour_angle_S = tk.IntVar()
        hour_angle_S.set(-6)
        sample_interval = tk.IntVar()
        sample_interval.set(90)
        sample_number = tk.IntVar()
        sample_number.set(600)
        hour_angle_E = tk.IntVar()
        hour_angle_E.set(6)
        self.fft = tk.BooleanVar()
        self.nufft = tk.BooleanVar()
        self.current_algorythm = tk.BooleanVar()
        self.fft.set(True)
        self.nufft.set(False)
        self.current_algorythm.set(True)
        self.source_use = False
        self.obserbatory_use = False
        default_export_ti = tk.StringVar()
        default_export_uv = tk.StringVar()
        default_export_gi = tk.StringVar()
        default_export_di = tk.StringVar()

        #	Input
        self.input_lat = tk.Entry(self.frame1, textvariable=latitude_str, width=12, justify='right')
        self.input_rad_dec = tk.Entry(self.frame1, textvariable=rad_dec_str, width=12, justify='right')
        self.input_h_angle_S = tk.Entry(self.frame1, width=12, justify='right')
        self.input_h_angle_S.config(textvariable=hour_angle_S)
        self.input_frecuency = tk.Entry(self.frame1, width=12, justify='right')
        self.input_frecuency.config(textvariable=sample_interval)
        self.input_s_number = tk.Entry(self.frame1, width=12, justify='right')
        self.input_s_number.config(textvariable=sample_number)
        self.input_h_angle_E = tk.Entry(self.frame1, width=12, justify='right')
        self.input_h_angle_E.config(textvariable=hour_angle_E)
        self.imput_export_ti = tk.Entry(self.frame4, width=15, justify='right')
        self.imput_export_ti.config(textvariable=default_export_ti)
        self.imput_export_uv = tk.Entry(self.frame5, width=15, justify='right')
        self.imput_export_uv.config(textvariable=default_export_uv)
        self.imput_export_gi = tk.Entry(self.frame6, width=15, justify='right')
        self.imput_export_gi.config(textvariable=default_export_gi)
        self.imput_export_di = tk.Entry(self.frame7, width=15, justify='right')
        self.imput_export_di.config(textvariable=default_export_di)

        #   Check Buttons
        checkbutton_style = tk.ttk.Style()
        checkbutton_style.configure('color.TCheckbutton', foreground='white', background='gray26')
        self.check_fft = tk.ttk.Checkbutton(self.frame1, text='FFT', variable=self.fft, style='color.TCheckbutton')
        self.check_fft.config(onvalue=True, offvalue=False, command=lambda: self.checkbuttonvalue(True))
        self.check_nufft = tk.ttk.Checkbutton(self.frame1, text='NUFFT', variable=self.nufft,
                                              style='color.TCheckbutton')
        self.check_nufft.config(onvalue=True, offvalue=False, command=lambda: self.checkbuttonvalue(False))

        #Dropdown List
        option_frequency = ["GHz", "MHz", "KHz", "Hz"]
        self.frequency_unit = tk.StringVar()
        self.frequency_unit.set(option_frequency[0])
        self.option_menu_freq = tk.OptionMenu(self.frame1, self.frequency_unit, *option_frequency)
        self.option_menu_freq.config(width=4, justify='right')
        OPTIONS = ["Jan", "Feb", "Mar"]
        variable = tk.StringVar()
        variable.set(OPTIONS[0])
        w = tk.OptionMenu(self.frame2, variable, *OPTIONS)
        w.config(width=4, justify='right')
        OPTIONS = ["Jan", "Feb", "Mar"]
        variable = tk.StringVar()
        variable.set(OPTIONS[0])
        w = tk.OptionMenu(self.frame2, variable, *OPTIONS)
        w.config(width=4, justify='right')
        OPTIONS = ["Jan", "Feb", "Mar"]
        variable = tk.StringVar()
        variable.set(OPTIONS[0])
        w = tk.OptionMenu(self.frame2, variable, *OPTIONS)
        w.config(width=4, justify='right')

        #   Buttons
        #	open file button source file
        self.import_button_S = tk.Button(self.frame2, text='Input file', width=10, height=1)
        self.import_button_S.config(command=self.charge_sky_image, bg="light blue")
        #	open file button observatory file
        self.import_button_O = tk.Button(self.frame2, text='Antenna config.', width=10, height=1)
        # self.import_button_O.config(command=self.open_file_O, bg="light blue")
        #	run button
        self.run_button = tk.Button(self.frame3, text='Run', width=10)
        self.run_button.config(command=self.run)  # self.draw_graphic_f, state=tk.DISABLED)
        #	export button
        self.export_button_ti = tk.Button(self.frame4, text='Export', width=10)
        self.export_button_ti.config(state=tk.DISABLED)
        self.export_button_uv = tk.Button(self.frame5, text='Export', width=10)
        self.export_button_uv.config(state=tk.DISABLED)
        self.export_button_gi = tk.Button(self.frame6, text='Export', width=10)
        self.export_button_gi.config(state=tk.DISABLED)
        self.export_button_di = tk.Button(self.frame7, text='Export', width=10)
        self.export_button_di.config(state=tk.DISABLED)

        # placing widgets
        self.frame1.grid(row=1, column=0, pady=20, padx=20, sticky=tk.N)
        self.frame2.grid(row=0, column=0, pady=20, padx=10, sticky=tk.N)
        self.frame3.pack(pady=0, padx=0)
        self.frame4.grid(row=0, column=0, pady=20, padx=10)
        self.frame5.grid(row=0, column=1, pady=20, padx=10)
        self.frame6.grid(row=0, column=2, pady=20, padx=10)
        self.frame7.grid(row=0, column=3, pady=20, padx=10)

        self.label_settings.grid(row=0, column=0, columnspan=4, pady=5, padx=0)
        self.label_latitude.grid(row=1, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.label_declination.grid(row=2, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.label_hour_angle_S.grid(row=3, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.label_hour_angle_E.grid(row=4, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.label_sample_frec.grid(row=5, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.label_sample_n.grid(row=6, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.label_algorythm.grid(row=7, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.input_lat.grid(row=1, column=1, columnspan=2, padx=10, pady=10)
        self.input_rad_dec.grid(row=2, column=1, columnspan=2, padx=10, pady=10)
        self.input_h_angle_S.grid(row=3, column=1, columnspan=2, padx=10, pady=10)
        self.input_h_angle_E.grid(row=4, column=1, columnspan=2, padx=10, pady=10)
        self.input_frecuency.grid(row=5, column=1, columnspan=2, padx=10, pady=10)
        self.input_s_number.grid(row=6, column=1, columnspan=2, padx=10, pady=10)
        self.check_fft.grid(row=7, column=1, columnspan=1, padx=10, pady=10)
        self.check_nufft.grid(row=7, column=2, columnspan=1, padx=10, pady=10)
        self.option_menu_freq.grid(row=5, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)

        self.import_button_S.grid(row=0, column=0, padx=20, pady=20)
        self.import_button_O.grid(row=0, column=1, padx=20, pady=20)
        self.label_file_antenna.grid(row=1, column=1, padx=20, pady=5)
        self.label_file_source.grid(row=1, column=0, padx=20, pady=5)
        self.run_button.grid(row=0, column=0, padx=0, pady=0)
        self.label_export_ti.grid(row=0, column=0, pady=20, padx=20)
        self.imput_export_ti.grid(row=1, column=0, pady=20, padx=20)
        self.export_button_ti.grid(row=2, column=0, pady=20, padx=20)
        self.label_export_uv.grid(row=0, column=0, pady=20, padx=20)
        self.imput_export_uv.grid(row=1, column=0, pady=20, padx=20)
        self.export_button_uv.grid(row=2, column=0, pady=20, padx=20)
        self.label_export_gi.grid(row=0, column=0, pady=20, padx=20)
        self.imput_export_gi.grid(row=1, column=0, pady=20, padx=20)
        self.export_button_gi.grid(row=2, column=0, pady=20, padx=20)
        self.label_export_di.grid(row=0, column=0, pady=20, padx=20)
        self.imput_export_di.grid(row=1, column=0, pady=20, padx=20)
        self.export_button_di.grid(row=2, column=0, pady=20, padx=20)


    def checkbuttonvalue(self, usefft: bool = None):
        if usefft:
            self.nufft.set(False)
            self.fft.set(True)
        else:
            self.nufft.set(True)
            self.fft.set(False)

    def create_sub_menu(self):
        self.my_menu = tk.Menu(self, bg="gray51")
        self.m1 = tk.Menu(self.my_menu)
        self.m1.add_command(label="Open source... ", command=self.run)
        self.my_menu.add_cascade(label="File", menu=self.m1)
        self.config(menu=self.my_menu)

    def run(self):
        # try:
        # get the values and transform change type from string to numerical values
        self.get_units()
        latitude_str = float(self.input_lat.get())
        rad_dec_str = float(self.input_rad_dec.get())
        hour_angle_S = float(self.input_h_angle_S.get())
        hour_angle_E = float(self.input_h_angle_E.get())
        sample_number = int(self.input_s_number.get())
        sample_interval = float(self.input_frecuency.get()) * self.get_units()
        usefft = bool(self.fft.get())

        # run the interferometer
        self.interferometer.run(latitude_str, rad_dec_str, hour_angle_S, hour_angle_E,
                                sample_number, sample_interval, usefft)

        m, n = np.shape(self.interferometer.sky_image)
        ulimit = self.interferometer.deltau
        vlimit = self.interferometer.deltav
        xlimit = self.interferometer.deltax
        ylimit = self.interferometer.deltay
        fft_image = np.log(abs(self.interferometer.fft_image) + 1)
        grid_image = np.log(abs(self.interferometer.gridded_vo) + 1)
        dirty_image = abs(self.interferometer.dirty_image)

        self.plotter.draw_plots_results(self.canvasB, fft_image, self.interferometer.visibilities, grid_image,
                                        dirty_image, [xlimit, ylimit], [ulimit, vlimit])

        # change to the plot view
        self.notebook.select(self.frameB)

        # except ValueError as e:
        # messagebox.showerror(message='error: "{}"'.format(e))
        # tk.messagebox.showwarning("Warning", "Fill in all the spaces with numerical values")

    def get_units(self):
        f = self.frequency_unit.get()
        if(f == "GHz"):
            return 1000000000
        if(f == "MHz"):
            return 1000000
        if(f == "KHz"):
            return 1000
        if(f == "Hz"):
            return 1

    def charge_sky_image(self, route: str = None):
        if route is None:
            route = tk.filedialog.askopenfilename()
        self.interferometer.read_image(route)
        self.filename_source.set(self.get_file_name(route))
        # draw sky image
        delta_X = c / (float(self.input_frecuency.get()) * self.interferometer.baseline.max_baseline)
        delta_X = delta_X / 7
        #self.plotter.draw_inputs(self.canvasA, self.interferometer.sky_image, self.interferometer.antenna_pos)

    def charge_antenna_config(self, route: str = None):
        if route is None:
            route = tk.filedialog.askopenfilename()
            self.interferometer.read_antenna_config(route)
            self.filename_antenna.set(self.get_file_name(route))
        # draw sky image
        # self.plotter.draw_inputs(self.canvasA, self.interferometer.sky_image, self.interferometer.antenna_pos)

    def get_file_name(self, filepath):
        filepath = filepath.split('/')
        filename = filepath[len(filepath) - 1]
        return filename
