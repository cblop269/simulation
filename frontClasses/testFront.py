# from tkinter import *
import tkinter as tk
from tkinter.filedialog import askopenfilename
import numpy as np
from astropy.io import fits

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from numba import jit
from mpl_toolkits.axes_grid1 import make_axes_locatable
from frontClasses.plotter import Plotter
from backClasses.interferometer import Interferometer


class Window(tk.Tk):

    def __init__(self):
        super().__init__()
        self.interferometer = Interferometer('/home/seba/Desktop/alma.C34-2.cfg')
        self.filename_antenna = tk.StringVar()
        self.filename_antenna.set(self.get_file_name('/home/seba/Desktop/alma.C34-2.cfg'))
        self.interferometer.read_image('/home/seba/Downloads/cameraman(1).fits')
        self.filename_source = tk.StringVar()
        self.filename_source.set(self.get_file_name('/home/seba/Downloads'))
        self.create_sub_menu()
        self.create_widgets()

    def create_widgets(self):
        #
        self.notebook = tk.ttk.Notebook(self)
        self.frameA = tk.Frame(self.notebook)
        self.notebook.add(self.frameA, text="Results")
        self.frameB = tk.Frame(self.notebook)
        self.notebook.add(self.frameB, text="Config")
        self.notebook.pack(padx=20, pady=20)  # place(relx=0.05, rely=0.05, relwidth=0.9, relheight=0.9)
        self.notebook.select(self.frameB)
        #
        self.plotter = Plotter()

        self.canvas = FigureCanvasTkAgg(self.plotter.figure, self.frameA)
        self.canvas.get_tk_widget().pack()
        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.update()
        # ///////////////////////////////////////////////////////////////////////////////////////////
        #   Frames
        self.frame1 = tk.Frame(self.frameB, highlightbackground="black", highlightthickness=1)
        self.frame2 = tk.Frame(self.frameB, highlightbackground="black", highlightthickness=0.5)
        self.frame3 = tk.Frame(self.frameB)
        #   labels
        self.label_configuration = tk.Label(self.frame1, text='Configuration', font=("Courier", 16))
        self.label_latitude = tk.Label(self.frame1, text='observatory latitude')
        self.label_declination = tk.Label(self.frame1, text='source declination')
        self.label_sample_frec = tk.Label(self.frame1, text='Sample frequency')
        self.label_sample_n = tk.Label(self.frame1, text='Dt')
        self.label_hour_angle_S = tk.Label(self.frame1, text='Hour Angle Start')
        self.label_hour_angle_E = tk.Label(self.frame1, text='Hour Angle End')
        self.label_algorythm = tk.Label(self.frame1, text='Algorythm')
        self.label_file_antenna = tk.Label(self.frame2, textvariable=self.filename_antenna)
        self.label_file_source = tk.Label(self.frame2, textvariable=self.filename_source)
        #   Variables
        latitude_str = tk.IntVar(value=-23.0234)
        rad_dec_str = tk.IntVar(value=-40)

        hour_angle_S = tk.IntVar()
        hour_angle_S.set(-6)
        sample_interval = tk.IntVar()
        sample_interval.set(90e9)
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

        #	Input
        self.input_lat = tk.Entry(self.frame1, textvariable=latitude_str, width=15, justify='right')
        self.input_rad_dec = tk.Entry(self.frame1, textvariable=rad_dec_str, width=15, justify='right')
        self.input_h_angle_S = tk.Entry(self.frame1, width=15, justify='right')
        self.input_h_angle_S.config(textvariable=hour_angle_S)
        self.input_frecuency = tk.Entry(self.frame1, width=15, justify='right')
        self.input_frecuency.config(textvariable=sample_interval)
        self.input_s_number = tk.Entry(self.frame1, width=15, justify='right')
        self.input_s_number.config(textvariable=sample_number)
        self.input_h_angle_E = tk.Entry(self.frame1, width=15, justify='right')
        self.input_h_angle_E.config(textvariable=hour_angle_E)

        #   Check Buttons
        self.check_fft = tk.ttk.Checkbutton(self.frame1, text='FFT', variable=self.fft)
        self.check_fft.config(onvalue=True, offvalue=False)
        # self.check_fft.config(command=self.unCheck)
        self.check_nufft = tk.ttk.Checkbutton(self.frame1, text='NUFFT', variable=self.nufft)
        self.check_nufft.config(onvalue=True, offvalue=False)
        # self.check_nufft.config(command=self.unCheck)

        #   Buttons
        #	open file button source file
        self.import_button_S = tk.Button(self.frame2, text='Source file', width=10)
        self.import_button_S.config(command=self.charge_sky_image, bg="light blue")
        #	open file button observatory file
        self.import_button_O = tk.Button(self.frame2, text='Observatory file', width=10)
        # self.import_button_O.config(command=self.open_file_O, bg="light blue")
        #	run button
        self.run_button = tk.Button(self.frame3, text='Run', width=10)
        self.run_button.config(command=self.run)  # self.draw_graphic_f, state=tk.DISABLED)
        #	export button
        self.export_button = tk.Button(self.frame3, text='Export', width=10)
        self.export_button.config(state=tk.DISABLED)

        # placing widgets
        self.frame1.place(relx=0.02, rely=0.02, relwidth=0.46, relheight=0.95)
        self.frame2.place(relx=0.52, rely=0.02, relwidth=0.46, relheight=0.68)
        self.frame3.place(relx=0.52, rely=0.7, relwidth=0.46, relheight=0.28)
        self.label_configuration.grid(row=0, column=1, columnspan=2, pady=5, padx=0)
        self.label_latitude.grid(row=1, column=0, columnspan=1, padx=5, sticky=tk.W)
        self.label_declination.grid(row=1, column=2, pady=5, sticky=tk.W)
        self.input_lat.grid(row=1, column=1, pady=5, padx=2)
        self.input_rad_dec.grid(row=1, column=3, pady=5, padx=2)
        self.label_sample_frec.grid(row=3, column=0, pady=5, padx=5, columnspan=2, sticky=tk.W)
        self.label_sample_n.grid(row=3, column=2, pady=10, sticky=tk.W)
        self.label_hour_angle_S.grid(row=2, column=0, pady=5, padx=5, columnspan=2, sticky=tk.W)
        self.label_hour_angle_E.grid(row=2, column=2, pady=5, sticky=tk.W)
        self.input_h_angle_S.grid(row=2, column=1, pady=5, padx=5)
        self.input_frecuency.grid(row=3, column=1, pady=5, padx=5)
        self.input_s_number.grid(row=3, column=3, pady=5, padx=5)
        self.input_h_angle_E.grid(row=2, column=3, pady=5, padx=5)
        self.label_algorythm.grid(row=4, column=0, pady=10)
        self.check_fft.grid(row=4, column=1, pady=10, padx=5)
        self.check_nufft.grid(row=4, column=2, pady=10, padx=10)
        self.import_button_S.place(relx=0.12, rely=0.1, relwidth=0.33, relheight=0.3)
        self.import_button_O.place(relx=0.55, rely=0.1, relwidth=0.33, relheight=0.3)
        self.label_file_antenna.place(relx=0.55, rely=0.5)
        self.label_file_source.place(relx=0.12, rely=0.5)
        self.run_button.place(relx=0.25, rely=0.3, relwidth=0.2, relheight=0.4)
        self.export_button.place(relx=0.55, rely=0.3, relwidth=0.2, relheight=0.4)

    def create_sub_menu(self):
        self.my_menu = tk.Menu(self)
        self.m1 = tk.Menu(self.my_menu)
        self.m1.add_command(label="Open source... ", command=self.run)
        self.my_menu.add_cascade(label="File", menu=self.m1)
        self.config(menu=self.my_menu)

    def run(self):
        #try:
        # get the values and transform change type from string to numerical values
        latitude_str = float(self.input_lat.get())
        rad_dec_str = float(self.input_rad_dec.get())
        hour_angle_S = float(self.input_h_angle_S.get())
        hour_angle_E = float(self.input_h_angle_E.get())
        sample_number = int(self.input_s_number.get())
        sample_interval = float(self.input_frecuency.get())

        # run the interferometer
        self.interferometer.run(latitude_str, rad_dec_str, hour_angle_S, hour_angle_E,
                                sample_number, sample_interval)
        # draw sky image
        delta_X = 1 / (2 * self.interferometer.visibilities.max_uv_coordinate)
        delta_X = delta_X / 7
        x_limit = (1 / 2) * delta_X
        y_limit = (1 / 2) * delta_X
        self.plotter.draw_image(1, self.interferometer.sky_image, self.canvas, x_limit, y_limit,
                                "Sky Image")

        # draw transform image
        m, n = np.shape(self.interferometer.fft_image)
        u_limit = (n / 2) * self.interferometer.delta_U
        v_limit = (m / 2) * self.interferometer.delta_U
        value = np.log(abs(self.interferometer.fft_image) + 1)
        self.plotter.draw_image(2, value, self.canvas, -u_limit, v_limit)

        # draw uv positions
        self.plotter.draw_uv_positions(self.interferometer.visibilities, self.canvas)

        # change to the plot view
        self.notebook.select(self.frameA)

        #except ValueError as e:
            # messagebox.showerror(message='error: "{}"'.format(e))
            #tk.messagebox.showwarning("Warning", "Fill in all the spaces with numerical values")

    def print_parameters(self):
        print(self.input_lat.get())
        print(self.input_rad_dec.get())
        print(self.input_h_angle_S.get())
        print(self.input_frecuency.get())
        print(self.input_s_number.get())
        print(self.input_h_angle_E.get())

    def charge_sky_image(self):
        route = tk.filedialog.askopenfilename()
        self.interferometer.read_image(route)
        self.filename_source.set(self.get_file_name(route))

    def get_file_name(self, filepath):
        filepath = filepath.split('/')
        filename = filepath[len(filepath) - 1]
        return filename
