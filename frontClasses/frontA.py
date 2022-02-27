import tkinter as tk
from pathlib import Path
from tkinter.filedialog import askopenfilename
import numpy as np
import astropy.units as u
from astropy.units import cds
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from numba import jit
from mpl_toolkits.axes_grid1 import make_axes_locatable
from frontClasses.plotter import Plotter
from backClasses.interferometer import Interferometer
from backClasses.filemanager import FileManager
from backClasses.gridder import Gridder
from backClasses.fouriertransformer import FT

class FrameA(tk.Frame):

    def __init__(self, notebook):#, interferometer):
        super().__init__(notebook)
        self.config(bg="gray26")
        #
        #self.interferometer = Interferometer('/home/seba/Desktop/alma.C34-2.cfg')
        self.filename_antenna = tk.StringVar()
        self.filename_source = tk.StringVar()
        self.filename_antenna.set(self.get_file_name('/home/seba/Desktop/alma.C34-2.cfg'))
        self.filename_source.set(self.get_file_name('/home/seba/Downloads/cameraman(1).fits'))
        #

        #self.charge_sky_image(interferometer, '/home/seba/Downloads/cameraman(1).fits')
        #self.charge_antenna_config(interferometer, '/home/seba/Desktop/alma.C34-2.cfg')
        #
        self.plotter = Plotter(1, 2, 9, 5, w=.5, l=0.08, b=0.01)
        self.colorbar = None

        self.canvasA = FigureCanvasTkAgg(self.plotter.figure, self)
        self.canvasA.get_tk_widget().grid(row=0, column=1, rowspan=2, pady=20, padx=20)

        #   Frames
        self.frame1 = tk.Frame(self, bg="gray26")
        self.frame2 = tk.Frame(self, bg="gray26")
        self.frame1.grid(row=1, column=0, pady=20, padx=20, sticky=tk.N)
        self.frame2.grid(row=0, column=0, pady=20, padx=10, sticky=tk.N)

        self.set_labels()
        self.set_variables()
        self.set_input_entries()
        self.set_check_buttons()
        self.set_dropdown_lists()
        self.set_buttons()

    def set_labels(self):
        #   labels
        self.lbl_settings = tk.Label(self.frame1, text='Simulation parameters', font=("Courier", 16), bg="gray26", fg="white")
        self.lbl_latitude = tk.Label(self.frame1, text='observatory latitude', bg="gray26", fg="white")
        self.lbl_declination = tk.Label(self.frame1, text='source declination', bg="gray26", fg="white")
        self.lbl_sample_frec = tk.Label(self.frame1, text='Sample frequency', bg="gray26", fg="white")
        self.lbl_sample_n = tk.Label(self.frame1, text='Sampling interval', bg="gray26", fg="white")
        self.lbl_hour_angle_S = tk.Label(self.frame1, text='Hour Angle Start', bg="gray26", fg="white")
        self.lbl_hour_angle_E = tk.Label(self.frame1, text='Hour Angle End', bg="gray26", fg="white")
        self.lbl_algorythm = tk.Label(self.frame1, text='Algorithm', bg="gray26", fg="white")
        self.lbl_file_antenna = tk.Label(self.frame2, textvariable=self.filename_antenna, bg="gray26", fg="white")
        self.lbl_file_source = tk.Label(self.frame2, textvariable=self.filename_source, bg="gray26", fg="white")
        self.lbl_settings.grid(row=0, column=0, columnspan=4, pady=5, padx=0)
        self.lbl_latitude.grid(row=1, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.lbl_declination.grid(row=2, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.lbl_hour_angle_S.grid(row=3, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.lbl_hour_angle_E.grid(row=4, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.lbl_sample_frec.grid(row=5, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.lbl_sample_n.grid(row=6, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.lbl_algorythm.grid(row=7, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.lbl_file_antenna.grid(row=1, column=1, padx=20, pady=5)
        self.lbl_file_source.grid(row=1, column=0, padx=20, pady=5)

    def set_variables(self):
        #   Variables
        self.latitude_str = tk.StringVar(value='-23.0234')
        self.rad_dec_str = tk.IntVar(value=-60)  # -40)
        self.hour_angle_S = tk.IntVar()
        self.hour_angle_S.set(-6)
        self.hour_angle_E = tk.IntVar()
        self.hour_angle_E.set(6)
        self.frequency_var = tk.IntVar()
        self.frequency_var.set(90)
        self.sample_interval = tk.IntVar()
        self.sample_interval.set(600)
        self.fft = tk.BooleanVar()
        self.nufft = tk.BooleanVar()
        self.current_algorythm = tk.BooleanVar()
        self.fft.set(False)
        self.nufft.set(True)
        self.current_algorythm.set(False)

    def set_input_entries(self):
        #	Input
        self.in_lat = tk.Entry(self.frame1, textvariable=self.latitude_str, width=12, justify='right')
        self.in_rad_dec = tk.Entry(self.frame1, textvariable=self.rad_dec_str, width=12, justify='right')
        self.in_h_angle_S = tk.Entry(self.frame1, width=12, justify='right')
        self.in_h_angle_S.config(textvariable=self.hour_angle_S)
        self.in_h_angle_E = tk.Entry(self.frame1, width=12, justify='right')
        self.in_h_angle_E.config(textvariable=self.hour_angle_E)
        self.in_frequency = tk.Entry(self.frame1, width=12, justify='right')
        self.in_frequency.config(textvariable=self.frequency_var)
        self.in_s_number = tk.Entry(self.frame1, width=12, justify='right')
        self.in_s_number.config(textvariable=self.sample_interval)
        self.in_lat.grid(row=1, column=1, columnspan=2, padx=10, pady=10)
        self.in_rad_dec.grid(row=2, column=1, columnspan=2, padx=10, pady=10)
        self.in_h_angle_S.grid(row=3, column=1, columnspan=2, padx=10, pady=10)
        self.in_h_angle_E.grid(row=4, column=1, columnspan=2, padx=10, pady=10)
        self.in_frequency.grid(row=5, column=1, columnspan=2, padx=10, pady=10)
        self.in_s_number.grid(row=6, column=1, columnspan=2, padx=10, pady=10)

    def set_check_buttons(self):
        #   Check Buttons
        checkbutton_style = tk.ttk.Style()
        checkbutton_style.configure('color.TCheckbutton', foreground='white', background='gray26')
        self.chk_fft = tk.ttk.Checkbutton(self.frame1, text='FFT', variable=self.fft, style='color.TCheckbutton')
        self.chk_fft.config(onvalue=True, offvalue=False, command=lambda: self.checkbuttonftused(True))
        self.chk_nufft = tk.ttk.Checkbutton(self.frame1, text='NUFFT', variable=self.nufft, style='color.TCheckbutton')
        self.chk_nufft.config(onvalue=True, offvalue=False, command=lambda: self.checkbuttonftused(False))
        self.chk_fft.grid(row=7, column=1, columnspan=1, padx=10, pady=10)
        self.chk_nufft.grid(row=7, column=2, columnspan=1, padx=10, pady=10)

    def set_dropdown_lists(self):
        # Dropdown List
        self.option_observatory = {" ° ": u.deg, "rad": u.rad}
        self.observatory_unit = tk.StringVar()
        self.observatory_unit.set(" ° ")
        self.dl_obs = tk.OptionMenu(self.frame1, self.observatory_unit, *self.option_observatory.keys())
        self.dl_obs.config(width=4, justify='right')

        self.option_antenna = {" ° ": u.deg, "rad": u.rad}
        self.antenna_unit = tk.StringVar()
        self.antenna_unit.set(" ° ")
        self.dl_ant = tk.OptionMenu(self.frame1, self.antenna_unit, *self.option_antenna.keys())
        self.dl_ant.config(width=4, justify='right')

        self.option_ha_start = {"HA": u.hourangle, "rad": u.rad, " ° ": u.deg}
        self.ha_start_unit = tk.StringVar()
        self.ha_start_unit.set("HA")
        self.dl_ha_s = tk.OptionMenu(self.frame1, self.ha_start_unit, *self.option_ha_start.keys())
        self.dl_ha_s.config(width=4, justify='right')

        self.option_ha_end = {"HA": u.hourangle, "rad": u.rad, " ° ": u.deg}
        self.ha_end_unit = tk.StringVar()
        self.ha_end_unit.set("HA")
        self.dl_ha_e = tk.OptionMenu(self.frame1, self.ha_end_unit, *self.option_ha_end.keys())
        self.dl_ha_e.config(width=4, justify='right')

        self.option_frequency = {"GHz": u.Hz * 1e9, "MHz": u.Hz * 1e6, "KHz": u.Hz * 1e3, "Hz": u.Hz}
        self.frequency_unit = tk.StringVar()
        self.frequency_unit.set("GHz")
        self.dl_freq = tk.OptionMenu(self.frame1, self.frequency_unit, *self.option_frequency.keys())
        self.dl_freq.config(width=4, justify='right')

        self.option_sample_interval = {"sec": 15 * u.deg / 3600, "min": 15 * u.deg / 60, "hr": 15 * u.deg}
        self.sample_interval_unit = tk.StringVar()
        self.sample_interval_unit.set("sec")
        self.dl_sampl = tk.OptionMenu(self.frame1, self.sample_interval_unit, *self.option_sample_interval.keys())
        self.dl_sampl.config(width=4, justify='right')

        self.dl_obs.grid(row=1, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.dl_ant.grid(row=2, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.dl_ha_s.grid(row=3, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.dl_ha_e.grid(row=4, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.dl_freq.grid(row=5, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.dl_sampl.grid(row=6, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)

    def set_buttons(self):
        #   Buttons
        # open file button source file
        self.btn_source = tk.Button(self.frame2, text='Input file', width=10, height=1)
        self.btn_source.config(command=self.charge_sky_image, bg="light blue")
        # open file button observatory file
        self.btn_observatory = tk.Button(self.frame2, text='Antenna config.', width=10, height=1)
        self.btn_observatory.config(command=self.charge_antenna_config, bg="light blue")

        # placing widgets
        self.btn_source.grid(row=0, column=0, padx=20, pady=20)
        self.btn_observatory.grid(row=0, column=1, padx=20, pady=20)

    def charge_sky_image(self, interferometer: Interferometer = None, route: str = None):
        print('Funco ', route)
        if route is None:
            route = tk.filedialog.askopenfilename()
            interferometer.read_image(route)
            self.filename_source.set(self.get_file_name(route))
            self.plotter.draw_inputs(self.canvasA, interferometer.sky_image, interferometer.antenna_pos)
        else:
            interferometer.read_image(route)
            self.filename_source.set(self.get_file_name(route))

    def charge_antenna_config(self, interferometer: Interferometer = None, route: str = None):
        print('Funco ', route)
        if route is None:
            route = tk.filedialog.askopenfilename()
            interferometer.read_antenna_config(route)
            self.filename_antenna.set(self.get_file_name(route))
            interferometer.compute_baselines()
            self.plotter.draw_inputs(self.canvasA, interferometer.sky_image, interferometer.antenna_pos)
        else:
            interferometer.read_antenna_config(route)
            self.filename_antenna.set(self.get_file_name(route))

    def get_file_name(self, filepath):
        filepath = filepath.split('/')
        filename = filepath[len(filepath) - 1]
        return filename

    def checkbuttonftused(self, usefft: bool = None):
        if usefft:
            self.nufft.set(False)
            self.fft.set(True)
        else:
            self.nufft.set(True)
            self.fft.set(False)

    def get_inputs(self):
        # get the respective value inputs and his units
        latitude = (float(self.in_lat.get()) - 90) * self.option_antenna[self.antenna_unit.get()]
        declination = float(self.in_rad_dec.get()) * self.option_observatory[self.observatory_unit.get()]
        ha_start = float(self.in_h_angle_S.get()) * self.option_ha_start[self.ha_start_unit.get()]
        ha_end = float(self.in_h_angle_E.get()) * self.option_ha_end[self.ha_end_unit.get()]
        sample_interval = int(self.in_s_number.get()) * self.option_sample_interval[self.sample_interval_unit.get()]
        frequency = float(self.in_frequency.get()) * self.option_frequency[self.frequency_unit.get()]
        usefft = bool(self.fft.get())

        return [latitude, declination, ha_start, ha_end, sample_interval, frequency, usefft]

    def draw_inputs(self, image, antenna_pos):
        if self.colorbar is not None:
            self.colorbar.remove()
        # Sky Image
        im = self.plotter.draw_image(self.plotter.ax[0], image, 'Sky Image', None, None)
        self.colorbar = self.plotter.figure.colorbar(im, ax=self.plotter.ax[0], orientation='horizontal', pad=0.1)
        # Antenna Configuration
        self.plotter.draw_scatter(self.plotter.ax[1], antenna_pos[:, 0], antenna_pos[:, 1], None, 'Antenna Configuration',
                                  'x (m)', 'y (m)', None, None, 1)
        self.canvasA.draw()
