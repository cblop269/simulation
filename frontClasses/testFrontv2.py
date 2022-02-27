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

class Window(tk.Tk):

    def __init__(self):
        super().__init__()
        self.config(bg="gray20")
        # self.minsize(1200, 600)
        #
        self.interferometer = Interferometer('/home/seba/Desktop/alma.C34-2.cfg')
        self.writer = FileManager()
        self.filename_antenna = tk.StringVar()
        self.filename_source = tk.StringVar()
        self.filename_antenna.set(self.get_file_name('/home/seba/Desktop/alma.C34-2.cfg'))
        self.filename_source.set(self.get_file_name('/home/seba/Downloads/cameraman(1).fits'))
        #
        self.create_sub_menu()
        self.create_widgets()

        self.charge_sky_image('/home/seba/Downloads/cameraman(1).fits')
        self.charge_antenna_config('/home/seba/Desktop/alma.C34-2.cfg')
        #
        self.plotter.draw_inputs(self.canvasA, self.interferometer.sky_image, self.interferometer.antenna_pos)

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
        self.frameD = tk.Frame(self.notebook, bg="gray26")
        self.frameE = tk.Frame(self.notebook, bg="gray26")
        self.frameF = tk.Frame(self.notebook, bg="gray26")
        self.frameG = tk.Frame(self.notebook, bg="gray26")
        self.frame_toolbar = tk.Frame(self.frameB, bg="gray26")
        self.notebook.add(self.frameA, text="Input")
        self.notebook.add(self.frameB, text="Output")
        self.notebook.add(self.frameD, text="Noise")
        self.notebook.add(self.frameC, text="Export")
        self.notebook.add(self.frameE, text="Visibilities")
        self.notebook.add(self.frameF, text="Gridded Image")
        self.notebook.add(self.frameG, text="Dirty Image")
        self.notebook.pack(padx=20, pady=20)
        self.notebook.select(self.frameA)
        #
        self.plotter = Plotter()
        self.canvasB = FigureCanvasTkAgg(self.plotter.figureB, self.frame_toolbar)
        self.canvasB.get_tk_widget().pack()
        self.frame_toolbar.grid(row=0, column=1, pady=10, padx=0, sticky=tk.NW)
        self.canvasA = FigureCanvasTkAgg(self.plotter.figureA, self.frameA)
        self.canvasA.get_tk_widget().grid(row=0, column=1, rowspan=2, pady=20, padx=20)
        self.canvasD = FigureCanvasTkAgg(self.plotter.figureD, self.frameD)
        self.canvasD.get_tk_widget().grid(row=0, column=1, rowspan=2, pady=20, padx=20)
        self.toolbar = NavigationToolbar2Tk(self.canvasB, self.frame_toolbar)
        self.toolbar.update()
        # ///////////////////////////////////////////////////////////////////////////////////////////
        self.canvasE = FigureCanvasTkAgg(self.plotter.figureE, self.frameE)
        self.canvasE.get_tk_widget().pack() #grid(row=0, column=1, rowspan=2, pady=20, padx=20)
        self.canvasF = FigureCanvasTkAgg(self.plotter.figureF, self.frameF)
        self.canvasF.get_tk_widget().pack() #grid(row=0, column=1, rowspan=2, pady=20, padx=20)
        self.canvasG = FigureCanvasTkAgg(self.plotter.figureG, self.frameG)
        self.canvasG.get_tk_widget().pack() #grid(row=0, column=1, rowspan=2, pady=20, padx=20)
        self.toolbarE = NavigationToolbar2Tk(self.canvasE, self.frameE)
        self.toolbarF = NavigationToolbar2Tk(self.canvasF, self.frameF)
        self.toolbarG = NavigationToolbar2Tk(self.canvasG, self.frameG)
        self.toolbarE.update()
        self.toolbarF.update()
        self.toolbarG.update()
        # ///////////////////////////////////////////////////////////////////////////////////////////

        #   Frames
        self.frame1 = tk.Frame(self.frameA, bg="gray26")
        self.frame2 = tk.Frame(self.frameA, bg="gray26")
        self.frame3 = tk.Frame(self, bg="gray26")
        self.frame4 = tk.Frame(self.frameC, bg="gray30")
        self.frame5 = tk.Frame(self.frameC, bg="gray30")
        self.frame6 = tk.Frame(self.frameC, bg="gray30")
        self.frame7 = tk.Frame(self.frameC, bg="gray30")
        self.frame8 = tk.Frame(self.frameD, bg="gray26")

        #   labels
        self.label_settings = tk.Label(self.frame1, text='Simulation parameters', font=("Courier", 16), bg="gray26",
                                       fg="white")
        self.label_latitude = tk.Label(self.frame1, text='observatory latitude', bg="gray26", fg="white")
        self.label_declination = tk.Label(self.frame1, text='source declination', bg="gray26", fg="white")
        self.label_sample_frec = tk.Label(self.frame1, text='Sample frequency', bg="gray26", fg="white")
        self.label_sample_n = tk.Label(self.frame1, text='Sampling interval', bg="gray26", fg="white")
        self.label_hour_angle_S = tk.Label(self.frame1, text='Hour Angle Start', bg="gray26", fg="white")
        self.label_hour_angle_E = tk.Label(self.frame1, text='Hour Angle End', bg="gray26", fg="white")
        self.label_algorythm = tk.Label(self.frame1, text='Algorithm', bg="gray26", fg="white")
        self.label_file_antenna = tk.Label(self.frame2, textvariable=self.filename_antenna, bg="gray26", fg="white")
        self.label_file_source = tk.Label(self.frame2, textvariable=self.filename_source, bg="gray26", fg="white")
        self.label_export_ti = tk.Label(self.frame4, text="Export Transform Image", bg="gray30", fg="white")
        self.label_export_uv = tk.Label(self.frame5, text="Export Visibilities", bg="gray30", fg="white")
        self.label_export_gi = tk.Label(self.frame6, text="Export Grided Image", bg="gray30", fg="white")
        self.label_export_di = tk.Label(self.frame7, text="Export Dirty Image", bg="gray30", fg="white")
        self.label_route_ti = tk.Label(self.frame4, text=str(Path.home() / "Downloads"), bg="gray30", fg="white")
        self.label_route_uv = tk.Label(self.frame5, text=str(Path.home() / "Downloads"), bg="gray30", fg="white")
        self.label_route_gi = tk.Label(self.frame6, text=str(Path.home() / "Downloads"), bg="gray30", fg="white")
        self.label_route_di = tk.Label(self.frame7, text=str(Path.home() / "Downloads"), bg="gray30", fg="white")
        self.label_settings2 = tk.Label(self.frame8, text="Settings", bg="gray26", fg="white", font=("Courier", 16))
        self.label_system_temperature = tk.Label(self.frame8, text="system temperature", bg="gray26", fg="white")
        self.label_integration_time = tk.Label(self.frame8, text="integration time", bg="gray26", fg="white")
        self.label_bandwidth = tk.Label(self.frame8, text="bandwidth", bg="gray26", fg="white")

        #   Variables
        latitude_str = tk.StringVar(value='-23.0234')
        rad_dec_str = tk.IntVar(value=-60)  # -40)
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
        self.fft.set(False)
        self.nufft.set(True)
        self.current_algorythm.set(False)
        self.use_noise = tk.BooleanVar()
        self.use_noise.set(True)
        default_export_ti = tk.StringVar()
        default_export_uv = tk.StringVar()
        default_export_gi = tk.StringVar()
        default_export_di = tk.StringVar()
        default_export_ti.set('TransformImage')
        default_export_uv.set('Visibilities')
        default_export_gi.set('GriddedImage')
        default_export_di.set('DirtyImage')
        system_temperature = tk.StringVar(value='200')
        integration_time = tk.StringVar(value='60')
        bandwidth = tk.StringVar(value='300')


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
        self.input_system_temperature = tk.Entry(self.frame8, width=15, justify='right')
        self.input_system_temperature.config(textvariable=system_temperature)
        self.input_integration_time = tk.Entry(self.frame8, width=15, justify='right')
        self.input_integration_time.config(textvariable=integration_time)
        self.input_bandwidth = tk.Entry(self.frame8, width=15, justify='right')
        self.input_bandwidth.config(textvariable=bandwidth)

        #   Check Buttons
        checkbutton_style = tk.ttk.Style()
        checkbutton_style.configure('color.TCheckbutton', foreground='white', background='gray26')
        self.check_fft = tk.ttk.Checkbutton(self.frame1, text='FFT', variable=self.fft, style='color.TCheckbutton')
        self.check_fft.config(onvalue=True, offvalue=False, command=lambda: self.checkbuttonftused(True))
        self.check_nufft = tk.ttk.Checkbutton(self.frame1, text='NUFFT', variable=self.nufft, style='color.TCheckbutton')
        self.check_nufft.config(onvalue=True, offvalue=False, command=lambda: self.checkbuttonftused(False))
        self.check_noise = tk.ttk.Checkbutton(self.frame8, text='Add Noise', variable=self.use_noise, style='color.TCheckbutton')
        self.check_noise.config(onvalue=True, offvalue=False, command=self.checkbuttonnoise)

        # Dropdown List
        self.option_observatory = {" ° ": u.deg, "rad": u.rad}
        self.observatory_unit = tk.StringVar()
        self.observatory_unit.set(" ° ")
        self.optionMenu_obs = tk.OptionMenu(self.frame1, self.observatory_unit, * self.option_observatory.keys())
        self.optionMenu_obs.config(width=4, justify='right')

        self.option_antenna = {" ° ": u.deg, "rad": u.rad}
        self.antenna_unit = tk.StringVar()
        self.antenna_unit.set(" ° ")
        self.optionMenu_ant = tk.OptionMenu(self.frame1, self.antenna_unit, *self.option_antenna.keys())
        self.optionMenu_ant.config(width=4, justify='right')

        self.option_ha_start = {"HA": u.hourangle, "rad": u.rad, " ° ": u.deg}
        self.ha_start_unit = tk.StringVar()
        self.ha_start_unit.set("HA")
        self.optionMenu_ha_s = tk.OptionMenu(self.frame1, self.ha_start_unit, * self.option_ha_start.keys())
        self.optionMenu_ha_s.config(width=4, justify='right')

        self.option_ha_end = {"HA": u.hourangle, "rad": u.rad, " ° ": u.deg}
        self.ha_end_unit = tk.StringVar()
        self.ha_end_unit.set("HA")
        self.optionMenu_ha_e = tk.OptionMenu(self.frame1, self.ha_end_unit, * self.option_ha_end.keys())
        self.optionMenu_ha_e.config(width=4, justify='right')

        self.option_frequency = {"GHz": u.Hz * 1e9, "MHz": u.Hz * 1e6, "KHz": u.Hz * 1e3, "Hz": u.Hz}
        self.frequency_unit = tk.StringVar()
        self.frequency_unit.set("GHz")
        self.optionMenu_freq = tk.OptionMenu(self.frame1, self.frequency_unit, *self.option_frequency.keys())
        self.optionMenu_freq.config(width=4, justify='right')

        self.option_sample_interval = {"sec": 15 * u.deg / 3600, "min": 15 * u.deg / 60, "hr": 15 * u.deg}
        self.sample_interval_unit = tk.StringVar()
        self.sample_interval_unit.set("sec")
        self.optionMenu_sampl = tk.OptionMenu(self.frame1, self.sample_interval_unit, *self.option_sample_interval.keys())
        self.optionMenu_sampl.config(width=4, justify='right')

        self.option_sysT = {"K": u.K}
        self.sysT_unit = tk.StringVar()
        self.sysT_unit.set("K")
        self.optionMenu_sysT = tk.OptionMenu(self.frame8, self.sysT_unit, *self.option_sysT.keys())
        self.optionMenu_sysT.config(width=4, justify='right')

        self.option_integr_time = {"sec": u.s, "min": u.min, "hr": u.hour}
        self.integr_time_unit = tk.StringVar()
        self.integr_time_unit.set("sec")
        self.optionMenu_integr_time = tk.OptionMenu(self.frame8, self.integr_time_unit, *self.option_integr_time.keys())
        self.optionMenu_integr_time.config(width=4, justify='right')

        self.option_bandwdith = {"GHz": u.Hz * 1e9, "MHz": u.Hz * 1e6, "KHz": u.Hz * 1e3, "Hz": u.Hz}
        self.bandwdith_unit = tk.StringVar()
        self.bandwdith_unit.set("GHz")
        self.optionMenu_bandwdith = tk.OptionMenu(self.frame8, self.bandwdith_unit, *self.option_bandwdith.keys())
        self.optionMenu_bandwdith.config(width=4, justify='right')

        #   Buttons
        # open file button source file
        self.import_button_S = tk.Button(self.frame2, text='Input file', width=10, height=1)
        self.import_button_S.config(command=self.charge_sky_image, bg="light blue")
        # open file button observatory file
        self.import_button_O = tk.Button(self.frame2, text='Antenna config.', width=10, height=1)
        self.import_button_O.config(command=self.charge_antenna_config, bg="light blue")
        # run button
        self.run_button = tk.Button(self.frame3, text='Run', width=10)
        self.run_button.config(command=self.run)  # self.draw_graphic_f, state=tk.DISABLED)
        # export button
        self.export_button_ti = tk.Button(self.frame4, text='Export', width=10)
        self.export_button_ti.config(state=tk.DISABLED, command=lambda: self.writer.write_array(
            self.interferometer.fft_image, self.label_route_ti.cget("text"), default_export_ti.get()))
        self.export_button_uv = tk.Button(self.frame5, text='Export', width=10)
        self.export_button_uv.config(state=tk.DISABLED, command=lambda: self.writer.write_visibilities(
            self.interferometer, self.label_route_uv.cget("text"), default_export_uv.get()))
        self.export_button_gi = tk.Button(self.frame6, text='Export', width=10)
        self.export_button_gi.config(state=tk.DISABLED, command=lambda: self.writer.write_array(
            self.interferometer.gridded_vo, self.label_route_gi.cget("text"), default_export_gi.get()))
        self.export_button_di = tk.Button(self.frame7, text='Export', width=10)
        self.export_button_di.config(state=tk.DISABLED, command=lambda: self.writer.write_image(
            self.interferometer.dirty_image, self.label_route_di.cget("text"), default_export_di.get()))

        # placing widgets
        self.frame1.grid(row=1, column=0, pady=20, padx=20, sticky=tk.N)
        self.frame2.grid(row=0, column=0, pady=20, padx=10, sticky=tk.N)
        self.frame3.pack(pady=0, padx=0)
        self.frame4.grid(row=0, column=0, pady=20, padx=10)
        self.frame5.grid(row=0, column=1, pady=20, padx=10)
        self.frame6.grid(row=0, column=2, pady=20, padx=10)
        self.frame7.grid(row=0, column=3, pady=20, padx=10)
        self.frame8.grid(row=0, column=0, pady=20, padx=10)

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
        self.optionMenu_obs.grid(row=1, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.optionMenu_ant.grid(row=2, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.optionMenu_ha_s.grid(row=3, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.optionMenu_ha_e.grid(row=4, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.optionMenu_freq.grid(row=5, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.optionMenu_sampl.grid(row=6, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)

        self.import_button_S.grid(row=0, column=0, padx=20, pady=20)
        self.import_button_O.grid(row=0, column=1, padx=20, pady=20)
        self.label_file_antenna.grid(row=1, column=1, padx=20, pady=5)
        self.label_file_source.grid(row=1, column=0, padx=20, pady=5)
        self.run_button.grid(row=0, column=0, padx=0, pady=0)
        self.label_export_ti.grid(row=0, column=0, pady=20, padx=20)
        self.label_route_ti.grid(row=1, column=0, pady=0, padx=20)
        self.imput_export_ti.grid(row=2, column=0, pady=20, padx=20)
        self.export_button_ti.grid(row=3, column=0, pady=20, padx=20)
        self.label_export_uv.grid(row=0, column=0, pady=20, padx=20)
        self.label_route_uv.grid(row=1, column=0, pady=0, padx=20)
        self.imput_export_uv.grid(row=2, column=0, pady=20, padx=20)
        self.export_button_uv.grid(row=3, column=0, pady=20, padx=20)
        self.label_export_gi.grid(row=0, column=0, pady=20, padx=20)
        self.label_route_gi.grid(row=1, column=0, pady=0, padx=20)
        self.imput_export_gi.grid(row=2, column=0, pady=20, padx=20)
        self.export_button_gi.grid(row=3, column=0, pady=20, padx=20)
        self.label_export_di.grid(row=0, column=0, pady=20, padx=20)
        self.label_route_di.grid(row=1, column=0, pady=0, padx=20)
        self.imput_export_di.grid(row=2, column=0, pady=20, padx=20)
        self.export_button_di.grid(row=3, column=0, pady=20, padx=20)

        self.check_noise.grid(row=1, column=1, columnspan=2, padx=10, pady=10)
        self.input_system_temperature.grid(row=4, column=1, columnspan=1, padx=10, pady=10)
        self.input_integration_time.grid(row=5, column=1, columnspan=1, padx=10, pady=10)
        self.input_bandwidth.grid(row=6, column=1, columnspan=1, padx=10, pady=10)
        self.label_settings2.grid(row=0, column=0, columnspan=2, padx=10, pady=10)
        self.label_system_temperature.grid(row=4, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.label_integration_time.grid(row=5, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.label_bandwidth.grid(row=6, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.optionMenu_sysT.grid(row=4, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.optionMenu_integr_time.grid(row=5, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.optionMenu_bandwdith.grid(row=6, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)

    def checkbuttonftused(self, usefft: bool = None):
        if usefft:
            self.nufft.set(False)
            self.fft.set(True)
        else:
            self.nufft.set(True)
            self.fft.set(False)

    def checkbuttonnoise(self):
        if self.use_noise.get():
            self.input_system_temperature.config(state='normal')
            self.input_integration_time.config(state='normal')
            self.input_bandwidth.config(state='normal')
        else:
            self.input_system_temperature.config(state='disabled')
            self.input_integration_time.config(state='disabled')
            self.input_bandwidth.config(state='disabled')

    def create_sub_menu(self):
        self.my_menu = tk.Menu(self, bg="gray51")
        self.m1 = tk.Menu(self.my_menu)
        self.m1.add_command(label="Open source... ", command=self.run)
        self.my_menu.add_cascade(label="File", menu=self.m1)
        self.config(menu=self.my_menu)

    def run(self):
        # try:

        # get the respective value inputs and his units
        latitude = (float(self.input_lat.get()) - 90) * self.option_antenna[self.antenna_unit.get()]
        declination = float(self.input_rad_dec.get()) * self.option_observatory[self.observatory_unit.get()]
        ha_start = float(self.input_h_angle_S.get()) * self.option_ha_start[self.ha_start_unit.get()]
        ha_end = float(self.input_h_angle_E.get()) * self.option_ha_end[self.ha_end_unit.get()]
        sample_interval = int(self.input_s_number.get()) * self.option_sample_interval[self.sample_interval_unit.get()]
        frequency = float(self.input_frecuency.get()) * self.option_frequency[self.frequency_unit.get()]
        usefft = bool(self.fft.get())

        # run the interferometer
        self.interferometer.run(latitude, declination, ha_start, ha_end, sample_interval, frequency, usefft)
        # if the noise option is being used, the sigma rms value will not be zero.
        if self.use_noise.get():
            # get the noise parameters
            system_temperature = (float(self.input_system_temperature.get())) * self.option_sysT[self.sysT_unit.get()]
            integration_time = float(self.input_integration_time.get()) * self.option_integr_time[self.integr_time_unit.get()]
            bandwidth = float(self.input_bandwidth.get()) * self.option_bandwdith[self.bandwdith_unit.get()]
            # exception in case of a division by zero
            if integration_time == 0 or bandwidth == 0:
                raise ValueError('invalid parameter, integration time and bandwidth must not be 0')
            # call the function of the interferometer
            self.interferometer.get_noise_level(system_temperature, integration_time, bandwidth)

        # fuctions to obtain the dirty image
        ft = FT()
        self.interferometer.add_noise()
        gridder = Gridder(self.interferometer.visibilities, self.interferometer.noise, len(self.interferometer.sky_image), 1, 0)
        self.interferometer.dirty_image = ft.transform_ifft(usefft, gridder.gridded_vo)
        #self.interferometer.inverse_transform_fft(usefft, gridder.gridded_vo)

        # get the inputs to draw plots
        deltau = self.interferometer.visibilities.deltau
        deltav = self.interferometer.visibilities.deltav
        deltax = self.interferometer.visibilities.deltax
        deltay = self.interferometer.visibilities.deltay

        grid_image = gridder.gridded_vo.value
        dirty_image = self.interferometer.dirty_image.value
        visibilities = np.empty([3, len(self.interferometer.visibilities.UVW[0])])
        visibilities[0] = self.interferometer.visibilities.UVW[0]
        visibilities[1] = self.interferometer.visibilities.UVW[1]
        visibilities[2] = np.log(np.abs(self.interferometer.visibilities.uv_value + self.interferometer.noise).value + 1)
        #fft_image = np.log(abs(self.interferometer.fft_image) + 1)
        grid_image = np.log(np.abs(grid_image) + 1)
        dirty_image = np.abs(dirty_image)

        # draw the plot
        self.plotter.draw_plots_results(self.canvasB, visibilities, grid_image,
                                        dirty_image, [deltax, deltay], [deltau, deltav])
        self.canvasE.draw()
        self.canvasF.draw()
        self.canvasG.draw()
        # change to the plot view
        self.notebook.select(self.frameB)

        # except ValueError as e:
        # messagebox.showerror(message='error: "{}"'.format(e))
        # tk.messagebox.showwarning("Warning", "Fill in all the spaces with numerical values")
        self.plotter.draw_noise(self.canvasD, self.interferometer.visibilities.UVW[:2],
                                self.interferometer.noise, self.interferometer.visibilities.uv_value,
                                deltav)
        self.export_button_ti.config(state=tk.NORMAL)
        self.export_button_uv.config(state=tk.NORMAL)
        self.export_button_gi.config(state=tk.NORMAL)
        self.export_button_di.config(state=tk.NORMAL)

    def charge_sky_image(self, route: str = None):
        if route is None:
            route = tk.filedialog.askopenfilename()
            self.interferometer.read_image(route)
            self.filename_source.set(self.get_file_name(route))
            self.plotter.draw_inputs(self.canvasA, self.interferometer.sky_image, self.interferometer.antenna_pos)
        else:
            self.interferometer.read_image(route)
            self.filename_source.set(self.get_file_name(route))

    def charge_antenna_config(self, route: str = None):
        if route is None:
            route = tk.filedialog.askopenfilename()
            self.interferometer.read_antenna_config(route)
            self.filename_antenna.set(self.get_file_name(route))
            self.interferometer.compute_baselines()
            self.plotter.draw_inputs(self.canvasA, self.interferometer.sky_image, self.interferometer.antenna_pos)
        else:
            self.interferometer.read_antenna_config(route)
            self.filename_antenna.set(self.get_file_name(route))

    def get_file_name(self, filepath):
        filepath = filepath.split('/')
        filename = filepath[len(filepath) - 1]
        return filename
