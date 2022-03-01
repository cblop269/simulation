import tkinter as tk
from pathlib import Path
from tkinter.filedialog import askopenfilename
import numpy as np
import astropy.units as u
from astropy.units import cds
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from numba import jit
from mpl_toolkits.axes_grid1 import make_axes_locatable
from backClasses.plotter import Plotter
from backClasses.interferometer import Interferometer
from backClasses.filemanager import FileManager
from backClasses.gridder import Gridder
from backClasses.fouriertransformer import FT

class FrameC(tk.Frame):

    def __init__(self, notebook):#, interferometer):
        super().__init__(notebook)
        self.config(bg="gray26")
        #
        self.plotter = Plotter(2, 1, 9, 6, w=.5, l=0.1, t=0.9)

        #   Frames
        self.canvasD = FigureCanvasTkAgg(self.plotter.figure, self)
        self.canvasD.get_tk_widget().grid(row=0, column=1, rowspan=2, pady=20, padx=20)
        self.frame8 = tk.Frame(self, bg="gray26")
        self.frame8.grid(row=0, column=0, pady=20, padx=10)

        self.set_labels()
        self.set_variables()
        self.set_input_entries()
        self.set_check_buttons()
        self.set_dropdown_lists()

    def set_labels(self):
        #   labels
        self.lbl_settings2 = tk.Label(self.frame8, text="Settings", bg="gray26", fg="white", font=("Courier", 16))
        self.lbl_system_temperature = tk.Label(self.frame8, text="system temperature", bg="gray26", fg="white")
        self.lbl_integration_time = tk.Label(self.frame8, text="integration time", bg="gray26", fg="white")
        self.lbl_bandwidth = tk.Label(self.frame8, text="bandwidth", bg="gray26", fg="white")
        self.lbl_settings2.grid(row=0, column=0, columnspan=2, padx=10, pady=10)
        self.lbl_system_temperature.grid(row=4, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.lbl_integration_time.grid(row=5, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.lbl_bandwidth.grid(row=6, column=0, columnspan=1, padx=10, pady=10, sticky=tk.W)

    def set_variables(self):
        #   Variables
        self.use_noise = tk.BooleanVar()
        self.use_noise.set(True)
        self.system_temperature = tk.StringVar(value='200')
        self.integration_time = tk.StringVar(value='60')
        self.bandwidth = tk.StringVar(value='300')

    def set_input_entries(self):
        #	Input
        self.in_system_temperature = tk.Entry(self.frame8, width=15, justify='right')
        self.in_system_temperature.config(textvariable=self.system_temperature)
        self.in_integration_time = tk.Entry(self.frame8, width=15, justify='right')
        self.in_integration_time.config(textvariable=self.integration_time)
        self.in_bandwidth = tk.Entry(self.frame8, width=15, justify='right')
        self.in_bandwidth.config(textvariable=self.bandwidth)
        self.in_system_temperature.grid(row=4, column=1, columnspan=1, padx=10, pady=10)
        self.in_integration_time.grid(row=5, column=1, columnspan=1, padx=10, pady=10)
        self.in_bandwidth.grid(row=6, column=1, columnspan=1, padx=10, pady=10)

    def set_check_buttons(self):
        #   Check Buttons
        checkbutton_style = tk.ttk.Style()
        checkbutton_style.configure('color.TCheckbutton', foreground='white', background='gray26')
        self.chk_noise = tk.ttk.Checkbutton(self.frame8, text='Add Noise', variable=self.use_noise, style='color.TCheckbutton')
        self.chk_noise.config(onvalue=True, offvalue=False, command=self.checkbuttonnoise)
        self.chk_noise.grid(row=1, column=1, columnspan=2, padx=10, pady=10)

    def set_dropdown_lists(self):
        # Dropdown List
        self.option_sysT = {"K": u.K}
        self.sysT_unit = tk.StringVar()
        self.sysT_unit.set("K")
        self.dl_sysT = tk.OptionMenu(self.frame8, self.sysT_unit, *self.option_sysT.keys())
        self.dl_sysT.config(width=4, justify='right')

        self.option_integr_time = {"sec": u.s, "min": u.min, "hr": u.hour}
        self.integr_time_unit = tk.StringVar()
        self.integr_time_unit.set("sec")
        self.dl_integr_time = tk.OptionMenu(self.frame8, self.integr_time_unit, *self.option_integr_time.keys())
        self.dl_integr_time.config(width=4, justify='right')

        self.option_bandwdith = {"GHz": u.Hz * 1e9, "MHz": u.Hz * 1e6, "KHz": u.Hz * 1e3, "Hz": u.Hz}
        self.bandwdith_unit = tk.StringVar()
        self.bandwdith_unit.set("GHz")
        self.dl_bandwdith = tk.OptionMenu(self.frame8, self.bandwdith_unit, *self.option_bandwdith.keys())
        self.dl_bandwdith.config(width=4, justify='right')
        self.dl_sysT.grid(row=4, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.dl_integr_time.grid(row=5, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)
        self.dl_bandwdith.grid(row=6, column=3, columnspan=1, padx=10, pady=10, sticky=tk.W)

    def checkbuttonnoise(self):
        if self.use_noise.get():
            self.use_noise.set(True)
            self.in_system_temperature.config(state='normal')
            self.in_integration_time.config(state='normal')
            self.in_bandwidth.config(state='normal')
        else:
            self.use_noise.set(False)
            self.in_system_temperature.config(state='disabled')
            self.in_integration_time.config(state='disabled')
            self.in_bandwidth.config(state='disabled')

    def get_noise_inputs(self):
        system_temperature = (float(self.in_system_temperature.get())) * self.option_sysT[self.sysT_unit.get()]
        integration_time = float(self.in_integration_time.get()) * self.option_integr_time[self.integr_time_unit.get()]
        bandwidth = float(self.in_bandwidth.get()) * self.option_bandwdith[self.bandwdith_unit.get()]
        return [system_temperature, integration_time, bandwidth]

    def draw_noise(self, positionUV, noise_value, value):
        self.plotter.ax[0].cla()
        self.plotter.ax[1].cla()
        self.plotter.ax[0].set_title('Visibilities Value')
        self.plotter.ax[0].set_xlabel('distance (λ)')
        self.plotter.ax[0].set_ylabel('value (Jy)')
        self.plotter.ax[1].set_title('Noise Values')
        self.plotter.ax[1].set_xlabel('distance (λ)')
        self.plotter.ax[1].set_ylabel('value (Jy)')
        # get the uv value
        noise_value = np.abs(noise_value)
        value = np.abs(value)
        # get the norm from u and v
        distance = np.linalg.norm(positionUV, axis=0)
        # draw plots
        self.plotter.ax[0].plot(distance[distance.argsort()], value[distance.argsort()], linewidth=0.2)
        self.plotter.ax[1].plot(distance[distance.argsort()], noise_value[distance.argsort()], linewidth=0.2)
        self.canvasD.draw()
