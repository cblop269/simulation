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
from frontClasses.frontA import FrameA
from frontClasses.frontB import FrameB
from frontClasses.frontC import FrameC
from frontClasses.frontD import FrameD
from frontClasses.frontE import FrameE

class Window(tk.Tk):

    def __init__(self):
        super().__init__()
        self.config(bg="gray20")
        # self.minsize(1200, 600)
        #
        self.interferometer = Interferometer('/home/seba/Desktop/alma.C34-2.cfg')
        '''self.filename_antenna = tk.StringVar()
        self.filename_source = tk.StringVar()
        self.filename_antenna.set(self.get_file_name('/home/seba/Desktop/alma.C34-2.cfg'))
        self.filename_source.set(self.get_file_name('/home/seba/Downloads/cameraman(1).fits'))'''
        #
        self.create_sub_menu()
        self.create_widgets()

        '''self.charge_sky_image('/home/seba/Downloads/cameraman(1).fits')
        self.charge_antenna_config('/home/seba/Desktop/alma.C34-2.cfg')
        #
        self.plotter.draw_inputs(self.canvasA, self.interferometer.sky_image, self.interferometer.antenna_pos)'''
        self.frameA.draw_inputs(self.interferometer.sky_image, self.interferometer.antenna_pos)

    def create_widgets(self):
        # Style
        style = tk.ttk.Style()
        style.theme_use('default')
        style.configure('TNotebook', background="gray20")
        style.configure('TNotebook.Tab', background="gray20")
        style.map("TNotebook.Tab", background=[("selected", "gray26")])
        #
        self.notebook = tk.ttk.Notebook(self)
        self.frameA = FrameA(self.notebook)
        self.frameA.charge_antenna_config(self.interferometer, '/home/seba/Desktop/alma.C34-2.cfg')
        self.frameA.charge_sky_image(self.interferometer, '/home/seba/Downloads/cameraman(1).fits')
        self.frameB = FrameB(self.notebook)
        self.frameNoise = FrameC(self.notebook)
        self.frameExport = FrameD(self.notebook)
        self.frameE = FrameE(self.notebook)
        self.frameF = FrameE(self.notebook)
        self.frameG = FrameE(self.notebook)
        self.notebook.add(self.frameA, text="Input")
        self.notebook.add(self.frameB, text="Output")
        self.notebook.add(self.frameNoise, text="Noise")
        self.notebook.add(self.frameExport, text="Export")
        self.notebook.add(self.frameE, text="Visibilities")
        self.notebook.add(self.frameF, text="Gridded Image")
        self.notebook.add(self.frameG, text="Dirty Image")
        self.notebook.pack(padx=20, pady=20)
        self.notebook.select(self.frameA)
        #
        # ///////////////////////////////////////////////////////////////////////////////////////////
        #   Frames
        self.frame3 = tk.Frame(self, bg="gray26")

        #   Buttons
        # run button
        self.run_button = tk.Button(self.frame3, text='Run', width=10)
        self.run_button.config(command=self.run)  # self.draw_graphic_f, state=tk.DISABLED)

        # placing widgets
        self.frame3.pack(pady=0, padx=0)
        self.run_button.grid(row=0, column=0, padx=0, pady=0)


    def create_sub_menu(self):
        self.my_menu = tk.Menu(self, bg="gray51")
        self.m1 = tk.Menu(self.my_menu)
        self.m1.add_command(label="Open source... ", command=self.run)
        self.my_menu.add_cascade(label="File", menu=self.m1)
        self.config(menu=self.my_menu)

    def run(self):
        # try:

        # get the respective value inputs and his units

        latitude, declination, ha_start, ha_end, sample_interval, frequency, usefft = self.frameA.get_inputs()

        # run the interferometer
        self.interferometer.run(latitude, declination, ha_start, ha_end, sample_interval, frequency, usefft)
        # if the noise option is being used, the sigma rms value will not be zero.
        if self.frameNoise.use_noise.get():
            # get the noise parameters
            system_temperature, integration_time, bandwidth = self.frameNoise.get_noise_inputs()
            # exception in case of a division by zero
            if integration_time == 0 or bandwidth == 0:
                raise ValueError('invalid parameter, integration time and bandwidth must not be 0')
            # call the function of the interferometer
            self.interferometer.get_noise_level(system_temperature, integration_time, bandwidth)
        else:
            self.interferometer.sigma = 0 * cds.Jy

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
        self.frameB.draw_plots_results(visibilities, grid_image, dirty_image, [deltax, deltay],
                                       [deltau, deltav])
        self.frameE.draw_plots_results(option=1, visibilities=visibilities)
        self.frameF.draw_plots_results(option=2, grid_image=grid_image, uv_delta=[deltau, deltav])
        self.frameG.draw_plots_results(option=3, dirty_image=dirty_image, xy_delta=[deltax, deltay])

        #self.canvasE.draw()
        #self.canvasF.draw()
        #self.canvasG.draw()
        # change to the plot view
        self.notebook.select(self.frameB)

        # except ValueError as e:
        # messagebox.showerror(message='error: "{}"'.format(e))
        # tk.messagebox.showwarning("Warning", "Fill in all the spaces with numerical values")
        self.frameNoise.draw_noise(self.interferometer.visibilities.UVW[:2], self.interferometer.noise,
                                   self.interferometer.visibilities.uv_value)
        self.frameExport.btn_uv.config(state=tk.NORMAL)
        self.frameExport.btn_gi.config(state=tk.NORMAL)
        self.frameExport.btn_di.config(state=tk.NORMAL)
