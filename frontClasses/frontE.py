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

class FrameE(tk.Frame):

    def __init__(self, notebook):#, interferometer):
        super().__init__(notebook)
        self.config(bg="gray26")
        #
        self.plotter = Plotter(1, 1, 7, 6)
        self.colorbar = None

        self.canvas = FigureCanvasTkAgg(self.plotter.figure, self)
        self.canvas.get_tk_widget().pack()  # grid(row=0, column=1, rowspan=2, pady=20, padx=20)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.update()

    def draw_plots_results(self, option, visibilities: np.ndarray = None, grid_image: np.ndarray = None,
                           dirty_image: np.ndarray = None, xy_delta: np.ndarray = None, uv_delta: np.ndarray = None):


        im = None

        # clear colorbars
        if self.colorbar is not None:
            self.colorbar.remove()

        if option == 1:
            horizontal_values = -visibilities[0]
            vertical_values = visibilities[1]
            uv_values = visibilities[2]
            limit = max(np.max(np.abs(visibilities[0])), np.max(np.abs(visibilities[1])))
            im = self.plotter.draw_scatter(self.plotter.ax, horizontal_values, vertical_values, uv_values, '', 'u (λ)',
                                           'v (λ)', limit, limit, 2)
        elif option == 2:
            # get deltas
            deltau = uv_delta[0]
            deltav = uv_delta[1]
            m, n = grid_image.shape
            limit_hor = deltau * n / 2
            limit_vert = deltav * m / 2
            im = self.plotter.draw_image(self.plotter.ax, grid_image, '', 'v (λ)', 'u (λ)', limit_hor, limit_vert)

        elif option == 3:
            # get deltas
            deltax = xy_delta[0]
            deltay = xy_delta[1]
            m, n = dirty_image.shape
            limit_hor = deltax * n / 2
            limit_vert = deltay * m / 2
            im = self.plotter.draw_image(self.plotter.ax, dirty_image, '', 'y (v)', 'x (v)', limit_hor, limit_vert)

        # draw colorbars
        self.colorbar = self.plotter.figure.colorbar(im, ax=self.plotter.ax, fraction=0.046, pad=0.04)
        self.colorbar.set_label('value applying log(|uv value| + 1)')
        self.canvas.draw()

