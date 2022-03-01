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

class FrameB(tk.Frame):

    def __init__(self, notebook):#, interferometer):
        super().__init__(notebook)
        self.config(bg="gray26")

        self.plotter = Plotter(1, 3, 14, 5, w=.5, l=0.08, b=0.01)
        self.colorbar_vi = None
        self.colorbar_gi = None
        self.colorbar_di = None

        #   Frames
        self.frame_toolbar = tk.Frame(self, bg="gray26")
        self.frame_toolbar.grid(row=0, column=1, pady=10, padx=0, sticky=tk.NW)
        #
        self.canvasB = FigureCanvasTkAgg(self.plotter.figure, self.frame_toolbar)
        self.canvasB.get_tk_widget().pack()

        self.toolbar = NavigationToolbar2Tk(self.canvasB, self.frame_toolbar)
        self.toolbar.update()

    def draw_plots_results(self, visibilities, grid_image, dirty_image, xy_delta, uv_delta):
        # get deltas
        deltax = xy_delta[0]
        deltay = xy_delta[1]
        deltau = uv_delta[0]
        deltav = uv_delta[1]

        # clear colorbars
        if self.colorbar_vi is not None and self.colorbar_gi is not None and self.colorbar_di is not None:
            self.colorbar_vi.remove()
            self.colorbar_gi.remove()
            self.colorbar_di.remove()

        # Visibilities
        horizontal_values = -visibilities[0]
        vertical_values = visibilities[1]
        uv_values = visibilities[2]
        limit = max(np.max(np.abs(visibilities[0])), np.max(np.abs(visibilities[1])))
        im2 = self.plotter.draw_scatter(self.plotter.ax[0], horizontal_values, vertical_values, uv_values,
                                        'Visibilities', 'u (λ)', 'v (λ)', limit, limit, 0.01)

        # Grid Image
        m, n = grid_image.shape
        limit_hor = deltau * n / 2
        limit_vert = deltav * m / 2
        im3 = self.plotter.draw_image(self.plotter.ax[1], grid_image, 'Gridded Image', 'v (λ)', 'u (λ)', limit_hor,
                                      limit_vert)

        # Dirty Image
        m, n = dirty_image.shape
        limit_hor = deltax * n / 2
        limit_vert = deltay * m / 2
        im4 = self.plotter.draw_image(self.plotter.ax[2], dirty_image, 'Dirty Image', 'y (v)', 'x (v)', limit_hor,
                                      limit_vert)

        # draw colorbars
        self.colorbar_vi = self.plotter.figure.colorbar(im2, ax=self.plotter.ax[0], orientation='horizontal', pad=0.15)
        self.colorbar_gi = self.plotter.figure.colorbar(im3, ax=self.plotter.ax[1], orientation='horizontal', pad=0.15)
        self.colorbar_di = self.plotter.figure.colorbar(im4, ax=self.plotter.ax[2], orientation='horizontal', pad=0.15)
        self.colorbar_vi.set_label('value applying log(|uv value| + 1)')
        self.colorbar_gi.set_label('value applying log(|uv value| + 1)')
        self.colorbar_di.set_label('value applying log(|uv value| + 1)')

        self.canvasB.draw()
