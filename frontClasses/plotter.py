import numpy as np
import matplotlib.colors as colors
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt


class Plotter:
    def __init__(self):
        self.figureA, self.axA = plt.subplots(1, 2, figsize=(9, 5))
        self.figureB, self.axB = plt.subplots(1, 4, figsize=(14, 5))
        # self.figure.subplots_adjust(hspace=1.5)
        self.figureA.subplots_adjust(wspace=.5, left=0.08, bottom=0.01)
        self.figureB.subplots_adjust(wspace=.5, left=0.08, bottom=0.01)
        # define colorbars
        self.cb_si = None
        self.cb_ti = None
        self.cb_vi = None
        self.cb_gi = None
        self.cb_di = None

    def draw_inputs(self, canvas, image, antenna_pos):
        if self.cb_si is not None:
            self.cb_si.remove()

        # Sky Image
        im = self.draw_image(self.figureA, self.axA[0], image, 'Sky Image', 'x (pixel)', 'y (pixel)')
        self.cb_si = self.figureA.colorbar(im, ax=self.axA[0], orientation='horizontal', pad=0.2)
        # Antenna Configuration
        self.draw_scatter(self.figureA, self.axA[1], antenna_pos[:, 0], antenna_pos[:, 1], None,
                          'Antenna Configuration', 'x (m)', 'y (m)')
        canvas.draw()

    def draw_plots_results(self, canvas, fft_image, visibilities, grid_image, dirty_image, xy_delta, uv_delta):
        # get deltas
        deltax = xy_delta[0]
        deltay = xy_delta[1]
        deltau = uv_delta[0]
        deltav = uv_delta[1]

        # clear colorbars
        if self.cb_ti is not None:
            self.cb_ti.remove()
            self.cb_vi.remove()
            self.cb_gi.remove()
            self.cb_di.remove()

        # FFT Image
        m, n = fft_image.shape
        limit_hor = deltau * n / 2
        limit_vert = deltav * m / 2
        im1 = self.draw_image(self.figureB, self.axB[0], fft_image, 'Transform Image', 'u (λ)', 'v (λ)', limit_hor,
                        limit_vert)

        # Visibilities
        horizontal_values = deltau * visibilities.UVW[0]
        vertical_values = deltav * visibilities.UVW[1]
        limit_hor = visibilities.max_uv_coordinate * deltau
        limit_vert = visibilities.max_uv_coordinate * deltav
        im2 = self.draw_scatter(self.figureB, self.axB[1], horizontal_values, vertical_values,
                          np.log(abs(visibilities.value) + 1), 'Visibilities', 'u (λ)', 'v (λ)', limit_hor, limit_vert)

        # Grid Image
        m, n = grid_image.shape
        limit_hor = deltau * n / 2
        limit_vert = deltav * m / 2
        im3 = self.draw_image(self.figureB, self.axB[2], grid_image, 'Gridded Image', 'v (λ)', 'u (λ)', limit_hor, limit_vert)

        # Dirty Image
        m, n = dirty_image.shape
        limit_hor = deltax * n / 2
        limit_vert = deltay * m / 2
        im4 = self.draw_image(self.figureB, self.axB[3], dirty_image, 'Dirty Image', 'y (v)', 'x (v)', limit_hor, limit_vert)

        # draw colorbars
        self.cb_ti = self.figureB.colorbar(im1, ax=self.axB[0], orientation='horizontal', pad=0.2)
        self.cb_vi = self.figureB.colorbar(im2, ax=self.axB[1], orientation='horizontal', pad=0.2)
        self.cb_gi = self.figureB.colorbar(im3, ax=self.axB[2], orientation='horizontal', pad=0.2)
        self.cb_di = self.figureB.colorbar(im4, ax=self.axB[3], orientation='horizontal', pad=0.2)

        canvas.draw()

    def draw_image(self, figure, subplot, image: np.ndarray, title: str = None, horiz_label: str = None,
                   vert_label: str = None, h_limit: np.float64 = None, v_limit: np.float64 = None):
        subplot.cla()
        subplot.set_title(title)
        subplot.set_xlabel(horiz_label)
        subplot.set_ylabel(vert_label)
        #
        extent = None
        origin = None

        if h_limit is not None and v_limit is not None:
            extent = [- h_limit, h_limit, - v_limit, v_limit]
            origin = None #'lower'

        im = subplot.imshow(image, extent=extent, origin=origin, aspect='equal')
        return im

    def draw_scatter(self, figure, subplot, horizontal_values: np.ndarray, vertical_values: np.ndarray,
                     value: np.ndarray = None, title: str = None, horiz_label: str = None, vert_label: str = None,
                     h_limit: np.float64 = None, v_limit: np.float64 = None):

        subplot.cla()
        subplot.set_title(title)
        subplot.set_xlabel(horiz_label)
        subplot.set_ylabel(vert_label)
        #
        border_space = 1.1
        s = 10
        if value is not None:
            s = 0.01

        im = subplot.scatter(horizontal_values, vertical_values, s=s, c=value)
        subplot.set_aspect('equal', adjustable='box')

        if h_limit is not None and v_limit is not None:
            subplot.set_xlim([-h_limit, h_limit])
            subplot.set_ylim([-v_limit, v_limit])
        return im
