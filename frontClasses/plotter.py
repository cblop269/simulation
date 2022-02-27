import numpy as np
import matplotlib.colors as colors
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt


class Plotter:
    def __init__(self, nrows, ncols, width, height, w:float=None, h:float=None, l:float=None, t:float=None, b:float=None):
        '''self.figureA, self.axA = plt.subplots(1, 2, figsize=(9, 5))
        self.figureB, self.axB = plt.subplots(1, 3, figsize=(14, 5))
        self.figureD, self.axD = plt.subplots(2, 1, figsize=(9, 6))
        self.figureE, self.axE = plt.subplots(1, 1, figsize=(7, 6))
        self.figureF, self.axF = plt.subplots(1, 1, figsize=(7, 6))
        self.figureG, self.axG = plt.subplots(1, 1, figsize=(7, 6))
        # self.figure.subplots_adjust(hspace=1.5)
        self.figureA.subplots_adjust(wspace=.5, left=0.08, bottom=0.01)
        self.figureB.subplots_adjust(wspace=.5, left=0.08, bottom=0.01)
        self.figureD.subplots_adjust(hspace=.5, left=0.1, top=0.9)'''
        self.figure, self.ax = plt.subplots(nrows, ncols, figsize=(width, height))
        self.figure.subplots_adjust(wspace=w, hspace=h, left=l, top=t, bottom=b)

        # define colorbars
        self.cb_si = None
        self.cb_vi = None
        self.cb_gi = None
        self.cb_di = None

    def draw_noise(self, canvas, positionUV, noise_value, value, deltauv):
        self.axD[0].cla()
        self.axD[1].cla()
        self.axD[0].set_title('Visibilities Value')
        self.axD[0].set_xlabel('distance (λ)')
        self.axD[0].set_ylabel('value (Jy)')
        self.axD[1].set_title('Noise Values')
        self.axD[1].set_xlabel('distance (λ)')
        self.axD[1].set_ylabel('value (Jy)')
        # get the uv value
        noise_value = np.abs(noise_value)
        value = np.abs(value)
        # get the norm from u and v
        distance = np.linalg.norm(positionUV, axis=0)
        # draw plots
        self.axD[0].plot(distance[distance.argsort()], value[distance.argsort()], linewidth=0.2)
        self.axD[1].plot(distance[distance.argsort()], noise_value[distance.argsort()], linewidth=0.2)
        #self.axD[1].bar(distance, value)
        canvas.draw()

    def draw_plots_results(self, canvas, visibilities, grid_image, dirty_image, xy_delta, uv_delta):
        # get deltas
        deltax = xy_delta[0]
        deltay = xy_delta[1]
        deltau = uv_delta[0]
        deltav = uv_delta[1]

        # clear colorbars
        if self.cb_vi is not None and self.cb_gi is not None and self.cb_di is not None:
            self.cb_vi.remove()
            self.cb_gi.remove()
            self.cb_di.remove()

        # Visibilities
        horizontal_values = -visibilities[0] #* deltau
        vertical_values = visibilities[1] #* deltav
        uv_values = visibilities[2]
        limit = max(np.max(np.abs(visibilities[0])), np.max(np.abs(visibilities[1])))
        self.draw_scatter(self.axB[0], horizontal_values, vertical_values, uv_values,
                                'Visibilities', 'u (λ)', 'v (λ)', limit, limit, 0.01)
        im2 = self.draw_scatter(self.axE, horizontal_values, vertical_values, uv_values, '', 'u (λ)', 'v (λ)', limit,
                          limit, 2)

        # Grid Image
        m, n = grid_image.shape
        limit_hor = deltau * n / 2
        limit_vert = deltav * m / 2
        im3 = self.draw_image(self.axB[1], grid_image, 'Gridded Image', 'v (λ)', 'u (λ)', limit_hor, limit_vert)
        self.draw_image(self.axF, grid_image, '', 'v (λ)', 'u (λ)', limit_hor, limit_vert)

        # Dirty Image
        m, n = dirty_image.shape
        limit_hor = deltax * n / 2
        limit_vert = deltay * m / 2
        self.draw_image(self.axB[2], dirty_image, 'Dirty Image', 'y (v)', 'x (v)', limit_hor, limit_vert)
        im4 = self.draw_image(self.axG, dirty_image, '', 'y (v)', 'x (v)', limit_hor, limit_vert)

        # draw colorbars
        self.cb_vi = self.figureB.colorbar(im2, ax=self.axE, fraction=0.046, pad=0.04)
        #self.cb_gi = self.figureB.colorbar(im3, ax=self.axF, fraction=0.046, pad=0.04)
        self.cb_di = self.figureB.colorbar(im4, ax=self.axG, fraction=0.046, pad=0.04)
        #self.cb_vi = self.figureB.colorbar(im2, ax=self.axB[0], orientation='horizontal', pad=0.15)
        self.cb_gi = self.figureB.colorbar(im3, ax=self.axB[1], orientation='horizontal', pad=0.15)
        #self.cb_di = self.figureB.colorbar(im4, ax=self.axB[2], orientation='horizontal', pad=0.15)
        self.cb_vi.set_label('value applying log(|uv value| + 1)')
        self.cb_gi.set_label('value applying log(|uv value| + 1)')
        self.cb_di.set_label('value applying log(|uv value| + 1)')

        canvas.draw()

    def draw_image(self, subplot, image: np.ndarray, title: str = None, horiz_label: str = None,
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
        subplot.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        subplot.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        return im

    def draw_scatter(self, subplot, horizontal_values: np.ndarray, vertical_values: np.ndarray,
                     value: np.ndarray = None, title: str = None, horiz_label: str = None, vert_label: str = None,
                     h_limit: np.float64 = None, v_limit: np.float64 = None, dot_size: float = None):
        subplot.cla()

        subplot.set_title(title)
        subplot.set_xlabel(horiz_label)
        subplot.set_ylabel(vert_label)
        #
        border_space = 1.1

        im = subplot.scatter(horizontal_values, vertical_values, s=dot_size, c=value)
        subplot.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        subplot.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        subplot.set_aspect('equal', 'box')

        if h_limit is not None and v_limit is not None:
            subplot.set_xlim([-h_limit * border_space, h_limit * border_space])
            subplot.set_ylim([-v_limit * border_space, v_limit * border_space])

        else:
            h_center = (np.max(horizontal_values.value) + np.min(horizontal_values.value)) / 2
            v_center = (np.max(vertical_values.value) + np.min(vertical_values.value)) / 2
            h_distance = (np.max(horizontal_values.value) - np.min(horizontal_values.value)) / 2
            v_distance = (np.max(vertical_values.value) - np.min(vertical_values.value)) / 2
            distance = max(h_distance, v_distance)
            subplot.set_xlim([h_center - (distance * border_space), h_center + (distance * border_space)])
            subplot.set_ylim([v_center - (distance * border_space), v_center + (distance * border_space)])

        return im
