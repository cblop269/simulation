import numpy as np
import matplotlib.colors as colors
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt


class Plotter:
    def __init__(self):
        self.figureA, self.axA = plt.subplots(1, 2, figsize=(9, 5))
        self.figureB, self.axB = plt.subplots(1, 4, figsize=(14, 5))
        # self.figureA = Figure(figsize=(9, 5), dpi=100)
        # self.figureB = Figure(figsize=(14, 5), dpi=100)
        # self.figure.subplots_adjust(hspace=1.5)
        self.figureA.subplots_adjust(wspace=.5, left=0.08, bottom=0.01)
        self.figureB.subplots_adjust(wspace=.5, left=0.08, bottom=0.01)

    def draw_inputs(self, canvas, image, antenna_pos, delta_X: np.float64 = None, delta_Y: np.float64 = None):

        # 1
        self.draw_image(self.figureA, self.axA[0], image, 'Sky Image', 'x (pixel)', 'y (pixel)')
        # 2
        self.draw_scatter(self.figureA, self.axA[1], antenna_pos[:, 0], antenna_pos[:, 1], None,
                          'Antenna configuration', 'x (m)', 'y (m)')
        '''self.axA[0].cla()
        self.axA[0].set_title('Sky Image')
        self.axA[0].set_xlabel('x(pixel)')
        self.axA[0].set_ylabel('y(pixel)')
        im_si = self.axA[0].imshow(image, vmin=0, vmax=np.max(image), aspect='equal', cmap='inferno')

        self.axA[1].cla()
        self.axA[1].set_title('Visibilities')
        self.axA[1].set_xlabel('x(pixel)')
        self.axA[1].set_ylabel('y(pixel)')
        value = None
        im_si = self.axA[1].scatter(antenna_pos[:, 0], antenna_pos[:, 1], vmin=0, vmax=np.max(value), cmap='inferno',
                                    c=value)
        # self.axA[1].set_xlim([x1 * 1.1, x2 * 1.1])
        # self.axA[1].set_ylim([y1 * 1.1, y2 * 1.1])
        self.axA[1].set_aspect('equal', adjustable='box')
        # self.ax[1].axis([x1, x2, y1, y2])'''
        #
        canvas.draw()

    def draw_plots_results(self, canvas, fft_image, visibilities, grid_image, dirty_image, xy_delta, uv_delta):
        #
        deltax = xy_delta[0]
        deltay = xy_delta[1]
        deltau = uv_delta[0]
        deltav = uv_delta[1]
        # 1
        m, n = fft_image.shape
        limit_hor = deltau * n / 2
        limit_vert = deltav * m / 2
        self.draw_image(self.figureB, self.axB[0], fft_image, 'Transform Image', 'u (λ)', 'v (λ)', limit_hor,
                        limit_vert)
        # 2
        limit_hor = deltau * visibilities.UVW[0]
        limit_vert = deltav * visibilities.UVW[1]
        self.draw_scatter(self.figureB, self.axB[1], limit_hor, limit_vert, np.log(abs(visibilities.value) + 1),
                          'Visibilities', 'u (λ)', 'v (λ)', visibilities.max_uv_coordinate * deltau,
                          visibilities.max_uv_coordinate * deltav)
        # 3
        m, n = grid_image.shape
        limit_hor = deltau * n / 2
        limit_vert = deltav * m / 2
        self.draw_image(self.figureB, self.axB[2], grid_image, 'Gridded Image', 'v (λ)', 'u (λ)', limit_hor, limit_vert)
        # 4
        m, n = dirty_image.shape
        limit_hor = deltax * n / 2
        limit_vert = deltay * m / 2
        self.draw_image(self.figureB, self.axB[3], dirty_image, 'Dirty image', 'y (v)', 'x (v)', limit_hor, limit_vert)

        canvas.draw()

    def draw_image(self, figure, subplot, image: np.ndarray, title: str = None, horiz_label: str = None,
                   vert_label: str = None, h_limit: np.float64 = None, v_limit: np.float64 = None):
        subplot.cla()
        subplot.set_title(title)
        subplot.set_xlabel(horiz_label)
        subplot.set_ylabel(vert_label)
        #
        extent = None
        if h_limit is not None and v_limit is not None:
            extent = [- h_limit, h_limit, - v_limit, v_limit]

        im = subplot.imshow(image, cmap='inferno', extent=extent, origin='lower', aspect='equal')
        figure.colorbar(im, ax=subplot, orientation='horizontal', pad=0.1)

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

        im = subplot.scatter(horizontal_values, vertical_values, s=s, c=value, cmap='inferno')
        subplot.set_aspect('equal', adjustable='box')
        figure.colorbar(im, ax=subplot, orientation='horizontal', pad=0.1)
        if h_limit is not None and v_limit is not None:
            print('1 ', h_limit)
            print('2 ', v_limit)
            subplot.set_xlim([-h_limit, h_limit])
            subplot.set_ylim([-v_limit, v_limit])
