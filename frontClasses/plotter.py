import numpy as np
import matplotlib.colors as colors
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable


class Plotter:
    def __init__(self):
        self.figureA = Figure(figsize=(9, 5), dpi=100)
        self.figureB = Figure(figsize=(13, 5), dpi=100)
        # self.figure.subplots_adjust(hspace=1.5)
        self.figureA.subplots_adjust(wspace=.5, left=0.08, bottom=0.01)
        self.figureB.subplots_adjust(wspace=.5, left=0.08, bottom=0.01)

    def draw_inputs(self, canvas, image, antenna_pos, delta_X: np.float64 = None, delta_Y: np.float64 = None):

        self.figureA.clear()
        sky_axe = self.figureA.add_subplot(1, 2, 1)
        antenna_axe = self.figureA.add_subplot(1, 2, 2)
        # a
        im = self.draw_image(sky_axe, image)
        sky_axe.set_title('Sky image', fontsize=10, pad=20)
        sky_axe.set_ylabel('y (pixel)', fontsize=10)
        sky_axe.set_xlabel('x (pixel)', fontsize=10)
        divider = make_axes_locatable(sky_axe)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        self.figureA.colorbar(im, ax=sky_axe, cax=cax)
        # b
        im2 = self.draw_scatter(antenna_axe, antenna_pos[:, 0], antenna_pos[:, 1])
        antenna_axe.set_title('Antenna config', fontsize=10, pad=20)
        antenna_axe.set_ylabel('y (m)', fontsize=10)
        antenna_axe.set_xlabel('x (m)', fontsize=10)
        #
        canvas.draw()

    def draw_plots_results(self, canvas, fft_image, visibilities, grid_image, dirty_image, xy_limits, uv_limits):

        self.figureB.clear()
        transform_axe = self.figureB.add_subplot(1, 4, 1)
        visibilities_axe = self.figureB.add_subplot(1, 4,  2)
        gridded_axe = self.figureB.add_subplot(1, 4,  3)
        dirty_axe = self.figureB.add_subplot(1, 4,  4)
        # 1
        m, n = np.shape(fft_image)
        u_limit = (n / 2) * uv_limits[0]
        v_limit = -(m / 2) * uv_limits[1]
        im = self.draw_image(transform_axe, fft_image, u_limit, v_limit)
        transform_axe.set_title('Transform image', fontsize=10, pad=20)
        transform_axe.set_ylabel('v (λ)', fontsize=10)
        transform_axe.set_xlabel('u (λ)', fontsize=10)
        divider = make_axes_locatable(transform_axe)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = self.figureB.colorbar(im, ax=transform_axe, cax=cax)#, orientation='horizontal')
        # 2
        im = self.draw_scatter(visibilities_axe, visibilities.UVW[0] / 1000, visibilities.UVW[1] / 1000,
                               np.log(abs(visibilities.value) + 1), visibilities.max_uv_coordinate / 1000)
        visibilities_axe.set_title('Visibilities', fontsize=10, pad=20)
        visibilities_axe.set_ylabel('v (λ)', fontsize=10)
        visibilities_axe.set_xlabel('u (λ)', fontsize=10)
        divider = make_axes_locatable(visibilities_axe)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = self.figureB.colorbar(im, ax=visibilities_axe, cax=cax)#, orientation='horizontal')
        canvas.draw()
        # 3
        m, n = np.shape(fft_image)
        u_limit = (n / 2) * uv_limits[0]
        v_limit = (m / 2) * uv_limits[1]
        im = self.draw_image(gridded_axe, grid_image, u_limit, v_limit)
        dirty_axe.invert_xaxis()
        gridded_axe.set_title('Gridded image', fontsize=10, pad=20)
        gridded_axe.set_ylabel('v (λ)', fontsize=10)
        gridded_axe.set_xlabel('u (λ)', fontsize=10)
        divider = make_axes_locatable(gridded_axe)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = self.figureB.colorbar(im, ax=gridded_axe, cax=cax)#, orientation='horizontal')
        canvas.draw()
        # 4
        m, n = np.shape(dirty_image)
        x_limit = (n / 2) * xy_limits[0]
        y_limit = (m / 2) * xy_limits[1]
        im = self.draw_image(dirty_axe, dirty_image)#, x_limit, y_limit)
        dirty_axe.set_title('Dirty image', fontsize=10, pad=20)
        dirty_axe.set_ylabel('y (v)', fontsize=10)
        dirty_axe.set_xlabel('x (v)', fontsize=10)
        divider = make_axes_locatable(dirty_axe)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = self.figureB.colorbar(im, ax=dirty_axe, cax=cax)#, orientation='horizontal')
        canvas.draw()

    def draw_image(self, plot, image: np.ndarray, x_limit: np.float64 = None,
                   y_limit: np.float64 = None, origin: str = None):
        #
        '''if 1 > 2 is None:
            origin = 'lower'
            color = colors.LogNorm(vmin=np.min(image), vmax=np.max(image))
            x_limit = x_limit / 1000000
            y_limit = y_limit / 1000000'''
        extent = None
        if x_limit is not None and y_limit is not None:
            extent = [- x_limit, x_limit, - y_limit, y_limit]

        min = np.min(image[np.nonzero(image)])
        print('vmin ', np.min(image), '- vmax ', np.max(image))
        if np.min(image) < 1:
            min = 1

        im = plot.imshow(image, cmap='inferno', extent=extent, origin=origin,
                         norm=None)#colors.LogNorm(vmin=min, vmax=np.max(image)))
        plot.set_aspect('equal')
        return im

    def draw_scatter(self, plot, horizontal_values: np.ndarray, vertical_values: np.ndarray, value: np.ndarray = None,
                   limit_axe: np.float64 = None):
        #
        border_space = 1.1
        '''if 1 > 2 is None:
            factor = 1000
            limit_axe = d.max_uv_coordinate * border_space / factor'''
        s = 10
        color = None
        if value is not None:
            s = 0.01
            # color = colors.LogNorm(vmin=np.min(value), vmax=np.max(value))

        im = plot.scatter(horizontal_values, vertical_values, s=s, c=value, cmap='inferno',
                          norm=color)
        # limit = max(np.max(abs(horizontal_values)), np.max(abs(vertical_values)))
        #plot.set_xlim([np.min(horizontal_values) * 1.1, np.max(horizontal_values) * 1.1])
        # plot.set_ylim([np.min(vertical_values) * 1.1, np.max(vertical_values) * 1.1])
        if limit_axe is not None:
            plot.set_xlim([-limit_axe, limit_axe])
            plot.set_ylim([-limit_axe, limit_axe])
        plot.set_aspect('equal', adjustable='box')
        return im

    '''
    def draw_image(self, position: int, image: np.ndarray, canvas, x_limit: np.float64,
                   y_limit: np.float64, title: str = None):
        #
        plot = self.figureB.add_subplot(2, 2, position)

        origin = None
        color = None
        horizontal_label = ""
        vertical_label = ""
        if title is None:
            title = "Transform Image"
            origin = 'lower'
            color = colors.LogNorm(vmin=np.min(image), vmax=np.max(image))
            print('min ti: ', np.min(image))
            print('max ti: ', np.max(image))
            horizontal_label = "u (Mλ)"
            vertical_label = "u (Mλ)"
            x_limit = x_limit / 1000000
            y_limit = y_limit / 1000000

        im = plot.imshow(image, cmap='inferno',
                         extent=[- x_limit, x_limit, - y_limit, y_limit],
                         origin=origin, norm=color)
        #
        plot.set_title(title, fontsize=20, pad=20)
        plot.set_ylabel(vertical_label, fontsize=10)
        plot.set_xlabel(horizontal_label, fontsize=10)
        plot.set_aspect('equal')
        divider = make_axes_locatable(plot)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        self.figureB.colorbar(im, ax=plot, cax=cax)
        canvas.draw()

    def draw_sky_image(self, position: int, image: np.ndarray, canvas, x_limit: np.float64,
                       y_limit: np.float64, title: str = None):
        #
        plot = self.figureA.add_subplot(1, 2, position)
        plot.clear()
        origin = None
        color = None
        horizontal_label = "x (ν)"
        vertical_label = "y (ν)"
        if title is None:
            title = "Transform Image"
            origin = 'lower'
            color = colors.LogNorm(vmin=np.min(image), vmax=np.max(image))
            print('min ti: ', np.min(image))
            print('max ti: ', np.max(image))
            horizontal_label = "u (Mλ)"
            vertical_label = "u (Mλ)"
            x_limit = x_limit / 1000000
            y_limit = y_limit / 1000000

        im = plot.imshow(image, cmap='inferno',
                         extent=[- x_limit, x_limit, - y_limit, y_limit],
                         origin=origin, norm=color)
        # configuration
        plot.set_title(title, fontsize=20, pad=20)
        plot.set_ylabel(vertical_label, fontsize=10)
        plot.set_xlabel(horizontal_label, fontsize=10)
        plot.set_aspect('equal')
        # create colorbar
        divider = make_axes_locatable(plot)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        self.figureA.colorbar(im, ax=plot, cax=cax)
        # show image
        canvas.draw()

    # @jit(forceobj=True)
    def draw_uv_positions(self,  d, canvas, figure, title: str = None):
        #
        border_space = 1.1
        plot = figure.add_subplot(1, 2, 2)

        if title is None:
            value = np.ones(len(d))
            title = "antenna config"
            horizontal_label = "x (m)"
            vertical_label = "y (m)"
            limit_axe = d.max() * border_space
            horizontal_value = d[:, 0]
            vertical_value = d[:, 1]
            im = plot.scatter(horizontal_value, vertical_value, s=1, cmap='inferno')

        else:
            value = d.value
            value = np.log(abs(value) + 1)
            factor = 1000
            horizontal_label = "u (Kλ)"
            vertical_label = "v (Kλ)"
            horizontal_value = -d.UVW[0] / factor
            vertical_value = d.UVW[1] / factor
            limit_axe = d.max_uv_coordinate * border_space / factor
            im = plot.scatter(horizontal_value, vertical_value, s=0.1, c=value, cmap='inferno',
                          norm=colors.LogNorm(vmin=np.min(value), vmax=np.max(value)))
            plot.set_xlim([-limit_axe, limit_axe])
            plot.set_ylim([-limit_axe, limit_axe])
            plot.invert_xaxis()
            divider = make_axes_locatable(plot)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            figure.colorbar(im, ax=plot, cax=cax)

        plot.set_title(title, fontsize=20, pad=20)
        plot.set_ylabel(vertical_label, fontsize=10)
        plot.set_xlabel(horizontal_label, fontsize=10)
        plot.set_aspect('equal')

        canvas.draw()'''
