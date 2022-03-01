import numpy as np
import matplotlib.colors as colors
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt


class Plotter:
    def __init__(self, nrows, ncols, width, height, w:float=None, h:float=None, l:float=None, t:float=None, b:float=None):
        self.figure, self.ax = plt.subplots(nrows, ncols, figsize=(width, height))
        self.figure.subplots_adjust(wspace=w, hspace=h, left=l, top=t, bottom=b)

    def draw_image(self, subplot, image: np.ndarray, title: str = None, horiz_label: str = None,
                   vert_label: str = None, h_limit: np.float64 = None, v_limit: np.float64 = None):
        """
        Function that plot an image, the image is an equally numpy array of N x N
        :param subplot: The subplot where the graph will be
        :param image: Array with the values
        :param title: of the plot
        :param horiz_label: The horizontal label of the plot
        :param vert_label: The vertical label of the plot
        :param h_limit: Horizontal axe limit
        :param v_limit: Vertical axe limit
        :return The reference of the plot
        """
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
        """
        Function that apply scatter function in the visibilities
        :param subplot: The subplot where the graph will be
        :param horizontal_values: Horizontal position, generally u axe
        :param vertical_values: Vertical position, generally v axe
        :param value: The value of a determined uv position
        :param title: of the plot
        :param horiz_label: The horizontal label of the plot
        :param vert_label: The vertical label of the plot
        :param h_limit: Horizontal axe limit
        :param v_limit: Vertical axe limit
        :param dot_size: The size of the points
        :return The reference of the plot
        """
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
