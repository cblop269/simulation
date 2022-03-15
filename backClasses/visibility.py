import numpy as np
from astropy.units import cds


class Visibility:
    deltay = 1
    deltax = -1
    deltav = 1
    deltau = -1
    def __init__(self, UVW: np.ndarray = None, frequency: np.ndarray = None, imagesize: int =None):
        self.uv_value_imag = None
        self.uv_value_real = None
        self.uv_value = None
        if isinstance(UVW, np.ndarray):
            self.UVW = UVW
        self.number_of_visibilities = len(UVW[0])
        self.max_uv_coordinate = self.__max_uv_coordinate()
        self.weight = np.ones(self.number_of_visibilities) * (cds.Jy / cds.Jy)
        if frequency is None:
            raise ValueError('frequency has not been initialized')
        else:
            self.frequency = frequency
        self.__calculate_deltas(imagesize)

    def __max_uv_coordinate(self):
        """
        Function that calculates the max uv position
        :return max uv position between axe u and axe v
        """
        max_u = np.max(np.abs(self.UVW[0]))
        max_v = np.max(np.abs(self.UVW[1]))
        max_uv = max(max_u, max_v)
        return max_uv

    def set_value(self, uv_value):
        """
        Function that set the uv values
        :param uv_value: the value of an uvw position
        """
        self.uv_value = uv_value
        self.uv_value_real = uv_value.real
        self.uv_value_imag = uv_value.imag

    def __calculate_deltas(self, imagesize):
        """
        Function that calculates and set the delta values
        :param imagesize, the number of pixels N in a image (of a N x N image)
        """
        epsilon = 1e-5
        maxuv = self.max_uv_coordinate + epsilon
        self.deltau = -2 * maxuv / (imagesize - 1)
        self.deltav = - self.deltau
        self.deltax = 1 / ((imagesize - 1) * self.deltau)
        self.deltay = 1 / ((imagesize - 1) * self.deltav)
