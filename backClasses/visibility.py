import numpy as np
from astropy.units import cds


class Visibility:

    def __init__(self, UVW: np.ndarray = None, frequency: np.ndarray = None, imagesize: int =None):
        self.deltay = None
        self.deltax = None
        self.deltav = None
        self.deltau = None
        self.uv_value_imag = None
        self.uv_value_real = None
        self.uv_value = None
        if isinstance(UVW, np.ndarray):
            self.UVW = UVW
        self.number_of_visibilities = len(UVW[0])
        self.max_uv_coordinate = self.max_uv_coordinate()
        self.weight = np.ones(self.number_of_visibilities)
        if frequency is None:
            raise ValueError('frequency has not been initialized')
        else:
            self.frequency = frequency
        self.calculate_deltas(imagesize)

    def max_uv_coordinate(self):  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        max_u = np.max(np.abs(self.UVW[0]))
        max_v = np.max(np.abs(self.UVW[1]))
        max_uv = max(max_u, max_v)
        return max_uv

    def set_value(self, uv_value):
        self.uv_value = uv_value
        self.uv_value_real = uv_value.real
        self.uv_value_imag = uv_value.imag

    def calculate_deltas(self, imagesize):
        epsilon = 1e-5
        maxuv = self.max_uv_coordinate + epsilon
        self.deltau = -2 * maxuv / (imagesize - 1)
        self.deltav = - self.deltau
        self.deltax = 1 / ((imagesize - 1) * self.deltau)
        self.deltay = 1 / ((imagesize - 1) * self.deltav)
