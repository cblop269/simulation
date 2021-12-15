import numpy as np


class Visibility:

    def __init__(self, UVW: np.ndarray = None, frequency: np.ndarray = None):
        if isinstance(UVW, np.ndarray):
            self.UVW = UVW
        self.number_of_visibilities = len(UVW[0])
        self.max_uv_coordinate = self.max_uv_coordinate()
        self.weight = np.ones(self.number_of_visibilities)
        if frequency is None:
            raise ValueError('frequency has not been initialized')
        else:
            self.frequency = frequency
        #self.value = 0
        #self.weight = 0
        #self.baselines = 0
        #self.noise = 0

    def max_uv_coordinate(self):  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        max_u = np.max(np.abs(self.UVW[0]))
        max_v = np.max(np.abs(self.UVW[1]))
        max_uv = max(max_u, max_v)
        return max_uv

    def set_value(self, uv_value):
        self.value = uv_value
        #print(self.value)