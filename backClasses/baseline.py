import numpy as np

class Baseline:

    def __init__(self, B, antenna_pos, antenna_radius):
        self.baselines = B
        self.antenna = antenna_pos
        self.max_baseline = self.__calculate_max_baseline()
        self.antenna_number = len(antenna_pos)
        self.antenna_area = self.__calculate_antenna_area(antenna_radius)

    def __calculate_max_baseline(self):
        """
        Calculate the max baseline
        :return the max baseline length
        """
        baseline_list = np.apply_along_axis(np.linalg.norm, 0, self.baselines)
        return np.max(baseline_list)

    def __calculate_antenna_area(self, antenna_radius):
        """
        Calculate the max baseline
        :return the max baseline length
        """
        #   save parameters
        antenna_area = np.pi * (antenna_radius ** 2)
        return antenna_area