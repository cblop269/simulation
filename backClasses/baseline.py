import numpy as np

class Baseline:

    def __init__(self, B, antpos):
        self.baselines = B
        self.antenna = antpos
        self.max_baseline = self.calculate_max()

    def search_antenna(self, baseline):
        """
        :param
        :return
        """
        pair = self.map[baseline[0], baseline[1], baseline[2]]
        antenna1 = self.antenna[pair[0]]
        antenna2 = self.antenna[pair[1]]
        return [antenna1, antenna2]

    def calculate_max(self):
        """
        Calculate the max baseline
        :return the max baseline length
        """
        baseline_list = np.apply_along_axis(np.linalg.norm, 0, self.baselines)
        return np.max(baseline_list)
