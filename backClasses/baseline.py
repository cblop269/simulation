import numpy as np


class Baseline:

    def __init__(self, B, antpos):
        self.baselines = B
        self.antenna = antpos
        #self.map = map_baseline
        self.max_baseline = self.calculate_max()

    def search_antenna(self, baseline):
        print(self.antenna)
        pair = self.map[baseline[0], baseline[1], baseline[2]]
        antenna1 = self.antenna[pair[0]]
        antenna2 = self.antenna[pair[1]]
        return [antenna1, antenna2]

    def calculate_max(self):
        baseline_list = (np.swapaxes(self.baselines, 1, 0))
        magnitude_list = np.empty([0])
        for i in baseline_list:
            magnitude_list = np.append(magnitude_list, np.linalg.norm(i))
        return np.max(magnitude_list)
