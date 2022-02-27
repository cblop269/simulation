import numpy as np
from astropy.units import cds

class Gridder:

    def __init__(self, visibilities, noise, imagesize:int=None, scheme:int=None, robust:int=None):
        self.gridded_vo = None
        self.gridded_image = None
        self.gridded_weights = np.zeros((imagesize, imagesize))
        # vis numpy array with -> [u v w Re Im We f]
        # scheme = 1, robust = 0
        # return gridded image
        # parameters
        u_pos = visibilities.UVW[0]
        v_pos = visibilities.UVW[1]
        uv_value = visibilities.uv_value + noise
        deltau = visibilities.deltau
        deltav = visibilities.deltav
        weight = visibilities.weight
        z = len(uv_value)
        self.pos_u_index = np.zeros(z).astype(int)
        self.pos_v_index = np.zeros(z).astype(int)

        # Gridding weights if ROBUST OR UNIFORM
        if scheme != 0:
            self.grid_robust(z, u_pos, v_pos, deltau, deltav, weight, imagesize)
            # Selecting a scheme for weights
            weight = self.scheme_weights(z, u_pos, v_pos, scheme, robust, weight)
        #
        self.grid_visibilities(z, u_pos, v_pos, deltau, deltav, weight, imagesize, uv_value)

    def grid_robust(self, z, u, v, deltau, deltav, weight, imagesize):
        for k in range(0, z):
            j = int(np.round(u[k] / deltau) + imagesize / 2)
            i = int(np.round(v[k] / deltav) + imagesize / 2)
            self.pos_v_index[k] = i
            self.pos_u_index[k] = j
            self.gridded_weights[i][j] += weight[k]

    def scheme_weights(self, z, u, v, scheme, robust, weight):

        average_weights = np.sum(self.gridded_weights ** 2) / np.sum(weight)
        for k in range(0, z):
            if scheme == 1:
                weight[k] /= self.gridded_weights[self.pos_v_index[k]][self.pos_u_index[k]]
            elif scheme == 2:
                weight[k] *= np.sqrt(u[k] ** 2 + v[k] ** 2)
            else:
                f2 = (5.0 * np.power(10.0, -robust)) ** 2 / average_weights
                weight[k] /= (1.0 + self.gridded_weights[self.pos_v_index[k]][self.pos_u_index[k]] * f2)
        return weight


    def grid_visibilities(self, z, u, v, deltau, deltav, weight, imagesize, uv_value):
        # Gridding visibilities and weights with a scheme
        gridded_weights = np.zeros((imagesize, imagesize)) * (cds.Jy / cds.Jy)
        gridded_Vo = np.zeros((imagesize, imagesize)) + 1.0j * np.zeros((imagesize, imagesize)) * cds.Jy

        for k in range(0, z):
            j = int(np.round(u[k] / deltau) + imagesize / 2)
            i = int(np.round(v[k] / deltav) + imagesize / 2)
            gridded_weights[i][j] += weight[k]
            gridded_Vo[i][j] += weight[k] * uv_value[k]

        # Check rows and columns where gridded weights are greater than 1
        rows, columns = np.where(gridded_weights > 0)
        # Divide the gridded visibilities by the gridded weights
        gridded_Vo[rows, columns] /= gridded_weights[rows, columns]

        self.gridded_image = gridded_weights
        self.gridded_vo = gridded_Vo
