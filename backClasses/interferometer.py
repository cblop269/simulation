import math
import sys
import numpy as np
import astropy.units as u
from time import time
from numba import jit
from scipy.constants import c
from math import floor
from astropy.io import fits
from scipy import interpolate
from pynufft import NUFFT
import finufft

from backClasses.visibility import Visibility
from backClasses.baseline import Baseline


class Interferometer:

    def __init__(self, antenna_route: str = None):  # , telescope_lat, source_decl, ha_start, ha_end, dt, obs_freq):
        self.sky_image = None
        self.fft_image = None
        self.delta_X = None
        self.delta_U = None
        self.visibilities = None
        if antenna_route is not None:
            # assignment of antenna position
            self.antenna_pos = self.read_antenna_config(antenna_route)
            # assignment of baseline
            baseline = self.compute_baselines()
            self.baseline = Baseline(baseline, self.antenna_pos)
        else:
            raise ValueError('antenna route has not been initialized')

    def read_antenna_config(self, route: str):
        try:
            file_observatory = np.loadtxt(route, skiprows=3, usecols=(0, 1, 2))
            return file_observatory
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            sys.exit()

    def compute_baselines(self):
        # antpos[] are XYZ coordinates in meters. XY aligned with celestial equator
        # Z pointing to NCP
        # antpos[1,:] -> East,  antpos[2, :] -> North, antpos[3,:] -> North Celestial Pole
        # B : baselines vectors [:, 1] = x, [:, 2] = y, [:, 3] = z
        # antpos[] antenna positions in meters [:, 1] = x, [:, 2] = y, [:, 3] = z

        # map_baseline = {}
        ant_set1 = np.repeat(self.antenna_pos, repeats=len(self.antenna_pos), axis=0)
        ant_set2 = np.concatenate([self.antenna_pos] * len(self.antenna_pos), axis=0)
        B2 = ant_set1 - ant_set2
        # take off the [0,0,0] baseline
        B2 = B2[B2[:, 0].argsort(kind='mergesort')]
        valid_baselines = int(len(B2) / 2) - int(len(self.antenna_pos) / 2)
        B2 = np.concatenate((B2[:valid_baselines, :], B2[valid_baselines + len(self.antenna_pos):, :]), axis=0)
        B2 = B2.transpose()
        return B2

    def compute_hour_angle(self, ha_start: float = None, ha_end: float = None, dt: int = None):
        # ha_rad[] vector with hour angle samplings in radians!!!
        # ha_star, ha_end in Hours (for instance -1 hora, 1 hora)
        # dt sampling interval in seconds (for instance 60 seconds)
        # check hour angle limits
        if ha_start < -12.0 or ha_start is None:
            ha_start = -12.0
        if ha_end > 12.0 or ha_end is None:
            ha_end = 12.0
        if dt is None:
            dt = 100
        # convert from degrees to radian
        ha_start_rad = 15 * ha_start * np.pi / 180.0
        ha_end_rad = 15 * ha_end * np.pi / 180.0
        dt_rad = 15 * dt * np.pi / (3600 * 180)
        ha_rad = np.arange(ha_start_rad, ha_end_rad, dt_rad)
        return ha_rad

    def create_uv_position(self, telescope_lat: float = None, source_decl: float = None,
                           ha: np.ndarray = None, obs_freq: float = None):
        # obs_freq is observing frequency in Hz (continuum ??)
        # telescope_lat in degrees. Example -23
        # source_decl Source declination in degrees. For instance 18
        decl_rad = 0
        lat = 0
        if telescope_lat is not None and source_decl is not None:
            decl_rad = source_decl * np.pi / 180
            lat = (-90 + telescope_lat) * np.pi / 180

        # Rotate around x to rise w -Dec degrees
        R2 = np.array([[1, 0, 0],
                        [0, np.cos((np.pi / 2) - decl_rad), -np.sin((np.pi / 2) - decl_rad)],
                        [0, np.sin((np.pi / 2) - decl_rad), np.cos((np.pi / 2) - decl_rad)]])

        # Rotate around x to correct for telescope latitude
        R4 = np.array([[1, 0, 0],
                        [0, np.cos(lat), -np.sin(lat)],
                        [0, np.sin(lat), np.cos(lat)]])

        lambda_num = c / obs_freq
        self.compute_uvw(lambda_num, R2, R4, ha, obs_freq)

    #@jit(forceobj=True)
    def compute_uvw(self, lambda_num: float = None, r2: np.ndarray = None, r4: np.ndarray = None,
                    ha: np.ndarray = None, obs_freq: float = None):

        # create a R3 for every HA
        R3 = np.array([np.cos(ha), -np.sin(ha), np.zeros([len(ha)]),
                       np.sin(ha), np.cos(ha), np.zeros([len(ha)]),
                       np.zeros([len(ha)]), np.zeros([len(ha)]), np.ones([len(ha)])])


        # re-order the array
        print('r3 ', R3.shape)
        R3 = np.swapaxes(R3, 1, 0)
        print('r3 ', R3.shape)
        # operate every element in the array
        uvw = np.apply_along_axis(self.apply_rotation_matrix, 1, R3, r4, r2, lambda_num)
        #uvw = r4[:, :, np.newaxis] @ R3 @ r2[:, :, np.newaxis]
        #uvw = uvw.transpose()
        #uvw = (uvw * self.baseline.baselines) / lambda_num

        # re-order the visibilities to be easier to scattering
        a, b, c = np.shape(uvw)
        '''print('original ', uvw.shape)
        print('uvw.v1: ', uvw)
        print('transpose ', uvw.transpose([1, 0, 2]).shape)
        print('uvw.v2: ', uvw.transpose([1, 0, 2]))
        print('reshape ', uvw.transpose([1, 0, 2]).reshape(b, a * c).shape)
        print('uvw.v3: ', uvw.transpose([1, 0, 2]).reshape(b, a * c))'''
        print('1111111 ', uvw)
        print('2222222 ', uvw.transpose([1, 0, 2]))
        uvw = uvw.transpose([1, 0, 2]).reshape(b, a * c)
        #print('reshape ', uvw.shape)
        self.visibilities = Visibility(uvw, obs_freq)

    def apply_rotation_matrix(self, matrix3: np.ndarray, matrix4: np.matrix, matrix2: np.matrix,
                              lambda_num: float):
        uvw = matrix4 @ matrix3.reshape(3, 3) @ matrix2
        uvw = uvw.transpose()
        uvw = (uvw @ self.baseline.baselines) / lambda_num
        return np.asarray(uvw)

    def read_image(self, route: str = None):
        try:
            with fits.open(route) as image:
                image_data = image[0].data
            self.sky_image = image_data
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            sys.exit()

    def write_image(self, data, route: str = None, ):
        fits.writeto('new1.fits', data)


    def transform_fft(self):
        # transform with fast fourier transform in 2 dimensions
        fft_image = np.fft.fft2(self.sky_image)
        self.fft_image = fft_image

    def transform_nufft(self):
        m, n = np.shape(self.sky_image)
        uvw = self.visibilities.UVW.transpose()
        uvw = uvw[:, :2]
        uvw = uvw[:, ::-1]
        uvw = uvw * np.pi / (np.max(np.abs(uvw)))  # convert to radians
        #
        NufftObj = NUFFT()
        Nd = (m, n)  # image size
        Kd = (1024, 1024)  # k-space size
        Jd = (10, 10)
        NufftObj.plan(uvw, Nd, Kd, Jd)

        #values = NufftObj.forward(self.sky_image)
        ###############################
        u = self.visibilities.UVW[0] * np.pi / (np.max(np.abs(self.visibilities.UVW)))
        v = self.visibilities.UVW[1] * np.pi / (np.max(np.abs(self.visibilities.UVW)))

        complex_sky_image = self.sky_image + 1j * np.zeros_like(self.sky_image)

        # the 2D transform outputs f array of shape (N1, N2)
        f = finufft.nufft2d2(v, u, complex_sky_image, eps=1e-12)

        self.visibilities.set_value(f)

    def inverse_transform_fft(self):
        dirty_image = np.fft.ifftshift(self.gridded_vo)
        dirty_image = np.fft.ifft2(dirty_image)
        dirty_image = np.fft.fftshift(dirty_image)
        self.dirty_image = dirty_image


    def bilinear_interpolation(self):
        # other variables
        N = max(len(self.fft_image[0]), len(self.fft_image))
        delta_X = 1 / (2 * self.visibilities.max_uv_coordinate)
        delta_X = delta_X / 7
        delta_U = 1 / (N * delta_X)

        # find the list of j and i equivalent to four points Q11 points
        list_j = np.floor(self.visibilities.UVW[0] / delta_U) + (N / 2)  # horizontal
        list_i = np.floor(self.visibilities.UVW[1] / delta_U) + (N / 2)  # vertical

        # then is applied  the linear interpolation
        l_alpha = (self.visibilities.UVW[0] - ((list_j - (N / 2)) * delta_U)) / delta_U
        l_beta = (self.visibilities.UVW[1] - ((list_i - (N / 2)) * delta_U)) / delta_U
        uv_value = np.apply_along_axis(self.calculate_uv_value, 0, [list_j, list_i, l_alpha, l_beta])

        # values are created for the visibilities object
        self.visibilities.set_value(uv_value)

    def calculate_uv_value(self, data: np.ndarray = None):
        # current u position = data[0]
        # current v position = data[1]
        # current alpha = data[2]
        # current beta = data[3])
        j = int(data[0])
        i = int(data[1])

        # linear interpolation in the horizontal-direction
        fxy1 = (data[2] * self.fft_image[i][j + 1]) + ((1 - data[2]) * self.fft_image[i][j])
        fxy2 = (data[2] * self.fft_image[i + 1][j + 1]) + ((1 - data[2]) * self.fft_image[i + 1][j])

        # linear interpolation in the vertical-direction
        uv_value = (data[3] * fxy2) + ((1 - data[3]) * fxy1)
        return uv_value

    @jit(forceobj=True)
    def fourier_series(self):
        m, n = np.shape(self.sky_image)
        uvw = self.visibilities.UVW.transpose()
        uvw = uvw[:, :2]

        # convert to radians
        uvw = uvw * np.pi / (np.max(np.abs(uvw)))
        #
        sky_image = self.sky_image.reshape(m * n)
        x = np.linspace(- n / 2, (n / 2) - 1, n)
        y = np.linspace(- m / 2, (m / 2) - 1, m)
        xg, yg = np.meshgrid(x, y)
        xg = xg.reshape(m * n)
        yg = yg.reshape(m * n)
        #
        result = np.apply_along_axis(self.test_function_2, 1, uvw, sky_image, yg, xg)
        # result = np.apply_along_axis(self.test_function_3, 1, uvw, self.sky_image, y, x, n, m)
        self.visibilities.set_value(result)

    @jit(forceobj=True)
    def test_function_2(self, uvw, sky_image, j_set, i_set):
        value = (uvw[0] * i_set) + (uvw[1] * j_set)
        e = np.cos(value) + 1j * np.sin(value)
        new_visibility = sky_image * (-e)
        new_visibility = np.sum(new_visibility)
        return new_visibility

    def test_function_3(self, uvw, sky_image, j_set, i_set, n, m):
        new_visibility = np.empty([0])
        for j in j_set:
            for i in i_set:
                value = (uvw[0] * i) + (uvw[1] * j)
                e = np.cos(value) + 1j * np.sin(value)
                result = sky_image[int(j + (m / 2))][int(i + (n / 2))] * (-e)
                new_visibility = np.append(new_visibility, result)
        return np.sum(new_visibility)

    def gridder(self, imagesize, scheme, robust):
        # vis numpy array with -> [u v w Re Im We f]
        # scheme = 1, robust = 0
        # return gridded image, deltau, max uv, delta x
        uvw = self.visibilities.UVW.transpose()
        uvw = uvw[:, :2]
        #uvw = uvw * np.pi / (np.max(np.abs(uvw)))  # convert to radians

        maxuv = np.amax(np.abs(uvw))
        epsilon = 1e-5
        maxuv = maxuv + epsilon

        z = len(self.visibilities.value)
        pos_u_index = np.zeros(z).astype(int)
        pos_v_index = np.zeros(z).astype(int)
        self.deltau = -2 * maxuv / (imagesize - 1)
        self.deltav = - self.deltau
        self.deltax = 1 / (imagesize * self.deltau)
        self.deltay = 1 / (imagesize * self.deltav)
        print('deltau ', self.deltau)
        print('deltax ', self.deltax)

        # Gridding weights if ROBUST OR UNIFORM
        gridded_weights = np.zeros((imagesize, imagesize))
        for k in range(0, z):
            j = int(np.round(uvw[k][0] / self.deltau) + imagesize / 2)
            i = int(np.round(uvw[k][1] / self.deltav) + imagesize / 2)
            pos_v_index[k] = i
            pos_u_index[k] = j
            gridded_weights[i][j] += self.visibilities.weight[k]
        # Selecting a scheme for weights
        self.scheme_weights(gridded_weights, z, pos_u_index, pos_v_index, scheme, robust)

        #
        # Gridding visibilities and weights with a scheme
        gridded_weights = np.zeros((imagesize, imagesize))
        gridded_Vo = np.zeros((imagesize, imagesize)) + 1.0j * np.zeros((imagesize, imagesize))

        for k in range(0, z):
            j = int(np.round(uvw[k][0] / self.deltau) + imagesize / 2)
            i = int(np.round(uvw[k][1] / self.deltav) + imagesize / 2)
            gridded_weights[i][j] += self.visibilities.weight[k]
            gridded_Vo[i][j] += self.visibilities.weight[k] * self.visibilities.value[k]

        # Check rows and columns where gridded weights are greater than 1
        rows, columns = np.where(gridded_weights > 0)
        # Divide the gridded visibilities by the gridded weights
        gridded_Vo[rows, columns] /= gridded_weights[rows, columns]

        self.gridded_image = gridded_weights
        self.gridded_vo = gridded_Vo

    def scheme_weights(self, gridded_weights, z, pos_u_index, pos_v_index, scheme, robust):
        # Selecting a scheme for weights
        average_weights = np.sum(gridded_weights ** 2) / np.sum(self.visibilities.weight)
        for k in range(0, z):
            if (scheme == 0):
                self.visibilities.weight[k] = self.visibilities.weight[k]
            elif (scheme == 1):
                self.visibilities.weight[k] /= gridded_weights[pos_v_index[k]][pos_u_index[k]]
            elif (scheme == 2):
                self.visibilities.weight[k] *= np.sqrt(self.visibilities.UVW[0][k] ** 2 + self.visibilities.UVW[1][k] ** 2)
            else:
                f2 = (5.0 * np.power(10.0, -robust)) ** 2 / average_weights
                self.visibilities.weight[k] /= (1.0 + gridded_weights[pos_v_index[k]][pos_u_index[k]] * f2)


    def run(self, telescope_lat: float = None, source_decl: float = None, ha_start: float = None,
            ha_end: float = None, dt: int = None, obs_freq: float = None, usefft: bool = None):
        # toma: la imagen de entrada(sky_image), el objeto visibilities(self)

        # uv positions
        ha = self.compute_hour_angle(ha_start, ha_end, dt)
        self.create_uv_position(telescope_lat, source_decl, ha, obs_freq)
        # image with fft
        self.transform_fft()
        self.fft_image = np.fft.fftshift(self.fft_image)

        # values of visibilities
        if usefft:
            self.bilinear_interpolation()

        else:
            self.transform_nufft()
            #self.fourier_series()
        self.gridder(len(self.fft_image), 1, 0)
        self.inverse_transform_fft()
