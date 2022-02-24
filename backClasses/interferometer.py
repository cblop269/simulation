import math
import sys
import numpy as np
import astropy.units as u
from astropy.units import cds
from time import time
from numba import jit
from scipy.constants import c
from scipy.constants import k
from astropy.io import fits
from pynufft import NUFFT
from pathlib import Path
import finufft

from backClasses.visibility import Visibility
from backClasses.baseline import Baseline


class Interferometer:

    def __init__(self, antenna_route: str = None):
        self.antenna_pos = None
        self.sigma = 0 * cds.Jy
        self.antenna_number = None
        self.antenna_area = None
        self.gridded_vo = None
        self.gridded_image = None
        self.sky_image = None
        self.fft_image = None
        self.visibilities = None
        if antenna_route is not None:
            # assignment of antenna position
            self.read_antenna_config(antenna_route)
            # assignment of baseline
            self.compute_baselines()
            #self.baseline = Baseline(baseline, self.antenna_pos)
        else:
            raise ValueError('antenna route has not been initialized')

    def read_antenna_config(self, route: str):
        try:
            # extract content from the route file
            file_observatory = np.loadtxt(route, skiprows=3, usecols=(0, 1, 2, 3))
            # get other parameters
            antenna_position = file_observatory[:, :3]
            antenna_number = len(file_observatory)
            antenna_radius = np.sum(file_observatory[:, 3]) / (2 * antenna_number)
            antenna_radius *= u.m
            # save parameters
            self.antenna_area = np.pi * (antenna_radius ** 2)
            self.antenna_number = antenna_number
            self.antenna_pos = antenna_position * u.m
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
        baseline = ant_set1 - ant_set2
        # take off the invalid [0,0,0] values of baseline (the distance of an antenna from itself)
        rows = np.where(np.sum(baseline, axis=1) == 0)
        baseline = np.delete(baseline, rows, 0)

        # Solo para pruebas unitarias
        #baseline = baseline[:1, :] # hay que remover esto <------------------

        # transpose the baselines array
        baseline = baseline.transpose()
        self.baseline = Baseline(baseline, self.antenna_pos)

    def compute_hour_angle(self, ha_start: float = None, ha_end: float = None, dt: int = None):
        # ha_rad[] vector with hour angle samplings in radians!!!
        # ha_star, ha_end in Hours (for instance -1 hora, 1 hora)
        # dt sampling interval in seconds (for instance 60 seconds)
        # check hour angle limits

        # take the units, to make the values dimensionless
        ha_unit = ha_start.unit
        ha_start = ha_start.value
        ha_end = ha_end.value
        dt = dt.value
        #
        ha_rad = np.arange(ha_start, ha_end, dt)
        ha_rad = ha_rad * ha_unit
        return ha_rad

    def create_uv_position(self, lat_rad: float = None, decl_rad: float = None,
                           ha: np.ndarray = None, obs_freq: float = None):
        # obs_freq is observing frequency in Hz (continuum ??)
        # telescope_lat in degrees. Example -23
        # source_decl Source declination in degrees. For instance 18
        '''decl_rad = 0
        lat = 0
        if telescope_lat is not None and source_decl is not None:
            decl_rad = source_decl * np.pi / 180
            lat = (-90 + telescope_lat) * np.pi / 180
        '''

        # Rotate around x to rise w -Dec degrees
        R2 = np.array([[1, 0, 0],
                       [0, np.cos(((np.pi / 2) * u.rad) - decl_rad), -np.sin(((np.pi / 2) * u.rad) - decl_rad)],
                       [0, np.sin(((np.pi / 2) * u.rad) - decl_rad), np.cos(((np.pi / 2) * u.rad) - decl_rad)]])

        # Rotate around x to correct for telescope latitude
        R4 = np.array([[1, 0, 0],
                       [0, np.cos(lat_rad), -np.sin(lat_rad)],
                       [0, np.sin(lat_rad), np.cos(lat_rad)]])

        # create a R3 for every HA
        R3 = np.array([np.cos(ha), -np.sin(ha), np.zeros([len(ha)]),
                       np.sin(ha), np.cos(ha), np.zeros([len(ha)]),
                       np.zeros([len(ha)]), np.zeros([len(ha)]), np.ones([len(ha)])])

        # re-order the array from "a matriz of array" to "array of matrix"
        R3 = np.swapaxes(R3, 1, 0)
        lambda_num = c * (u.m / u.s) / obs_freq
        # operate every element in the array
        uvw = np.apply_along_axis(self.apply_rotation_matrix, 1, R3, R4, R2, lambda_num)
        # re-order the visibilities to be easier to scattering, transposing to get:
        #   u positions array, v positions array and w positions array
        m, n, l = np.shape(uvw)
        uvw = uvw.transpose([1, 0, 2]).reshape(n, m * l)
        imagesize = len(self.sky_image)
        self.visibilities = Visibility(uvw, obs_freq, imagesize)

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

    def write_image(self, data, route: str = None, name: str = None):
        print('------>', route, '-------->', name)
        '''data_folder = Path("source_data/text_files/")
        file_to_open = data_folder / "raw_data.txt"
        f = open(file_to_open)
        print(f.read())'''
        fits.writeto(route + name + '.fits', data)

    def transform_fft(self):
        # transform with fast fourier transform in 2 dimensions
        complex_sky_image = self.sky_image + 1j * np.zeros_like(self.sky_image)
        fft_image = np.fft.fft2(complex_sky_image)
        self.fft_image = fft_image

    def transform_nufft(self):
        # Transform from lambda to 2pi range
        u = self.visibilities.UVW[0] * np.pi / self.visibilities.max_uv_coordinate
        v = self.visibilities.UVW[1] * np.pi / self.visibilities.max_uv_coordinate

        complex_sky_image = self.sky_image + 1j * np.zeros_like(self.sky_image)

        # the 2D transform outputs f array of shape (N1, N2)
        f = finufft.nufft2d2(v, u, complex_sky_image, eps=1e-12)
        self.visibilities.set_value(f * cds.Jy)

    def inverse_transform_fft(self, usefft):
        dirty_image = np.fft.ifftshift(self.gridded_vo)
        dirty_image = np.fft.ifft2(dirty_image)
        if not usefft:
            dirty_image = np.fft.fftshift(dirty_image)
        self.dirty_image = dirty_image

    def bilinear_interpolation(self):
        u = self.visibilities.UVW[0]
        v = self.visibilities.UVW[1]
        # other variables
        N = max(len(self.fft_image[0]), len(self.fft_image))
        deltau = - self.visibilities.delatu

        # find the list of j and i equivalent to one of the four points, Q11 points
        list_j = np.floor(u / deltau) + ((N - 1) / 2)  # horizontal
        list_i = np.floor(v / deltau) + ((N - 1) / 2)  # vertical

        # then is applied  the linear interpolation
        l_alpha = (u - ((list_j - ((N - 1) / 2)) * deltau)) / deltau
        l_beta = (v - ((list_i - ((N - 1) / 2)) * deltau)) / deltau
        uv_value = np.apply_along_axis(self.calculate_uv_value, 0, [list_j, list_i, l_alpha, l_beta])

        # values are created for the visibilities object
        self.visibilities.set_value(uv_value * cds.Jy)

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
        self.visibilities.set_value(result * cds.Jy)

    @jit(forceobj=True)
    def test_function_2(self, uvw, sky_image, j_set, i_set):
        value = (uvw[0] * i_set) + (uvw[1] * j_set)
        e = np.cos(value) + 1j * np.sin(value)
        new_visibility = sky_image * (e)
        new_visibility = np.sum(new_visibility)
        return new_visibility

    def gridder(self, imagesize, scheme, robust):
        # vis numpy array with -> [u v w Re Im We f]
        # scheme = 1, robust = 0
        # return gridded image, deltau, max uv, delta x
        uvw = self.visibilities.UVW.transpose()
        uvw = uvw[:, :2]
        uv_value = self.visibilities.uv_value + self.noise

        z = len(uv_value)
        pos_u_index = np.zeros(z).astype(int)
        pos_v_index = np.zeros(z).astype(int)
        deltau = self.visibilities.deltau
        deltav = self.visibilities.deltav

        # Gridding weights if ROBUST OR UNIFORM
        gridded_weights = np.zeros((imagesize, imagesize))
        for k in range(0, z):
            j = int(np.round(uvw[k][0] / deltau) + imagesize / 2)
            i = int(np.round(uvw[k][1] / deltav) + imagesize / 2)
            pos_v_index[k] = i
            pos_u_index[k] = j
            gridded_weights[i][j] += self.visibilities.weight[k]

        # Selecting a scheme for weights
        self.scheme_weights(gridded_weights, z, pos_u_index, pos_v_index, scheme, robust)

        # Gridding visibilities and weights with a scheme
        gridded_weights = np.zeros((imagesize, imagesize)) * (cds.Jy / cds.Jy)
        gridded_Vo = np.zeros((imagesize, imagesize)) + 1.0j * np.zeros((imagesize, imagesize)) * cds.Jy

        for k in range(0, z):
            j = int(np.round(uvw[k][0] / deltau) + imagesize / 2)
            i = int(np.round(uvw[k][1] / deltav) + imagesize / 2)
            gridded_weights[i][j] += self.visibilities.weight[k]
            gridded_Vo[i][j] += self.visibilities.weight[k] * uv_value[k]

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
            if scheme == 0:
                self.visibilities.weight[k] = self.visibilities.weight[k]
            elif scheme == 1:
                self.visibilities.weight[k] /= gridded_weights[pos_v_index[k]][pos_u_index[k]]
            elif scheme == 2:
                self.visibilities.weight[k] *= np.sqrt(
                    self.visibilities.UVW[0][k] ** 2 + self.visibilities.UVW[1][k] ** 2)
            else:
                f2 = (5.0 * np.power(10.0, -robust)) ** 2 / average_weights
                self.visibilities.weight[k] /= (1.0 + gridded_weights[pos_v_index[k]][pos_u_index[k]] * f2)

    def run(self, telescope_lat: float = None, source_decl: float = None, ha_start: float = None,
            ha_end: float = None, dt: int = None, obs_freq: float = None, usefft: bool = None):
        # toma: la imagen de entrada(sky_image), el objeto visibilities(self)

        # convert units
        telescope_lat = telescope_lat.to(u.rad)
        source_decl = source_decl.to(u.rad)
        ha_start = ha_start.to(u.rad)
        ha_end = ha_end.to(u.rad)
        dt = dt.to(u.deg)
        dt = dt.to(u.rad)
        obs_freq = obs_freq.to(u.Hz)

        # uv positions
        ha = self.compute_hour_angle(ha_start, ha_end, dt)
        self.create_uv_position(telescope_lat, source_decl, ha, obs_freq)
        # values of visibilities
        if usefft:
            # image with fft
            self.transform_fft()
            self.fft_image = np.fft.fftshift(self.fft_image)
            self.bilinear_interpolation()
            #self.fourier_series()

        else:
            self.transform_nufft()
            # self.fourier_series()

        # add noise
        if rms_parameters is not None:
            self.get_noise_level(rms_parameters[0], rms_parameters[1], rms_parameters[2])
            self.add_noise()
            if rms_parameters[1] == 0 or rms_parameters[2] == 0:
                raise ValueError('invalid parameter, integration time and bandwidth must not be 0')

        self.gridder(len(self.fft_image), 1, 0)
        self.inverse_transform_fft()

    def get_noise_level(self, system_temperature, integration_time, bandwidth):

        sqrt_root = self.antenna_number * (self.antenna_number - 1)
        sqrt_root = math.sqrt(sqrt_root * integration_time * bandwidth)
        self.sigma = 2 * (k * u.J / u.K) * system_temperature / (self.antenna_area * sqrt_root)
        self.sigma = self.sigma.to(u.Jy)

    def add_noise(self):
        shape_values = self.visibilities.uv_value.shape
        real_noise = np.random.normal(0, self.sigma.value, shape_values) * self.sigma.unit
        imaginary_noise = np.random.normal(0, self.sigma.value, shape_values) * self.sigma.unit * 1j
        self.noise = real_noise + imaginary_noise
