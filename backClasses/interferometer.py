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
from backClasses.fouriertransformer import FT


class Interferometer:

    def __init__(self, antenna_route: str = None):
        self.antenna_pos = None
        self.sigma = 0 * cds.Jy
        self.antenna_number = None
        self.antenna_area = None
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

    @jit(forceobj=True)
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

    @jit
    def bilinear_interpolation(self):
        u = self.visibilities.UVW[0]
        v = self.visibilities.UVW[1]
        # other variables
        N = max(len(self.fft_image[0]), len(self.fft_image))
        deltau = - self.visibilities.deltau

        # find the list of j and i equivalent to one of the four points, Q11 points
        list_j = np.floor(u / deltau) + ((N - 1) / 2)  # horizontal
        list_i = np.floor(v / deltau) + ((N - 1) / 2)  # vertical

        # then is applied  the linear interpolation
        l_alpha = (u - ((list_j - ((N - 1) / 2)) * deltau)) / deltau
        l_beta = (v - ((list_i - ((N - 1) / 2)) * deltau)) / deltau
        uv_value = np.apply_along_axis(self.calculate_uv_value, 0, [list_j, list_i, l_alpha, l_beta])

        # values are created for the visibilities object
        self.visibilities.set_value(uv_value * cds.Jy)

    @jit
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

    #@jit(forceobj=True)

    def run(self, telescope_lat: float = None, source_decl: float = None, ha_start: float = None,
            ha_end: float = None, dt: int = None, obs_freq: float = None, usefft: bool = None):
        # toma: la imagen de entrada(sky_image), el objeto visibilities(self)

        ft = FT()
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
            self.fft_image = ft.transform_fft(self.sky_image) # self.transform_fft()
            self.bilinear_interpolation()
            #self.fourier_series()
            # delatu
        else:
            u_pos = self.visibilities.UVW[0]
            v_pos = self.visibilities.UVW[1]
            max_uv = self.visibilities.max_uv_coordinate
            self.visibilities.set_value(ft.transform_nufft(self.sky_image, u_pos, v_pos, max_uv))


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
