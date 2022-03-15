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

from backClasses.visibility import Visibility
from backClasses.baseline import Baseline
from backClasses.fouriertransformer import FT

class Interferometer:
    sigma = 0 * cds.Jy
    sky_image = None
    fft_image = None
    dirty_image = None
    visibilities = None
    def __init__(self, antenna_route: str = None):
        self.baseline = None
        self.noise = None
        if antenna_route is not None:
            # assignment of antenna parameters
            antenna_config = self.read_antenna_config(antenna_route)
            # assignment of baseline
            self.compute_baselines(antenna_config)
        else:
            raise ValueError('antenna route has not been initialized')

    def read_antenna_config(self, route: str):
        """
        Read a configuration of antennas with extension '.cfg', and return the data
        :param route: The route of the antenna configuration file (.cfg)
        :return: The position in y, v, w and the diameter od dish
        """
        try:
            # extract content from the route file
            file_observatory = np.loadtxt(route, skiprows=3, usecols=(0, 1, 2, 3))
            return file_observatory
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))


    @jit(forceobj=True)
    def compute_baselines(self, antenna_config:np.ndarray=None):
        """
        Compute all the baselines for the antenna positions
        :return: The baseline results
        """
        antenna_pos = antenna_config[:, :3] * u.m
        # subtract all the positions
        ant_set1 = np.repeat(antenna_pos, repeats=len(antenna_pos), axis=0)
        ant_set2 = np.concatenate([antenna_pos] * len(antenna_pos), axis=0)
        baseline = ant_set1 - ant_set2
        # take off the invalid [0,0,0] values of baseline (the distance of an antenna from itself)
        rows = np.where(np.sum(baseline, axis=1) == 0)
        baseline = np.delete(baseline, rows, 0)

        # Only for unitery tests
        #baseline = baseline[:1, :] # remove this <------------------

        # transpose the baselines array
        baseline = baseline.transpose()
        antenna_radius = np.sum(antenna_config[:, 3]) * u.m / (2 * len(antenna_config))
        self.baseline = Baseline(baseline, antenna_pos, antenna_radius)

    def __compute_hour_angle(self, ha_start: float = None, ha_end: float = None, dt: int = None):
        """
        Calculate an array with all the Hour angles along the sample time
        :param ha_start: Hour angle of start
        :param ha_end: Hour angle of end
        :param dt: The sample interval
        :return: The array with all the Hour angles
        """
        # take the units, to make the values dimensionless
        ha_unit = ha_start.unit
        ha_start = ha_start.value
        ha_end = ha_end.value
        dt = dt.value
        #
        ha_rad = np.arange(ha_start, ha_end, dt)
        ha_rad = ha_rad * ha_unit
        return ha_rad

    def __create_uv_position(self, lat_rad: float = None, decl_rad: float = None, ha: np.ndarray = None,
                             obs_freq: float = None):
        """
        Function that calculates the uvw positions
        :param lat_rad: latitude of the telescope in rad units
        :param decl_rad: Declination of the source in rad units
        :param ha: array of Hour angles in rad units
        :param obs_freq: Observation frequency
        :return The uvw positions
        """
        # Rotate around x to rise w -Dec degrees
        R2 = np.array([[1, 0, 0],
                       [0, np.cos(((np.pi / 2) * u.rad) - decl_rad), -np.sin(((np.pi / 2) * u.rad) - decl_rad)],
                       [0, np.sin(((np.pi / 2) * u.rad) - decl_rad), np.cos(((np.pi / 2) * u.rad) - decl_rad)]])

        # Rotate around x to correct for telescope latitude
        R4 = np.array([[1, 0, 0],
                       [0, np.cos(lat_rad - ((np.pi / 2) * u.rad)), -np.sin(lat_rad - ((np.pi / 2) * u.rad))],
                       [0, np.sin(lat_rad - ((np.pi / 2) * u.rad)), np.cos(lat_rad - ((np.pi / 2) * u.rad))]])

        # create a R3 for every HA
        R3 = np.array([np.cos(ha), -np.sin(ha), np.zeros([len(ha)]),
                       np.sin(ha), np.cos(ha), np.zeros([len(ha)]),
                       np.zeros([len(ha)]), np.zeros([len(ha)]), np.ones([len(ha)])])

        # re-order the array from "a matriz of array" to "array of matrix"
        R3 = np.swapaxes(R3, 1, 0)
        lambda_num = c * (u.m / u.s) / obs_freq
        # operate every element in the array
        uvw = np.apply_along_axis(self.__apply_rotation_matrix, 1, R3, R4, R2, lambda_num)
        # re-order the visibilities to be easier to scattering, transposing to get:
        #   u positions array, v positions array and w positions array
        m, n, l = np.shape(uvw)
        uvw = uvw.transpose([1, 0, 2]).reshape(n, m * l)
        return uvw

    def __apply_rotation_matrix(self, matrix3: np.ndarray, matrix4: np.matrix, matrix2: np.matrix,
                              lambda_num: float):
        """
        Function that calculates the projection application of matrix
        :param matrix3: Rotation matrix of HA
        :param matrix4: Rotation matrix of latitude
        :param matrix2: Rotation matrix of declination
        :param lambda_num: The wavelength
        :return The uv position of the projection
        """
        uvw = matrix4 @ matrix3.reshape(3, 3) @ matrix2
        uvw = uvw.transpose()
        uvw = (uvw @ self.baseline.baselines) / lambda_num
        return np.asarray(uvw)

    def read_image(self, route: str = None):
        """
        Read a file a '.fits' file with the sky image, and define it as attribute
        :param route: The route of the sky image file (.cfg)
        """
        try:
            with fits.open(route) as image:
                image_data = image[0].data
            self.sky_image = image_data
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))

    @jit
    def bilinear_interpolation(self, u_pos, v_pos):
        """
        Function that aplies the bilinear interpolation to the visibilities
        :param u_pos: the u positions
        :param v_pos: the v positions
        :return The uv position interpolated
        """
        # other variables
        N = max(len(self.fft_image[0]), len(self.fft_image))
        deltau = - self.visibilities.deltau

        # find the list of j and i equivalent to one of the four points, Q11 points
        list_j = np.floor(u_pos / deltau) + ((N - 1) / 2)  # horizontal
        list_i = np.floor(v_pos / deltau) + ((N - 1) / 2)  # vertical

        # then is applied  the linear interpolation
        l_alpha = (u_pos - ((list_j - ((N - 1) / 2)) * deltau)) / deltau
        l_beta = (v_pos - ((list_i - ((N - 1) / 2)) * deltau)) / deltau
        uv_value = np.apply_along_axis(self.__interpolate_uv_value, 0, [list_j, list_i, l_alpha, l_beta])

        return uv_value * cds.Jy

    @jit
    def __interpolate_uv_value(self, data: np.ndarray = None):
        """
        Function that calculates the bilinear interpolation of an uv value
        :param data: current u position = data[0]
                     current v position = data[1]
                     current alpha = data[2]
                     current beta = data[3])
        :return The uv value interpolated
        """
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
        """
        Function that makes the interometer run, and save the visibilities as an attribute
        :param telescope_lat: latitude of the telescope
        :param source_decl: Declination of the source
        :param ha_start: Hour angle of start
        :param ha_end: Hour angle of end
        :param dt: The sample interval
        :param obs_freq: Observation frequency
        :param usefft: Boolean value, if True fft will be used, if False nufft will be used
        """
        # create the object fourier transformer
        ft = FT()
        # convert units
        telescope_lat = telescope_lat.to(u.rad)
        source_decl = source_decl.to(u.rad)
        ha_start = ha_start.to(u.rad)
        ha_end = ha_end.to(u.rad)
        dt = dt.to(u.deg)
        dt = dt.to(u.rad)
        obs_freq = obs_freq.to(u.Hz)

        # Calculate the uv positions
        ha = self.__compute_hour_angle(ha_start, ha_end, dt)
        uvw = self.__create_uv_position(telescope_lat, source_decl, ha, obs_freq)
        imagesize = len(self.sky_image)
        self.visibilities = Visibility(uvw, obs_freq, imagesize)

        # Calculate the values of visibilities
        u_pos = self.visibilities.UVW[0]
        v_pos = self.visibilities.UVW[1]
        if usefft:
            # with fft
            self.fft_image = ft.transform_fft(self.sky_image) # self.transform_fft()
            uv_value = self.bilinear_interpolation(u_pos, v_pos)
            self.visibilities.set_value(uv_value)
        else:
            # with nufft
            max_uv = self.visibilities.max_uv_coordinate
            uv_value = ft.transform_nufft(self.sky_image, u_pos, v_pos, max_uv)
            self.visibilities.set_value(uv_value)


    def get_noise_level(self, sys_T, intg_t, bw):
        """
        Recalculate the sigma atributte
        :param sys_T: system_temperature
        :param intg_t: integration_time
        :param bw: Bandwidth
        """
        baseline = self.baseline
        sqrt_root = baseline.antenna_number * (baseline.antenna_number - 1)
        sqrt_root = math.sqrt(sqrt_root * intg_t * bw)
        self.sigma = 2 * (k * u.J / u.K) * sys_T / (baseline.antenna_area * sqrt_root)
        self.sigma = self.sigma.to(u.Jy)

    def add_noise(self):
        """
        create an array of noise for every visibility
        """
        shape_values = self.visibilities.uv_value.shape
        real_noise = np.random.normal(0, self.sigma.value, shape_values) * self.sigma.unit
        imaginary_noise = np.random.normal(0, self.sigma.value, shape_values) * self.sigma.unit * 1j
        self.noise = real_noise + imaginary_noise

    def set_dirty_image(self, dirty_image:np.ndarray = None):
        self.dirty_image = dirty_image
