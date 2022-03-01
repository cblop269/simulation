import finufft
import numpy as np
from astropy.units import cds
from numba import jit


class FT:

    def transform_fft(self, sky_image):
        """
        Apliccation of the fft
        :param sky_image: the input image of the object
        :return A array with the transformed values of the image
        """
        # transform with fast fourier transform in 2 dimensions
        complex_sky_image = sky_image + 1j * np.zeros_like(sky_image)
        fft_image = np.fft.fft2(complex_sky_image)
        fft_image = np.fft.fftshift(fft_image)
        return fft_image

    def transform_nufft(self, sky_image, u_pos, v_pos, max_uv):
        """
        Apliccation of the nufft
        :param sky_image: the input image of the object
        :param u_pos: the u positions
        :param v_pos: the v positions
        :param max_uv: the max uv position value
        :return A array with the transformed values of the image
        """
        # Transform from lambda to 2pi range
        u = u_pos * np.pi / max_uv
        v = v_pos * np.pi / max_uv

        complex_sky_image = sky_image + 1j * np.zeros_like(sky_image)
        # the 2D transform outputs f array of shape (N1, N2)
        f = finufft.nufft2d2(v, u, complex_sky_image, eps=1e-12)
        return f * cds.Jy

    def transform_ifft(self, usefft, gridded_vo):
        """
        Apliccation of the inverse fft
        :param usefft: A boolean value to decide if the shift of the image will be used or not
        :param gridded_vo: The values of the gridden visibilities
        :return The dirty image reconstructed
        """
        dirty_image = np.fft.ifftshift(gridded_vo)
        dirty_image = np.fft.ifft2(dirty_image)
        if not usefft:
            dirty_image = np.fft.fftshift(dirty_image)
        return dirty_image

    @jit  # (forceobj=True)
    def transform_fs(self, sky_image, u_pos, v_pos, max_uv):
        """
        Aplication of the Fourie Series to get the uv values
        :param sky_image: the input image of the object
        :param u_pos: the u positions
        :param v_pos: the v positions
        :param max_uv: the max uv position value
        :return A array with the transformed values of the image
        """
        m, n = np.shape(sky_image)

        uv = np.concatenate((u_pos, v_pos), axis=0)
        uv = uv.transpose()

        # convert to radians
        uv = uv * np.pi / max_uv
        #
        sky_image = sky_image.reshape(m * n)
        x = np.linspace(- n / 2, (n / 2) - 1, n)
        y = np.linspace(- m / 2, (m / 2) - 1, m)
        xg, yg = np.meshgrid(x, y)
        xg = xg.reshape(m * n)
        yg = yg.reshape(m * n)
        #
        result = np.apply_along_axis(self.calculate_uv_value_FS, 1, uv, sky_image, xg, yg)
        return result * cds.Jy

    @jit
    def calculate_uv_value_FS(self, uvw, sky_image, j_set, i_set):
        """
        Sub-function of the fourier series, calculate a summation for one value
        :param uvw: The uvw positions
        :param sky_image: the input image of the object
        :param j_set: the x positions in the image
        :param i_set: the y positions in the image
        :return A uv value with the Fourier series
        """
        value = (uvw[0] * j_set) + (uvw[1] * i_set)
        e = np.cos(value) + 1j * np.sin(value)
        new_visibility = sky_image * e
        new_visibility = np.sum(new_visibility)
        return new_visibility