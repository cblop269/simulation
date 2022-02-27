import finufft
import numpy as np
from astropy.units import cds
from numba import jit


class FT:
    def __init__(self):
        value = 0

    def transform_fft(self, sky_image):
        # transform with fast fourier transform in 2 dimensions
        complex_sky_image = sky_image + 1j * np.zeros_like(sky_image)
        fft_image = np.fft.fft2(complex_sky_image)
        fft_image = np.fft.fftshift(fft_image)
        return fft_image

    def transform_nufft(self, sky_image, u_pos, v_pos, max_uv):
        # Transform from lambda to 2pi range
        u = u_pos * np.pi / max_uv
        v = v_pos * np.pi / max_uv

        complex_sky_image = sky_image + 1j * np.zeros_like(sky_image)
        # the 2D transform outputs f array of shape (N1, N2)
        f = finufft.nufft2d2(v, u, complex_sky_image, eps=1e-12)
        return f * cds.Jy

    def transform_ifft(self, usefft, gridded_vo):
        dirty_image = np.fft.ifftshift(gridded_vo)
        dirty_image = np.fft.ifft2(dirty_image)
        if not usefft:
            dirty_image = np.fft.fftshift(dirty_image)
        return dirty_image

    @jit  # (forceobj=True)
    def transform_fs(self, sky_image, u_pos, v_pos, max_uv):
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
        value = (uvw[0] * j_set) + (uvw[1] * i_set)
        e = np.cos(value) + 1j * np.sin(value)
        new_visibility = sky_image * e
        new_visibility = np.sum(new_visibility)
        return new_visibility