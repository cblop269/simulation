import sys
import astropy.units as u
import numpy as np
from astropy.io import fits


class FileManager:

    def __init__(self):
        value = 0

    def read_image(self, route: str = None):
        try:
            with fits.open(route) as image:
                print(image.info())
                image_data = image[0].data
            self.sky_image = image_data
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            sys.exit()

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

    def write_image(self, data, route: str = None, name: str = None):
        if type(data) is u.quantity.Quantity:
            data = data.value

        real = data.real
        imag = data.imag
        hdul = fits.HDUList([fits.PrimaryHDU(real), fits.ImageHDU(imag)]) #hdul = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(real), fits.ImageHDU(imag)])

        hdul.writeto(route + '/' + name + '.fits')

    def write_array(self, data, route: str = None, name: str = None):
        # header
        # delta u: delta v: number of antenna: system temperature: integration time:
        # value real: value imaginary: weight: noise real: noise imaginary: frequency:
        if type(data) is u.quantity.Quantity:
            data = data.value
        np.savetxt(route + '/' + name + '.npy', data)

    def write_visibilities(self, interferometer, route: str = None, name: str = None):
        # header
        # delta u: delta v: number of antenna: system temperature: integration time:
        # value real: value imaginary: weight: noise real: noise imaginary: frequency:
        deltau = 'Delta U:' + str(interferometer.deltau)
        deltav = 'Delta V:' + str(interferometer.deltav)
        antenna_number = 'N° Antenna:' + str(interferometer.antenna_number)
        syst_temp = 'System Temperature:' + str(interferometer.antenna_number)
        integration_time = 'Integration Time:' + str(interferometer.antenna_number)
        header1 = '#' + deltau + '\t' + deltav + '\t' + antenna_number + '\t' + syst_temp + '\t' + integration_time\
                  + '\n'
        header2 = '#' + 'U' + '\t' + 'V' + '\t' + 'W' + '\t' 'Value Re' + '\t' + 'Value Im' + '\t' + 'Weight'\
                  + '\t' + 'Noise Re' + '\t' + 'Noise Im' + '\t' + 'Freq' + '\n'
        # open file
        myfile = open(route + '/' + name + '.txt', 'w')
        myfile.write(header1)
        myfile.write(header2)
        for z in range(0, len(interferometer.visibilities.value)):
            u = str(interferometer.visibilities.UVW[0][z])
            v = str(interferometer.visibilities.UVW[1][z])
            w = str(interferometer.visibilities.UVW[2][z])
            valueR = str((interferometer.visibilities.value.value[z]).real)
            valueI = str((interferometer.visibilities.value.value[z]).imag)
            weight = str(1)
            noiseR = str((interferometer.noise.value[z]).real)
            noiseI = str((interferometer.noise.value[z]).imag)
            frequency = str(interferometer.visibilities.frequency)
            line = u + '\t' + v + '\t' + w + '\t' + valueR + '\t' + valueI + '\t' + weight + '\t' + noiseR + '\t'\
                   + noiseI + '\t' + frequency + '\n'
            myfile.write(line)
        myfile.close()
