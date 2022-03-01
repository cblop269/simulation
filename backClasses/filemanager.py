import sys
import astropy.units as u
import numpy as np
from astropy.io import fits


class FileManager:

    def write_image(self, data, route: str = None, name: str = None):
        """
        Write the dirty image in a fits file
        :param data: The data from the image
        :param route: The route in where the file was created
        :param name: name of the file
        """
        if type(data) is u.quantity.Quantity:
            data = data.value

        real = data.real
        imag = data.imag
        abs = np.abs(real + imag)
        #hdul = fits.HDUList([fits.PrimaryHDU(real), fits.ImageHDU(imag)])
        hdul = fits.HDUList([fits.PrimaryHDU(abs), fits.ImageHDU(real), fits.ImageHDU(imag)])

        hdul.writeto(route + '/' + name + '.fits')

    def write_array(self, data, route: str = None, name: str = None):
        """
        Write the numpy array with the gridded values in a .npy file
        :param data: The data with the values
        :param route: The route in where the file was created
        :param name: name of the file
        """
        # header
        # delta u: delta v: number of antenna: system temperature: integration time:
        # value real: value imaginary: weight: noise real: noise imaginary: frequency:
        if type(data) is u.quantity.Quantity:
            data = data.value
        np.savetxt(route + '/' + name + '.npy', data)

    def write_visibilities(self, interferometer, route: str = None, name: str = None):
        """
        Write a .txt file with the visibilities
        :param interferometer: The interferometer with the parameters
        :param route: The route in where the file was created
        :param name: name of the file
        """
        # header
        # delta u: delta v: number of antenna: system temperature: integration time:
        # value real: value imaginary: weight: noise real: noise imaginary: frequency:
        deltau = 'Delta U:' + str(interferometer.visibilities.deltau)
        deltav = 'Delta V:' + str(interferometer.visibilities.deltav)
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
        for z in range(0, len(interferometer.visibilities.uv_value)):
            u = str(interferometer.visibilities.UVW[0][z])
            v = str(interferometer.visibilities.UVW[1][z])
            w = str(interferometer.visibilities.UVW[2][z])
            valueR = str((interferometer.visibilities.uv_value.value[z]).real)
            valueI = str((interferometer.visibilities.uv_value.value[z]).imag)
            weight = str(interferometer.visibilities.weight[z])
            noiseR = str((interferometer.noise.value[z]).real)
            noiseI = str((interferometer.noise.value[z]).imag)
            frequency = str(interferometer.visibilities.frequency)
            line = u + '\t' + v + '\t' + w + '\t' + valueR + '\t' + valueI + '\t' + weight + '\t' + noiseR + '\t'\
                   + noiseI + '\t' + frequency + '\n'
            myfile.write(line)
        myfile.close()
