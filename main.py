from tkinter import *
from tkinter.filedialog import askopenfilename
import numpy as np

from backClasses.interferometer import Interferometer
from frontClasses.testFrontv2 import Window
from frontClasses.frontWindow import App

if __name__ == '__main__':
    # interferometer = Interferometer('/home/seba/Desktop/hola.txt', -23.0234, -40, -6, 6, 600, 90e9)
    # window = App()
    # window.mainloop()
    window = Window()
    window.mainloop()
