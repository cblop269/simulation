import tkinter as tk
from pathlib import Path
from tkinter.filedialog import askopenfilename
import numpy as np
import astropy.units as u
from astropy.units import cds
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from numba import jit
from mpl_toolkits.axes_grid1 import make_axes_locatable
from frontClasses.plotter import Plotter
from backClasses.interferometer import Interferometer
from backClasses.filemanager import FileManager
from backClasses.gridder import Gridder
from backClasses.fouriertransformer import FT

class FrameD(tk.Frame):

    def __init__(self, notebook):#, interferometer):
        super().__init__(notebook)
        self.config(bg="gray26")
        #
        self.writer = FileManager()
        #   Frames
        self.frame5 = tk.Frame(self, bg="gray30")
        self.frame6 = tk.Frame(self, bg="gray30")
        self.frame7 = tk.Frame(self, bg="gray30")
        self.frame5.grid(row=0, column=1, pady=20, padx=10)
        self.frame6.grid(row=0, column=2, pady=20, padx=10)
        self.frame7.grid(row=0, column=3, pady=20, padx=10)

        self.set_labels()
        self.set_variables()
        self.set_input_entries()
        self.set_buttons(None, None)

    def set_labels(self):
        #   labels
        self.lbl_export_uv = tk.Label(self.frame5, text="Export Visibilities", bg="gray30", fg="white")
        self.lbl_export_gi = tk.Label(self.frame6, text="Export Grided Image", bg="gray30", fg="white")
        self.lbl_export_di = tk.Label(self.frame7, text="Export Dirty Image", bg="gray30", fg="white")
        self.lbl_route_uv = tk.Label(self.frame5, text=str(Path.home() / "Downloads"), bg="gray30", fg="white")
        self.lbl_route_gi = tk.Label(self.frame6, text=str(Path.home() / "Downloads"), bg="gray30", fg="white")
        self.lbl_route_di = tk.Label(self.frame7, text=str(Path.home() / "Downloads"), bg="gray30", fg="white")
        self.lbl_export_uv.grid(row=0, column=0, pady=20, padx=20)
        self.lbl_route_uv.grid(row=1, column=0, pady=0, padx=20)
        self.lbl_export_gi.grid(row=0, column=0, pady=20, padx=20)
        self.lbl_route_gi.grid(row=1, column=0, pady=0, padx=20)
        self.lbl_export_di.grid(row=0, column=0, pady=20, padx=20)
        self.lbl_route_di.grid(row=1, column=0, pady=0, padx=20)

    def set_variables(self):
        #   Variables
        self.default_export_uv = tk.StringVar()
        self.default_export_gi = tk.StringVar()
        self.default_export_di = tk.StringVar()
        self.default_export_uv.set('Visibilities')
        self.default_export_gi.set('GriddedImage')
        self.default_export_di.set('DirtyImage')

    def set_input_entries(self):
        #	Input
        self.imput_export_uv = tk.Entry(self.frame5, width=15, justify='right')
        self.imput_export_uv.config(textvariable=self.default_export_uv)
        self.imput_export_gi = tk.Entry(self.frame6, width=15, justify='right')
        self.imput_export_gi.config(textvariable=self.default_export_gi)
        self.imput_export_di = tk.Entry(self.frame7, width=15, justify='right')
        self.imput_export_di.config(textvariable=self.default_export_di)
        self.imput_export_uv.grid(row=2, column=0, pady=20, padx=20)
        self.imput_export_gi.grid(row=2, column=0, pady=20, padx=20)
        self.imput_export_di.grid(row=2, column=0, pady=20, padx=20)


    def set_buttons(self, interferometer, gridder):
        self.btn_uv = tk.Button(self.frame5, text='Export', width=10)
        self.btn_uv.config(state=tk.DISABLED, command=lambda: self.writer.write_visibilities(interferometer, self.label_route_uv.cget("text"), self.default_export_uv.get()))
        self.btn_gi = tk.Button(self.frame6, text='Export', width=10)
        self.btn_gi.config(state=tk.DISABLED, command=lambda: self.writer.write_array(gridder.gridded_vo, self.label_route_gi.cget("text"), self.default_export_gi.get()))
        self.btn_di = tk.Button(self.frame7, text='Export', width=10)
        self.btn_di.config(state=tk.DISABLED, command=lambda: self.writer.write_image(interferometer.dirty_image, self.label_route_di.cget("text"), self.default_export_di.get()))
        self.btn_uv.grid(row=3, column=0, pady=20, padx=20)
        self.btn_gi.grid(row=3, column=0, pady=20, padx=20)
        self.btn_di.grid(row=3, column=0, pady=20, padx=20)