
import tkinter as tk
from tkinter import ttk
#from back import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import os


class App(tk.Tk):#tk.Frame):
		
	def __init__(self, master=None):
		super().__init__()
		
		self.minsize(1000, 400)#810, 500)
		self.my_menu = tk.Menu(self)
		self.m1 = tk.Menu(self.my_menu)
		self.m2 = tk.Menu(self.my_menu)
		self.m3 = tk.Menu(self.my_menu)
		self.m1.add_command(label="Open observatory", command= self.open_file_O)
		self.m1.add_command(label="Open source", command= self.open_file_S)
		self.m2.add_command(label="option 2-a", command= lambda:print("HOLA"))
		self.m2.add_command(label="option 2-b", command= lambda:print("HOLA"))
		self.m3.add_command(label="option 3-a", command= lambda:print("HOLA"))
		self.m3.add_command(label="option 3-b", command= lambda:print("HOLA"))
		self.my_menu.add_cascade(label="File", menu=self.m1)
		self.my_menu.add_cascade(label="option 2", menu=self.m2)
		self.my_menu.add_cascade(label="option 3", menu=self.m3)
		self.config(menu=self.my_menu)

		self.style=ttk.Style()
		self.style.theme_use("alt")
		
		self.master = master
		self.window = tk.Frame(self, width=810, height=500)
		self.window.place(relwidth=1, relheight=1)
		self.create_widgets()
	
	def unCheck(self):
		if(self.current_algorythm.get()):
			self.nufft.set(True)
			self.fft.set(False)
			self.current_algorythm.set(False)
		else:
			self.nufft.set(False)
			self.fft.set(True)
			self.current_algorythm.set(True)
		print(self.nufft.get())
		print(self.fft.get())
	
	def read_data(self):
		print(self.latitude_str.get()) 
		print(self.longitude_str.get())
		print(self.rad_asc_str.get())
		print(self.rad_dec_str.get())
		print(self.sample_time.get())
		print(self.sample_interval.get())
		print(self.sample_number.get())
		print(self.hour_angle.get())
		print(self.current_algorythm.get())
		print(self.frecuency_list.get())
		prueba()
	
	def draw_graphic_f(self):
		self.canvas = FigureCanvasTkAgg(draw_vi(self.source_directory), master = self.frameE)
		self.canvas.draw()

		self.export_button.config(state = tk.NORMAL)
		self.toolbar = NavigationToolbar2Tk(self.canvas, self.frameE)
		self.toolbar.update()
		self.canvas._tkcanvas.place(relx=0, rely=0, relwidth=1, relheight=0.9)#pack()#side=tk.TOP, fill=tk.BOTH, expand=True)
		self.notebook.select(self.frameE)
		# frame A
		self.canvas_V = FigureCanvasTkAgg(draw_observatory(self.observatory_directory), master = self.frameA)
		self.canvas_V.draw()
		self.toolbar_V = NavigationToolbar2Tk(self.canvas_V, self.frameA)
		self.toolbar_V.update()
		self.canvas_V._tkcanvas.place(relx=0, rely=0, relwidth=0.6, relheight=0.9)

	def open_file_O(self):
		self.observatory_directory = tk.filedialog.askopenfilename()
		self.import_button_O.config(bg="orange")
		self.obserbatory_use = True
		if(self.source_use==True and self.obserbatory_use==True):
			self.run_button_out.config(state = tk.NORMAL)
			self.run_button.config(state = tk.NORMAL)
		self.canvas_O = FigureCanvasTkAgg(draw_observatory(self.observatory_directory), master = self.frameC)
		self.canvas_O.draw()
		self.toolbar_O = NavigationToolbar2Tk(self.canvas_O, self.frameC)
		self.toolbar_O.update()
		self.canvas_O._tkcanvas.place(relx=0, rely=0, relwidth=0.6, relheight=0.9)
		
	def open_file_S(self):
		self.source_directory = tk.filedialog.askopenfilename()
		print(self.source_directory)
		self.import_button_S.config(bg="orange")
		self.source_use = True
		if(self.source_use==True and self.obserbatory_use==True):
			self.run_button_out.config(state = tk.NORMAL)
			self.run_button.config(state = tk.NORMAL)
		self.canvas_S = FigureCanvasTkAgg(draw_source(self.source_directory), master = self.frameB)
		self.canvas_S.draw()
		self.toolbar_S = NavigationToolbar2Tk(self.canvas_S, self.frameB)
		self.toolbar_S.update()
		self.canvas_S._tkcanvas.place(relx=0, rely=0, relwidth=0.6, relheight=0.9)
		
	def create_widgets(self):
		
		#frame 0------------------------------------------------------------------------------
		self.frame0 = tk.Frame(self.window, padx=10, pady=10, bd=5)
		self.frame0.place(relwidth=1, relheight=1)
		
		#	creating components
		self.notebook = ttk.Notebook(self.frame0)
		self.frameA = tk.Frame(self.notebook)#, width=400, height=400)
		self.frameA.config(padx=10, pady=10)
		self.frameA.config(highlightbackground="black", highlightthickness=0.5)
		self.frameB = tk.Frame(self.notebook)#, width=440, height=440)
		self.frameB.config(padx=10, pady=10)#, bg='blue')
		self.frameC = tk.Frame(self.notebook)#, width=440, height=440)
		self.frameC.config(padx=10, pady=10)#, bg='red')
		self.frameD = tk.Frame(self.notebook)#, width=440, height=440)
		self.frameD.config(padx=10, pady=10)#, bg='red')
		self.frameE = tk.Frame(self.notebook)#, width=440, height=440)
		self.frameE.config(padx=10, pady=10)#, bg='red')
		#	adding conponentes to notebook and the frame
		self.notebook.add(self.frameE, text="Graphics")
		self.notebook.add(self.frameA, text="Visibility")
		self.notebook.add(self.frameB, text="Phantom")
		self.notebook.add(self.frameC, text="Observatory")
		self.notebook.add(self.frameD, text="Configuration")
		self.notebook.place(relwidth=1, relheight=1)
		self.notebook.select(self.frameD)
		
		
		#frame 1------------------------------------------------------------------------------
		self.frame1 = tk.Frame(self.frameD, highlightbackground="black", highlightthickness=1, )
		self.frame1.place(relx=0, rely=0, relwidth=0.48, relheight=0.48)
		#	labels
		self.label_position = tk.Label(self.frame1, text = 'Position')
		self.label_position.grid(row=0, column=0, columnspan=1, sticky=tk.W)
		self.label_observatory = tk.Label(self.frame1, text = 'observatory')
		self.label_observatory.grid(row=1, column=1, pady=10, sticky=tk.W)
		self.label_source = tk.Label(self.frame1, text = 'source')
		self.label_source.grid(row=2, column=1, pady=5, sticky=tk.W)
		#	values inside the boxes
		self.latitude_str = tk.StringVar(value = 'latitude') 
		self.longitude_str = tk.StringVar(value = 'longitude')
		self.rad_asc_str = tk.StringVar(value = 'rad asc') 
		self.rad_dec_str = tk.StringVar(value = 'rad dec') 
		#	boxes
		self.input_lat = tk.Entry(self.frame1, textvariable = self.latitude_str, width=10)
		self.input_lat.grid(row = 1, column = 2, pady=10, padx=2)
		self.input_long = tk.Entry(self.frame1, textvariable = self.longitude_str, width=10)
		self.input_long.grid(row = 1, column = 3, pady=10, padx=2)
		self.input_rad_asc = tk.Entry(self.frame1, textvariable = self.rad_asc_str)
		self.input_rad_asc.config(width=10)
		self.input_rad_asc.grid(row = 2, column = 2, pady=5, padx=2)
		self.input_rad_dec = tk.Entry(self.frame1, textvariable = self.rad_dec_str)
		self.input_rad_dec.config(width=10)
		self.input_rad_dec.grid(row = 2, column = 3, pady=5, padx=2)
		
		
		#frame 2--------------------------------------------------------------------------------
		self.frame2 = tk.Frame(self.frameD, highlightbackground="black", highlightthickness=0.5)
		self.frame2.place(relx=0,rely=0.52,relwidth=0.48, relheight=0.48)
		#	label
		self.label_observatory = tk.Label(self.frame2, text='Observation time')
		self.label_observatory.grid(row=0, column=0, pady=5, columnspan=2, sticky=tk.W)
		self.label_sample_frec = tk.Label(self.frame2, text='Sample frec')
		self.label_sample_frec.grid(row=2, column=1, pady=5, columnspan=2, sticky=tk.W)
		self.label_sample_n = tk.Label(self.frame2, text='Sample N')
		self.label_sample_n .grid(row=2, column=4, pady=10, sticky=tk.W)
		self.label_sample_time = tk.Label(self.frame2, text='Sampling time')
		self.label_sample_time.grid(row=1, column=1, pady=5, columnspan=2, sticky=tk.W)
		self.label_hour_angle = tk.Label(self.frame2, text='Hour A')
		self.label_hour_angle.grid(row=1, column=4, pady=5, sticky=tk.W)
		#	numerical values
		self.sample_time = tk.IntVar()
		self.sample_time.set(0.0000)
		self.sample_interval = tk.IntVar()
		self.sample_interval.set(0.0000)
		self.sample_number = tk.IntVar()
		self.sample_number.set(0.0000)
		self.hour_angle = tk.IntVar()
		self.hour_angle.set(0.0000)
		#	spinboxes
		#self.spin_s_time = tk.Spinbox(self.frame2, width=5, from_=0, to=10, increment=0.1)
		self.spin_s_time = tk.Entry(self.frame2, width=5)
		self.spin_s_time.config(textvariable = self.sample_time)
		self.spin_s_time.grid(row = 1, column = 3, pady=10, padx=5)
		
		#self.spin_interval =tk.Spinbox(self.frame2, width=5, from_=0, to=10, increment=0.1)
		self.spin_interval = tk.Entry(self.frame2, width=5)
		self.spin_interval.config(textvariable = self.sample_interval)
		self.spin_interval.grid(row = 2, column = 3, pady=5, padx=5)
		
		#self.spin_s_number =tk.Spinbox(self.frame2, width=5, from_=0, to=10, increment=0.1)
		self.spin_s_number = tk.Entry(self.frame2, width=5)
		self.spin_s_number.config(textvariable = self.sample_number)
		self.spin_s_number.grid(row = 2, column = 5, pady=5, padx=5)
		
		#self.spin_h_angle = tk.Spinbox(self.frame2, width=5, from_=0, to=10, increment=0.1)
		self.spin_h_angle = tk.Entry(self.frame2, width=5)
		self.spin_h_angle.config(textvariable = self.hour_angle)
		self.spin_h_angle.grid(row = 1, column = 5, pady=5, padx=5)
		
		
		#frame 3------------------------------------------------------------------------------
		self.frame3 = tk.Frame(self.frameD, highlightbackground="black", highlightthickness=0.5)
		self.frame3.place(relx=0.52, rely=0.35,relwidth=0.48, relheight=0.31)
		#	label
		self.label_algorythm = tk.Label(self.frame3, text='Algorythm')
		self.label_algorythm.grid(row = 0, column = 0, pady=10)
		self.label_frecuency = tk.Label(self.frame3, text='Frecuency')
		self.label_frecuency.grid(row = 1, column = 0, pady=5, padx=10)
		#	checks fft & nufft
		self.fft = tk.BooleanVar()
		self.nufft = tk.BooleanVar()
		self.current_algorythm = tk.BooleanVar()
		self.fft.set(True)
		self.nufft.set(False)
		self.current_algorythm.set(True)
		self.check_fft = ttk.Checkbutton(self.frame3, text='FFT', variable=self.fft)
		self.check_fft.config(onvalue=True, offvalue=False)
		self.check_fft.config(command = self.unCheck)
		self.check_fft.grid(row = 0, column = 1, pady=10, padx=10)
		self.check_nufft = ttk.Checkbutton(self.frame3, text='NUFFT', variable=self.nufft)
		self.check_nufft.config(onvalue=True, offvalue=False)
		self.check_nufft.config(command = self.unCheck)
		self.check_nufft.grid(row = 0, column = 2, pady=10, padx=10)
		#	box to frecuency
		self.frecuency_list = tk.Entry(self.frame3, textvariable = 'enter frecuency')
		self.frecuency_list.grid(row=1, column=1, pady=5, padx=10, columnspan=2)
		#self.frecuency_list.config(sticky=tk.W+tk.E, )
		
		
		#frame 4--------------------------------------------------------------------------------
		self.frame4 = tk.Frame(self.frameD, highlightbackground="black", highlightthickness=0.5)
		self.frame4.place(relx=0.52, rely=0.7, relwidth=0.48, relheight=0.30)
		#	run button
		self.run_button = tk.Button(self.frame4, text = 'Run', width=10)
		self.run_button.config(command = self.draw_graphic_f, state=tk.DISABLED)
		self.run_button.place(relx=0.25, rely=0.3, relwidth=0.2, relheight=0.4)#grid(row = 0, column = 0)
		#	export button
		self.export_button = tk.Button(self.frame4, text = 'Export', width=10)
		self.export_button.config(state = tk.DISABLED)
		self.export_button.place(relx=0.55, rely=0.3, relwidth=0.2, relheight=0.4)#grid(row = 0, column = 1)
		
		
		#frame 5--------------------------------------------------------------------------------
		self.frame5 = tk.Frame(self.frameD, highlightbackground="black", highlightthickness=0.5)
		self.frame5.place(relx=0.52, rely=0, relwidth=0.48, relheight=0.31)
		
		self.source_use = False
		self.obserbatory_use = False
		#self.source_use.set(False)
		#self.obserbatory_use.set(False)
		
		#	open file button source file
		self.import_button_S = tk.Button(self.frame5, text = 'Source file', width=10)
		self.import_button_S.config(command = self.open_file_S, bg="light blue")
		self.import_button_S.place(relx=0.12, rely=0.1, relwidth=0.33, relheight=0.6)
		#	open file button observatory file
		self.import_button_O = tk.Button(self.frame5, text = 'Observatory file', width=10)
		self.import_button_O.config(command = self.open_file_O, bg="light blue")
		self.import_button_O.place(relx=0.55, rely=0.1, relwidth=0.33, relheight=0.6)
		
		#Frame 6---------------------------------------------------------------------------------
		self.frame6 = tk.Frame(self.frameB, highlightbackground="black", highlightthickness=1, )
		self.frame6.place(relx=0.6, rely=0.05, relwidth=0.4, relheight=0.8)
		#	labels
		self.label_position_B = tk.Label(self.frame6, text = 'Position')
		self.label_position_B.grid(row=0, column=0, columnspan=1, sticky=tk.W)
		#	boxes
		self.input_rad_asc_B = tk.Entry(self.frame6, textvariable = self.rad_asc_str)
		self.input_rad_asc_B.config(width=30)
		self.input_rad_asc_B.grid(row = 2, column = 2, pady=5, padx=2)
		self.input_rad_dec_B = tk.Entry(self.frame6, textvariable = self.rad_dec_str)
		self.input_rad_dec_B.config(width=30)
		self.input_rad_dec_B.grid(row = 3, column = 2, pady=5, padx=2)
		
		#Frame 7---------------------------------------------------------------------------------
		self.frame7 = tk.Frame(self.frameC, highlightbackground="black", highlightthickness=1, )
		self.frame7.place(relx=0.6, rely=0.05, relwidth=0.4, relheight=0.8)
		#	labels
		self.label_position_B = tk.Label(self.frame7, text = 'Position')
		self.label_position_B.grid(row=0, column=0, columnspan=1, sticky=tk.W)
		#	boxes
		self.input_lat = tk.Entry(self.frame7, textvariable = self.latitude_str, width=10)
		self.input_lat.config(width=30)
		self.input_lat.grid(row = 2, column = 2, pady=5, padx=2)
		self.input_long = tk.Entry(self.frame7, textvariable = self.longitude_str, width=10)
		self.input_long.config(width=30)
		self.input_long.grid(row = 3, column = 2, pady=5, padx=2)
		
		#Frame 8---------------------------------------------------------------------------------
		self.frame8 = tk.Frame(self.frameA, highlightbackground="black", highlightthickness=1, )
		self.frame8.place(relx=0.6, rely=0.05, relwidth=0.4, relheight=0.8)
		#	labels
		self.label_observatory_A = tk.Label(self.frame8, text='Observation time')
		self.label_observatory_A.grid(row=0, column=0, pady=5, columnspan=2)
		
		self.label_sample_frec_A = tk.Label(self.frame8, text='Sample frec')
		self.label_sample_frec_A.grid(row=2, column=0, pady=5, columnspan=1, padx=5)
		self.label_sample_n_A = tk.Label(self.frame8, text='Sample N')
		self.label_sample_n_A.grid(row=2, column=2, pady=10, padx=5)
		self.label_sample_time_A = tk.Label(self.frame8, text='Sampling time')
		self.label_sample_time_A.grid(row=1, column=0, pady=5, columnspan=1, padx=5)
		self.label_hour_angle_A = tk.Label(self.frame8, text='Hour A')
		self.label_hour_angle_A.grid(row=1, column=2, pady=5, padx=10)
		#	spinboxes
		self.spin_s_time_A = tk.Entry(self.frame8, width=5)
		self.spin_s_time_A.config(textvariable = self.sample_time)
		self.spin_s_time_A.grid(row = 1, column = 1, pady=10, padx=5)
		
		self.spin_interval_A = tk.Entry(self.frame8, width=5)
		self.spin_interval_A.config(textvariable = self.sample_interval)
		self.spin_interval_A.grid(row = 2, column = 1, pady=5, padx=5)
		
		self.spin_s_number_A = tk.Entry(self.frame8, width=5)
		self.spin_s_number_A.config(textvariable = self.sample_number)
		self.spin_s_number_A.grid(row = 2, column = 3, pady=5, padx=5)
		
		self.spin_h_angle_A = tk.Entry(self.frame8, width=5)
		self.spin_h_angle_A.config(textvariable = self.hour_angle)
		self.spin_h_angle_A.grid(row = 1, column = 3, pady=5, padx=5)
		
		
		#	label
		self.label_algorythm_A = tk.Label(self.frame8, text='Algorythm')
		self.label_algorythm_A.grid(row = 3, column = 0, pady=10, padx=10)
		self.label_frecuency_A = tk.Label(self.frame8, text='Frecuency')
		self.label_frecuency_A.grid(row = 4, column = 0, pady=5, padx=10)
		#	checks fft & nufft
		self.check_fft_A = ttk.Checkbutton(self.frame8, text='FFT', variable=self.fft)
		self.check_fft_A.config(onvalue=True, offvalue=False)
		self.check_fft_A.config(command = self.unCheck)
		self.check_fft_A.grid(row = 3, column = 1, pady=10, padx=10)
		self.check_nufft_A = ttk.Checkbutton(self.frame8, text='NUFFT', variable=self.nufft)
		self.check_nufft_A.config(onvalue=True, offvalue=False)
		self.check_nufft_A.config(command = self.unCheck)
		self.check_nufft_A.grid(row = 3, column = 2, pady=10, padx=10)
		#	box to frecuency
		self.frecuency_list_A = tk.Entry(self.frame8, textvariable = 'enter frecuency')
		self.frecuency_list_A.grid(row=4, column=1, pady=5, padx=10, columnspan=3)
		
		
		# OUT of frames---------------------------------------------------------------------------
		self.run_button_out = tk.Button(self.window, text = 'Run', width=10)
		self.run_button_out.config(command = self.draw_graphic_f, state = tk.DISABLED)
		self.run_button_out.place(relx=0.85, rely=0)#grid(row = 0, column = 0)

