#! /usr/bin/env python
### originally developped to test tkinter GUI display, now function to show particle stacks.
### input can be either mrcs or star file.
### usage ./showstack.py  stack.mrcs
###       ./showstack.py  stack.star

import tkinter as tk
from tkinter import ttk
import numpy as np
from scipy.ndimage import rotate, shift
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import mrcfile
import sys
from scipy.ndimage import gaussian_filter
from starparse import StarFile




def rotate_image(image_array, psi=0, x=0, y=0, order=3 ):
	if x != 0 or y != 0:
		image_array = shift ( image_array, (y,x), mode='wrap' )
	if psi != 0:
		image_array = rotate ( image_array, -psi, axes=(0,1), mode='constant', reshape=False, order=order)
	return image_array



class ShowStack:
	def __init__(self, filename):
		self.filename = filename
		if self.filename.endswith('star'):
			self.starfile = StarFile( self.filename )
			self.optics_df = self.starfile.optics_df
			self.particles_df = self.starfile.particles_df
			self.particle_number = self.particles_df.shape[0]
			self.pixel_size = float(self.optics_df.loc[0,'_rlnImagePixelSize'])
		if self.filename.endswith('mrcs'):
			self.stackfile = filename
			self.particle_number = 10000
			self.pixel_size = 1.856
		self.zslice = 1
		self.angle = 0.0
		self.sigma = 0

		self.start_window()

	def readstack(self, filename, zslice=1):

		if filename.endswith('mrcs'):
			with mrcfile.mmap(filename, mode='r') as imagestack:
				#this is more than one particle stacks
				#image_array = imagestack.data[zslice-1]

				#hacking mode
				image_array = imagestack.data
			return image_array

		elif filename.endswith('star'):
			image = self.particles_df.loc[ zslice-1, '_rlnImageName'].split('@')
			with mrcfile.mmap(image[1], mode='r') as imagestack:
				image_slice = int(image[0])-1
				image_array = imagestack.data[image_slice]
			x_shift = float( self.particles_df.loc[ self.zslice-1, '_rlnOriginXAngst']) / self.pixel_size
			y_shift = float( self.particles_df.loc[ self.zslice-1, '_rlnOriginYAngst']) / self.pixel_size
			center_y, center_x = np.array(image_array.shape) // 2
			self.x, self.y = center_x - x_shift, center_y - y_shift
			self.psi = self.particles_df.loc[ self.zslice-1, '_rlnAnglePsi']
			self.angle_str.set("{:.2f}".format(self.psi))
			return image_array
		
		
	def start_window(self):
		self.root = tk.Tk()
		self.root.title("Interactive Plot")
		self.fig, self.ax = plt.subplots(figsize=(5, 3))
		self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
		self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
		self.canvas.draw()

		self.frame_controls = ttk.Frame(self.root)
		self.frame_controls.pack(side=tk.TOP, padx=5, pady=5)

		self.label_slice = ttk.Label(self.frame_controls, text="slice:")
		self.label_slice.pack(side=tk.LEFT, padx=(0, 5))
		self.zslice_text = tk.StringVar()
		self.zslice_text.set(str(self.zslice))
		self.entry_slice = self.create_entry( self.frame_controls, self.zslice_text)

		self.button_random = tk.Button(self.frame_controls, text='R', command=self.random_slice, padx=1, pady=1)
		self.button_random.pack(side=tk.LEFT, padx=(0,1))
		self.button_left = tk.Button(self.frame_controls, text="pre", command=self.previous_slice, padx=1, pady=1)
		self.button_left.pack(side=tk.LEFT, padx=(0,1))
		self.button_right = tk.Button(self.frame_controls, text="next", command=self.next_slice, padx=1, pady=1)
		self.button_right.pack(side=tk.LEFT, padx=(0,1))


		self.label_angle = ttk.Label(self.frame_controls, text="angle:")
		self.label_angle.pack(side=tk.LEFT, padx=(10,5))
		self.angle_text = tk.StringVar()
		self.entry_angle = self.create_entry( self.frame_controls, self.angle_text)
		self.angle_str = tk.StringVar()
		angle_fromstar = ttk.Label(self.frame_controls, textvariable=self.angle_str)
		angle_fromstar.pack(side=tk.LEFT, padx=(10,5))

		self.label_sigma = ttk.Label(self.frame_controls, text='sigma:')
		self.label_sigma.pack(side=tk.LEFT, padx=(10,5))
		self.entry_sigma = self.create_entry( self.frame_controls)

		self.show_center_yes = tk.BooleanVar()
		self.show_center_yes.set(False)
		self.C_center = ttk.Checkbutton(self.frame_controls, text='center', command=self.update_plot, variable=self.show_center_yes)
		self.C_center.pack(side=tk.LEFT, padx=1)		
		self.rotate_yes = tk.BooleanVar()
		self.rotate_yes.set(False)
		self.C_rotate = ttk.Checkbutton(self.frame_controls, text='rotate', command=self.update_plot, variable=self.rotate_yes)
		self.C_rotate.pack(side=tk.LEFT, padx=1)

		self.button_update = ttk.Button(self.frame_controls, text="Update", command=self.update_variable)
		self.button_update.pack(side=tk.LEFT, padx=(5, 0))

		self.root.mainloop()

	def next_slice(self):
		self.zslice += 1
		self.update_plot()
	
	def previous_slice(self):
		if self.zslice > 1:
			self.zslice -= 1
		else:
			print(' already the first slice ')
		self.update_plot()

	def random_slice(self):
		# the max is under modify.
		self.zslice = np.random.randint(1, self.particle_number)
		self.update_plot()

	def update_variable(self):
		try:
			zslice = int(self.entry_slice.get())
		except:
			zslice = 1
		try:
			angle = float(self.entry_angle.get())
		except:
			angle = 0.0
		try:
			sigma = int(self.entry_sigma.get())
		except:
			sigma = 0

		self.zslice = zslice
		self.angle = angle
		self.sigma = sigma
		self.update_plot()

	def update_plot(self):
		self.zslice_text.set(str(self.zslice))
		try:
			self.ax.clear()
			image_array = self.readstack(self.filename, self.zslice)
			image_array = gaussian_filter(image_array, self.sigma)
			if self.rotate_yes.get():
				image_array = rotate_image(image_array, psi=self.psi)
			else:
				image_array = rotate_image(image_array, psi=self.angle)
			self.ax.imshow(image_array, origin='lower', cmap='gray')
			if self.show_center_yes.get():
				self.ax.scatter(self.x, self.y, color='red')
			self.canvas.draw()
		except ValueError:
			pass
	

	def create_entry(self, parent, default_text=None):
		default_text = default_text or tk.StringVar()
		entry = ttk.Entry(parent, width=10, justify=tk.CENTER, textvariable=default_text)
		entry.pack(side=tk.LEFT)
		entry.bind('<Return>', lambda event: self.update_variable())
		entry.bind('<KP_Enter>', lambda event: self.update_variable())
		return entry

if __name__ == '__main__':
	ShowStack(sys.argv[1])




