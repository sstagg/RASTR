import tkinter as tk
import tkinter.scrolledtext as st
from tkinter import ttk
from matplotlib import pyplot
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from src.common.mrc_utils import readslice, rotate_image
import cupy as cp
import numpy as np
from cupyx.scipy.signal import correlate
from cupyx.scipy.ndimage import gaussian_filter
from src.common.volume_utils import low_pass_filter, image_bin
from src.scripts.diameTR import find_peak

# only used in psi finding optimization for display purpose
def extract_sub_array(image, length, width, angle):
	center_y, center_x = cp.array(image.shape) / 2 
	angle_rad = cp.deg2rad(angle)
	y,x = cp.indices(image.shape)

	x_shifted = x - center_x
	y_shifted = y - center_y

	# width define the distance of pixel away from box axis
	width_rotated = x_shifted * cp.sin(angle_rad) + y_shifted * cp.cos(angle_rad)
	# length define the distance of the projection of pixel on axis to center of box, 
	length_rotated = -x_shifted * cp.cos(angle_rad) + y_shifted * cp.sin(angle_rad)

	# Define the mask for the sub-array by creating boolen array
	mask = (cp.abs(width_rotated) <= width // 2) & (cp.abs(length_rotated) <= length // 2)

	# Create an empty sub-array with the same shape as the image
	sub_array = cp.zeros_like(image)

	# Assign the values from the image to the sub-array using the mask
	sub_array[mask] = image[mask]

	return sub_array

### psi 40 is the same as psi 220 for this script's purpose, so we bring them all to 0-180
def angle_within180( angle ):
	return (angle + 360.0) % 180



def angle_within360( angle ):
	if angle < 0.0:
		angle += 360.0
		return angle_within360( angle )
	elif angle > 360.0:
		angle -= 360.0
		return angle_within360( angle )
	else:
		return angle

### GUI window to manually find the optimal parameters
### sigma and minimum gap
class optimize_diameter_parameter:
	def __init__(self, particles_df, pixel_size ):
		self.particles_df = particles_df
		self.particle_number = self.particles_df.shape[0]
		self.zslice = 1
		self.angle = 0
		self.sigma = 0
		self.min_gap = 0
		self.pixel_size = pixel_size
		
	def update_parameter(self):
		try:
			self.zslice = int(self.entry_slice.get())
		except: self.zslice = 1

		try:
			self.angle_str.set("{:4.2f}".format(self.particles_df.loc[ self.zslice-1, '_rlnAnglePsi']))
		except: self.angle_str.set(0.0)
		try:
			self.angle = float(self.entry_angle.get())
		except: 
			self.angle = self.particles_df.loc[ self.zslice-1, '_rlnAnglePsi']
		try:
			self.sigma = int(self.entry_sigma.get())
		except: self.sigma = 0
		try:
			self.lowpass = float(self.entry_lowpass.get())
		except: self.lowpass = 0.0
		try:
			self.min_gap = float(self.entry_min_gap.get())
		except: self.min_gap = 0.0

	def create_entry(self, parent, default_text=None):
		default_text = default_text or tk.StringVar()
		entry = ttk.Entry(parent, width=10, textvariable=default_text)
		entry.pack(side=tk.LEFT)
		entry.bind('<Return>', lambda event: self.update_plot())
		entry.bind('<KP_Enter>', lambda event: self.update_plot())
		return entry

	def update_plot(self):
		self.update_parameter()
		self.update_figure()


	def update_figure(self):
		zslice = self.zslice
		angle = self.angle
		sigma = self.sigma
		min_gap = self.min_gap
		lowpass = self.lowpass

		particle_line = self.particles_df.loc[zslice-1]
		image = particle_line['_rlnImageName'].split('@')
		image_array = readslice( image )
		image_array = gaussian_filter( image_array, sigma)
		if lowpass > 0:
			image_array = low_pass_filter( image_array, resolution=lowpass, pixel_size=self.pixel_size)
		image_array = rotate_image(image_array, psi=angle)
		image_array_vertical = rotate_image(image_array, psi=90.0)
		image_1d = cp.sum( image_array, axis=1)
		image_1d = image_1d - image_1d.min()
		peak_one, peak_two = find_peak(image_1d.get(), min_gap=min_gap)
	
		self.ax1.clear()
		self.ax1.imshow(image_array_vertical.get(), origin='lower', cmap='gray')
		self.ax2.clear()
		self.ax2.plot(image_1d.get())
		self.ax2.set_xlim(0,image_1d.shape[0]-1)
		self.ax2.set_box_aspect(0.5)
		self.ax2.scatter(peak_one, image_1d[peak_one].get(), color='red')
		self.ax2.scatter(peak_two, image_1d[peak_two].get(), color='red')

		self.add_log('image:%8i, peaks %s, %s, diameter: %.2f pixels' %(self.zslice, peak_one, peak_two, abs(peak_one - peak_two)))

		self.canvas.draw()

	def previous_slice(self):
		if self.zslice > 1:
			self.zslice -= 1
		self.slice_default.set(str(self.zslice))
		self.update_plot()

	def next_slice(self):
		if self.zslice < self.particle_number:
			self.zslice += 1
		self.slice_default.set(str(self.zslice))
		self.update_plot()
	def random_slice(self):
		self.zslice = np.random.randint(1, self.particle_number)
		self.slice_default.set(str(self.zslice))
		self.update_plot()

	def save_variable(self):
		self.update_parameter()
		self.optimized_sigma = self.sigma
		self.optimized_min_gap = self.min_gap
		self.optimized_lowpass = self.lowpass
		self.ax1.clear()
		self.ax2.clear()
		self.root.destroy()

	def add_log(self, log_text):
		self.log_area.configure(state='normal')
		self.log_area.insert(tk.END, log_text + '\n')
		self.log_area.configure(state='disabled')
		self.log_area.see(tk.END)


	def startwindow(self):
		### create tkinter window
		self.root = tk.Tk()
		self.root.title("diameter finding optimization")
	
		### Create the input box and button 
		frame_controls = ttk.Frame(self.root)
		frame_controls.pack(side=tk.TOP, padx=5, pady=5)
		self.frame_controls = frame_controls
		label_slice = ttk.Label(frame_controls, text="image:")
		label_slice.pack(side=tk.LEFT, padx=(0, 5))
		self.slice_default = tk.StringVar()
		self.slice_default.set('1')
		self.entry_slice = self.create_entry( self.frame_controls, self.slice_default)
		self.button_random = tk.Button(self.frame_controls, text="R", command=self.random_slice, padx=1, pady=1)
		self.button_random.pack(side=tk.LEFT, padx=(0,1))
		self.button_left = tk.Button(self.frame_controls, text="\u25C0", command=self.previous_slice, padx=1, pady=1)
		self.button_left.pack(side=tk.LEFT, padx=(0,1))
		self.button_right = tk.Button(self.frame_controls, text="\u25B6", command=self.next_slice, padx=1, pady=1)
		self.button_right.pack(side=tk.LEFT, padx=(0,1))



		label_angle = ttk.Label(frame_controls, text="angle:")
		label_angle.pack(side=tk.LEFT, padx=(10,5))
		self.angle_str = tk.StringVar()
		self.entry_angle = self.create_entry( self.frame_controls)
		angle_fromstar = ttk.Label(frame_controls, textvariable=self.angle_str)
		angle_fromstar.pack(side=tk.LEFT, padx=(10,5))

		label_sigma = ttk.Label(frame_controls, text='sigma:')
		label_sigma.pack(side=tk.LEFT, padx=(10,5))
		sigma_default = tk.StringVar()
		sigma_default.set('0')
		self.entry_sigma = self.create_entry( self.frame_controls )

		label_lowpass = ttk.Label(frame_controls, text='lowpass:')
		label_lowpass.pack(side=tk.LEFT, padx=(10,5))
		lowpass_default = tk.StringVar()
		lowpass_default.set('0')
		self.entry_lowpass = self.create_entry( self.frame_controls)

		label_min_gap = ttk.Label(frame_controls, text='gap:')
		label_min_gap.pack(side=tk.LEFT, padx=(10,5))
		self.entry_min_gap = self.create_entry( self.frame_controls )

		button_update = ttk.Button(frame_controls, text="Update", command=self.update_plot)
		button_update.pack(side=tk.LEFT, padx=(5, 0))

		button_save = ttk.Button(frame_controls, text="Save", command=self.save_variable)
		button_save.pack(side=tk.LEFT)

		### create log area
		self.log_area = st.ScrolledText( self.root, wrap=tk.WORD, width=60, height=10)
		self.log_area.pack(padx=10, pady=10)
		self.log_area.configure(state='disabled')

		### create Matplotlib plot
		image = self.particles_df.loc[ 0, '_rlnImageName'].split('@')
		image_array = readslice( image )
		image_1d = cp.sum(image_array, axis=1)
	
		self.fig, (self.ax1, self.ax2) = pyplot.subplots(2, 1, figsize=(4,6), gridspec_kw={'height_ratios':[2,1]}, layout="constrained")
		self.ax1.imshow(image_array.get(), origin='lower', cmap='gray')
		self.ax2.plot(image_1d.get())
		self.ax2.set_xlim(0, image_1d.shape[0]-1)
		self.ax2.set_box_aspect(0.5)
	
		### Display the plot in the tkinter window
		self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
		self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
		self.canvas.draw()
	
		### start Tkinter event loop
		self.root.mainloop()
	
	def report(self):
		return self.optimized_sigma, self.optimized_min_gap, self.optimized_lowpass, True

class optimize_angle_finding:
	def __init__(self, starfile):
		self.particles_df = starfile.particles_df
		self.optics_df = starfile.optics_df
		self.pixel_size = float(self.optics_df.loc[0,'_rlnImagePixelSize'])
		self.width = 2
		self.length = int(self.optics_df.loc[0,'_rlnImageSize'])
		self.fft_size = int(self.optics_df.loc[0,'_rlnImageSize'])
	def create_entry(self, parent, default_text=None):
		default_text = default_text or tk.StringVar()
		entry = ttk.Entry(parent, width=10, textvariable=default_text)
		entry.pack(side=tk.LEFT)
		entry.bind('<Return>', lambda event: self.update_plot())
		entry.bind('<KP_Enter>', lambda event: self.update_plot())
		return entry

	def get_fft(self, array, fft_size):
		correlation_array = array

		### the fft of auto correlation == square of fft of raw image. The speed to get the same final fft is similar.
		fft = cp.fft.fft2( correlation_array )
		fft_real = cp.abs( fft )
		real_shift = cp.fft.fftshift( fft_real )

		center_x, center_y = real_shift.shape[1]//2, real_shift.shape[0]//2
		real_shift[center_y, center_x] = 0.0
		half_box = fft_size / 2
		real_shift = real_shift[center_y-half_box:center_y+half_box , center_x-half_box:center_x+half_box]
		return real_shift

	def get_radon_array(self, array, angles, length=100, width=4 ):
		thetas = angles
		sums = cp.array([cp.sum(extract_sub_array(array, length, width, theta)) for theta in thetas])
		return [thetas.get(), sums.get()]

	# For making autocorrelation array higher contrast
	def normalize(self, array):
		array[array< array.mean()] = 0
		array = array - array.min()
		center_x, center_y = array.shape[1] //2, array.shape[0] //2
		array[center_y-5:center_y+5,center_x-5:center_x+5] = 0.0
		return array

	def update_array(self):
		image_array = self.image_array
		if self.sigma != 0:
			image_array = gaussian_filter( image_array, self.sigma)
			
		image_array = cp.pad(image_array, pad_width=self.pad, mode='constant', constant_values=0)
		
		if self.bin_factor != 1:
			image_array = image_bin( image_array, self.bin_factor)
		#print(image_array.shape)
			
		self.image_array = image_array
		self.autocorrelation_array = correlate ( image_array, image_array, mode='full', method='fft')
		center_x, center_y = self.autocorrelation_array.shape[1] //2, self.autocorrelation_array.shape[0] //2
		self.autocorrelation_array[center_y-1:center_y+1,center_x-1:center_x+1] = 0.0
		
		self.fft_array = self.get_fft( self.autocorrelation_array, self.fft_size)
		angles = cp.arange(0, 360, 2.0)
		self.autocorrelation_radon = self.get_radon_array(self.autocorrelation_array, angles, self.length, self.width)
		self.fft_radon = self.get_radon_array(self.fft_array, angles, self.length, self.width)
		
		if self.strategy_var.get() == 'fft':
			angle = self.fft_radon[0][np.argmax(self.fft_radon[1])] + 90
		elif self.strategy_var.get() == 'real':
			angle = self.autocorrelation_radon[0][np.argmax(self.autocorrelation_radon[1])]
		else: angle = 0
		
		angle = angle_within180(angle)
		self.rotated_image_array = rotate_image( self.image_array, psi=angle)
		
		
		

	def update_plot(self):
		self.update_parameter()
		self.update_figure()


	def update_figure(self):

		self.slice_default.set(str(self.slice))
		image = self.particles_df.loc[ self.slice-1, '_rlnImageName'].split('@')
		self.image_array = readslice( image )
		self.update_array()
		
		self.ax_image.imshow( self.rotated_image_array.get(), origin='lower', cmap='gray')
		self.ax1.imshow(self.normalize(self.autocorrelation_array.get()), origin='lower', cmap='gray')
		self.ax1.set_title('real')
		self.ax2.imshow(self.fft_array.get(), origin='lower', cmap='gray')
		self.ax2.set_title('fft')
		self.ax3.clear()
		self.ax3.plot( self.autocorrelation_radon[0], self.autocorrelation_radon[1] )
		self.ax3.set_xlabel(str(self.autocorrelation_radon[0][np.argmax(self.autocorrelation_radon[1])]))
		self.ax4.clear()
		self.ax4.plot( self.fft_radon[0], self.fft_radon[1] )
		self.ax4.set_xlabel(str(self.fft_radon[0][np.argmax(self.fft_radon[1])]))
		self.fig.tight_layout()
		self.canvas.draw()



	def update_parameter(self):
		try:
			self.slice = int(self.entry_slice.get())
		except: self.slice = 1

		try:
			self.fft_size = int(self.entry_fftsize.get())
		except: self.fft_size = 240

		try:
			self.width = int(self.entry_width.get())
		except: self.width = 2

		try:
			self.length = int(self.entry_length.get())
		except: self.length = 100
		
		try: self.pad = int(self.entry_pad.get())
		except: self.pad = 0
		
		try: self.sigma = int(self.entry_sigma.get())
		except: self.sigma = 0
		
		try: self.bin_factor = int(self.entry_bin.get())
		except: self.bin_factor = 1
		
	def save_variable(self):
		self.strategy = self.strategy_var.get()
		self.update_parameter()
		self.root.destroy()

	def report(self):
		return {
			'strategy': self.strategy,
			'width': self.width,
			'length': self.length,
			'fft_size': self.fft_size,
			'pad': self.pad,
			'sigma': self.sigma,
			'bin_factor': self.bin_factor,
		}
	
	def previous_slice(self):
		if self.slice > 1:
			self.slice -= 1
		self.update_figure()

	def next_slice(self):
		if self.slice < self.particles_df.shape[0]:
			self.slice += 1
		self.update_figure()
	
	def random_slice(self):
		self.slice = np.random.randint(1, self.particles_df.shape[0])
		self.update_figure()


	def start_window(self):
		### create tkinter window
		self.root = tk.Tk()
		self.root.title("psi finding optimization")

		### Create the input box and button
		self.frame_controls = ttk.Frame(self.root)
		self.frame_controls.pack(side=tk.TOP, padx=5, pady=5)

		label_slice = ttk.Label(self.frame_controls, text="image:")
		label_slice.pack(side=tk.LEFT, padx=(0, 5))
		self.slice_default = tk.StringVar()
		self.slice_default.set('1')
		self.entry_slice = self.create_entry( self.frame_controls, self.slice_default)
		self.button_random = tk.Button(self.frame_controls, text="R", command=self.random_slice, padx=1, pady=1)
		self.button_random.pack(side=tk.LEFT, padx=(0,1))
		self.button_left = tk.Button(self.frame_controls, text="\u25C0", command=self.previous_slice, padx=1, pady=1)
		self.button_left.pack(side=tk.LEFT, padx=(0,1))
		self.button_right = tk.Button(self.frame_controls, text="\u25B6", command=self.next_slice, padx=1, pady=1)
		self.button_right.pack(side=tk.LEFT, padx=(0,1))



		label_fftsize = ttk.Label(self.frame_controls, text="fft size:")
		label_fftsize.pack(side=tk.LEFT, padx=(10,5))
		fftsize_default = tk.StringVar()
		fftsize_default.set(str(self.fft_size*2))
		self.entry_fftsize = self.create_entry( self.frame_controls, fftsize_default )

		label_width = ttk.Label(self.frame_controls, text='width:')
		label_width.pack(side=tk.LEFT, padx=(10,5))
		width_default = tk.StringVar()
		width_default.set('2')
		self.entry_width = self.create_entry( self.frame_controls, width_default )

		label_length = ttk.Label(self.frame_controls, text='length:')
		label_length.pack(side=tk.LEFT, padx=(10,5))
		length_default = tk.StringVar()
		length_default.set(str(self.fft_size))
		self.entry_length = self.create_entry( self.frame_controls, length_default )

		label_pad = ttk.Label(self.frame_controls, text='pad:')
		label_pad.pack(side=tk.LEFT, padx=(10,5))
		default_pad = tk.StringVar()
		default_pad.set(str(self.fft_size//2))
		self.entry_pad = self.create_entry( self.frame_controls, default_pad)

		label_sigma = ttk.Label(self.frame_controls, text='sigma:')
		label_sigma.pack(side=tk.LEFT, padx=(10,5))
		default_sigma = tk.StringVar()
		default_sigma.set('0')
		self.entry_sigma = self.create_entry( self.frame_controls, default_sigma)
		
		label_bin = ttk.Label(self.frame_controls, text='bin')
		label_bin.pack(side=tk.LEFT, padx=(10,5))
		default_bin = tk.StringVar()
		default_bin.set('1')
		self.entry_bin = self.create_entry( self.frame_controls, default_bin)
		
		self.strategy_var = tk.StringVar()
		self.strategy_var.set('fft')
		ttk.Label(self.frame_controls, text='choose strategy:', justify=tk.LEFT).pack( side=tk.LEFT, padx=(10,10))
		self.R1 = ttk.Radiobutton(self.frame_controls, text='real', variable=self.strategy_var, value='real', command=self.update_figure)
		self.R1.pack( side=tk.LEFT, anchor = tk.W)
		self.R2 = ttk.Radiobutton(self.frame_controls, text='fft', variable=self.strategy_var, value='fft', command=self.update_figure)
		self.R2.pack( side=tk.LEFT, anchor = tk.W)

		button_update = ttk.Button(self.frame_controls, text="Update", command=self.update_plot)
		button_update.pack(side=tk.LEFT, padx=(10, 0))

		button_save = ttk.Button(self.frame_controls, text="Save", command=self.save_variable)
		button_save.pack(side=tk.LEFT)

		self.fig = pyplot.figure(figsize=(8,4))
		self.ax_image = pyplot.subplot(1,2,1)
		self.ax1 = pyplot.subplot(2,4,3)
		self.ax2 = pyplot.subplot(2,4,4)
		self.ax3 = pyplot.subplot(2,4,7)
		#self.ax3.set_yticks(())
		self.ax4 = pyplot.subplot(2,4,8)
		#self.ax4.set_yticks(())
		self.fig.tight_layout()
		#self.ax1, self.ax2, self.ax3, self.ax4 = axes.flatten()

		### Display the plot in the tkinter window
		self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
		self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
		self.canvas.draw()

		self.update_plot()

		### start Tkinter event loop
		self.root.mainloop()



class choose_threshold_window:
	def __init__(self, particles_df, label_size=16 ):
		self.particles_df = particles_df
		self.label_size = label_size
		self.startwindow()

	def startwindow(self):
		self.root = tk.Tk()
		self.root.title("choose diameter")

		### Create the input box and button
		frame_controls = ttk.Frame(self.root)
		frame_controls.pack(side=tk.TOP, padx=8, pady=8)


		label_min = ttk.Label(frame_controls, text="min:")
		label_min.pack(side=tk.LEFT, padx=(0, 5))
		self.entry_min = ttk.Entry(frame_controls, width=10)
		self.entry_min.pack(side=tk.LEFT)

		label_max = ttk.Label(frame_controls, text='max:')
		label_max.pack(side=tk.LEFT, padx=(10, 5))
		self.entry_max = ttk.Entry(frame_controls, width=10)
		self.entry_max.pack(side=tk.LEFT)

		# an option to whether or not showing the defocus, default off
		self.show_defocus = tk.BooleanVar()
		self.show_defocus.set(False)
		self.checkbutton = ttk.Checkbutton(frame_controls, text='show defocus', variable=self.show_defocus)
		self.checkbutton.pack(side=tk.LEFT, padx=(10, 0))

		self.button_save = ttk.Button(frame_controls, text="Save", command=self.save)
		self.button_save.pack(side=tk.LEFT, padx=(5, 0))


		self.fig, self.ax1 = pyplot.subplots(figsize=(6, 4))
		self.ax2 = self.ax1.twinx()
		self.fig.subplots_adjust(bottom=0.1)
		self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
		self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1, pady=5)

		# Bin width 2 angstrom
		min_val = self.particles_df['_rlndiameTR'].min()
		max_val = self.particles_df['_rlndiameTR'].max()
		bins = np.arange(min_val, max_val + 2, 2)

		self.ax1.hist(self.particles_df['_rlndiameTR'], bins=bins, color='blue', edgecolor='black')
		self.ax1.set_title('Diameter distribution')
		self.ax1.set_ylabel('#particles', fontsize=self.label_size)
		self.ax1.set_xlabel('diameter(A)', fontsize=self.label_size)
		self.ax1.grid(True)
		self.fig.tight_layout()


		# Calculate and plot mean defocus values for each bin
		bin_centers = (bins[:-1] + bins[1:]) / 2
		mean_defocus = []
		for i in range(len(bins)-1):
			mask = (self.particles_df['_rlndiameTR'] >= bins[i]) & \
				   (self.particles_df['_rlndiameTR'] < bins[i+1])
			mean_def = self.particles_df.loc[mask, '_rlnDefocusU'].mean()
			mean_defocus.append(mean_def if not np.isnan(mean_def) else 0)


 		# Plot defocus values on secondary y-axis
		self.ax2.plot(bin_centers, mean_defocus, color='red', marker='o', 
					 linestyle='-', label='Mean Defocus')
		self.ax2.set_ylabel('Defocus (Å)', fontsize=self.label_size, color='red')
		self.ax2.tick_params(axis='y', labelcolor='red')

		self.canvas.mpl_connect('button_press_event', self.on_mouse_move)

		self.coords_label = ttk.Label(self.root, text="")
		self.coords_label.pack()

		self.canvas.draw()
		self.root.mainloop()

	def update_plot(self):
		if self.show_defocus.get():
			self.ax2.set_visible(True)
		else:
			self.ax2.set_visible(False)
		
	def save(self):
		try:
			min_diameter = float(self.entry_min.get())
		except: min_diameter = 0.0
		try:
			max_diameter = float(self.entry_max.get())
		except: max_diameter = 10000.0
		self.min_diameter = min_diameter
		self.max_diameter = max_diameter
		self.thresholding()
		self.root.destroy()

	def thresholding(self):
		df = self.particles_df
		self.particles_df = df.loc[(df['_rlndiameTR'] >= self.min_diameter) & (df['_rlndiameTR'] <= self.max_diameter)]
	
	def on_mouse_move(self, event):
		if event.inaxes is not self.ax1 and event.inaxes is not self.ax2:
			return

		x, y = event.xdata, event.ydata
		if event.inaxes is self.ax1:
			self.coords_label['text'] = f'x={x:.2f}, y={y:.2f} (particles)'
		else:
			self.coords_label['text'] = f'x={x:.2f}, y={y:.2f} (defocus)'
