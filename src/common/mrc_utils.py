

from src.common.volume_utils import rotate_image
import cupy as cp
import mrcfile
from scipy import fftpack
import numpy as np


### Input image should be a two element list with slice number first, image file name second.
### eg. [1, 'image.mrcs']
def readslice( image ):
	with mrcfile.mmap(image[1], mode='r') as imagestack:
		zslice = int(image[0])-1
		image_array = imagestack.data[zslice]
		image_array = cp.asarray(image_array)
		
		if image_array.ndim < 2:
			#this happens when there is only one particle in a micrograph
			image_array = cp.asarray(imagestack.data)
	return image_array


def get_power_spectrum( image_array):
	boxsize = image_array.shape[0]
	# normalize input image
	image_array = (image_array - cp.mean(image_array)) / cp.std(image_array)
	# pad image	
	image_array = cp.pad(image_array, pad_width=boxsize, mode='constant', constant_values=0)
	
	fft = cp.fft.fft2(image_array)
	fft = cp.fft.fftshift(fft)
	power_spectrum = cp.abs(fft) ** 2
	center_y, center_x = power_spectrum.shape[0]//2, power_spectrum.shape[1]//2
	power_spectrum[center_y, center_x] = 0.0
	
	return power_spectrum

def get_average_power(particles_df):
	num_particles = particles_df.shape[0]
	
	for index, row in particles_df.iterrows():
		image = row['_rlnImageName'].split('@')
		psi = float(row['_rlnAnglePsi'])
		image_array = readslice(image)
		image_array = rotate_image(image_array, psi=psi, order=1)
		power_spectrum = get_power_spectrum(image_array)

		if index == 0:
			average_power = power_spectrum
		else:
			average_power += power_spectrum
		print(f"Processed {index+1}/{num_particles} particles", end='\r')
	average_power /= particles_df.shape[0]
	average_power = rotate_image(average_power, psi=-90, x=0, y=0, order=3)
	return average_power


# CTF correction by phase flipping
def correct_image(image, pixel_size, voltage, defocus_u, defocus_v, ast_angle, output_path=None):
	"""
	Apply CTF correction to a cryo-EM image.
	
	Parameters:
	-----------
	image_path : str
		Path to input MRC file
	pixel_size : float
		Pixel size in Angstroms
	voltage : float
		Microscope voltage in kV
	defocus_u : float
		Major axis defocus in Angstroms
	defocus_v : float
		Minor axis defocus in Angstroms
	ast_angle : float
		Astigmatism angle in radians
	output_path : str, optional
		Path to save corrected image
		
	Returns:
	--------
	ndarray
		Corrected image
	"""
	# Read the image
	image = image
	
	# Calculate FFT of the image
	fft = fftpack.fft2(image)
	
	# Calculate CTF
	ctf = create_astigmatic_ctf(image.shape,defocus_u, defocus_v, ast_angle,pixel_size, voltage)
	
	# Apply phase flipping correction
	corrected_fft = fft * np.sign(ctf)
	
	# Convert back to real space
	corrected_image = -cp.real(fftpack.ifft2(corrected_fft))
	
	return corrected_image

def create_astigmatic_ctf(size, defocus_u, defocus_v, ast_angle, pixel_size,voltage,Cs=2.7, amplitude_contrast=0.07):
	"""
	Calculate the Contrast Transfer Function (CTF).
	
	Parameters:
	-----------
	size : tuple
		Image size in pixels (height, width)
	pixel_size : float
		Pixel size in Angstroms
	voltage : float
		Microscope voltage in kV
	defocus_u : float
		Major axis defocus in Angstroms
	defocus_v : float
		Minor axis defocus in Angstroms
	ast_angle : float
		Astigmatism angle in radians
	Cs : float, optional
		Spherical aberration in mm
	amplitude_contrast : float, optional
		Amplitude contrast ratio
	
	Returns:
	--------
	ndarray
		2D array containing the CTF values
	"""
	# Convert units
	wavelength = 12.2639 / np.sqrt(voltage * 1000 + 0.97845 * voltage**2)  # electron wavelength
	Cs = Cs * 1e7  # convert mm to Angstroms
	
	# Create frequency arrays
	ny, nx = size
	y = cp.fft.fftfreq(ny, d=pixel_size)
	x = cp.fft.fftfreq(nx, d=pixel_size)
	X, Y = np.meshgrid(x, y)
	
	# Convert to polar coordinates
	k2 = X**2 + Y**2  # spatial frequency squared
	ang = cp.arctan2(Y, X) - ast_angle
	
	# Calculate defocus for each point
	defocus = 0.5 * (defocus_u + defocus_v + (defocus_u - defocus_v) * cp.cos(2 * ang))
	
	# Calculate phase shift
	chi = cp.pi * wavelength * k2 * (defocus - 0.5 * wavelength**2 * k2 * Cs)
	
	# Calculate CTF with amplitude contrast
	ctf = -cp.sqrt(1 - amplitude_contrast**2) * cp.sin(chi) - amplitude_contrast * cp.cos(chi)
	
	return ctf

