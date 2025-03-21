#! /usr/bin/env python

'''
test for ctf correction
'''

import mrcfile
import os
import cupy as cp
from matplotlib import pyplot
from common.starparse import StarFile
from scipy.ndimage import gaussian_filter,rotate, shift
import numpy as np
from scipy import fftpack
import mrcfile
import matplotlib.pyplot as plt

def calculate_ctf(size, pixel_size, voltage, defocus_u, defocus_v, ast_angle, Cs=2.7, amplitude_contrast=0.07):
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
	y = fftpack.fftfreq(ny, d=pixel_size)
	x = fftpack.fftfreq(nx, d=pixel_size)
	X, Y = np.meshgrid(x, y)
	
	# Convert to polar coordinates
	k2 = X**2 + Y**2  # spatial frequency squared
	ang = np.arctan2(Y, X) - ast_angle
	
	# Calculate defocus for each point
	defocus = 0.5 * (defocus_u + defocus_v + (defocus_u - defocus_v) * np.cos(2 * ang))
	
	# Calculate phase shift
	chi = np.pi * wavelength * k2 * (defocus - 0.5 * wavelength**2 * k2 * Cs)
	
	# Calculate CTF with amplitude contrast
	ctf = -np.sqrt(1 - amplitude_contrast**2) * np.sin(chi) - amplitude_contrast * np.cos(chi)
	
	return ctf

def readslice( image ):
	with mrcfile.mmap(image[1], mode='r') as imagestack:
		zslice = int(image[0])-1
		image_array = imagestack.data[zslice]
		image_array = np.asarray(image_array)
		
		if image_array.ndim < 2:
			#this happens when there is only one particle in a micrograph
			image_array = np.asarray(imagestack.data)
	return image_array


def correct_image(image_path, pixel_size, voltage, defocus_u, defocus_v, ast_angle, output_path=None):
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
	image = readslice(image_path)
	
	# Calculate FFT of the image
	fft = fftpack.fft2(image)
	
	# Calculate CTF
	ctf = calculate_ctf(image.shape, pixel_size, voltage, defocus_u, defocus_v, ast_angle)
	
	# Apply phase flipping correction
	corrected_fft = fft * np.sign(ctf)
	
	# Convert back to real space
	corrected_image = np.real(fftpack.ifft2(corrected_fft))
	
	# Save the corrected image if output path is provided
	if output_path:
		with mrcfile.new(output_path, overwrite=True) as mrc:
			mrc.set_data(corrected_image.astype(np.float32))
	
	return corrected_image

def rotate_image( image, psi=0, x=0, y=0, order=1 ):
	if x != 0 or y != 0:
		image = shift ( image, (y,x), mode='wrap' )
	if psi != 0:
		image = rotate ( image, -psi, axes=(0,1), mode='constant', reshape=False, order=order )
	return image

def visualize_correction(original_image, corrected_image, ctf=None, psi=0):
	"""
	Visualize the original and corrected images, and optionally the CTF.
	
	Parameters:
	-----------
	original_image : ndarray
		Original image
	corrected_image : ndarray
		CTF-corrected image
	ctf : ndarray, optional
		Contrast Transfer Function
	"""
	if ctf is None:
		fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
	else:
		fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
	
	# Plot original image
	original_fft = fftpack.fftshift(fftpack.fft2(original_image))
	original_fft = np.abs(original_fft)
	original_image = gaussian_filter(original_image, 4)
	original_image = rotate_image( original_image, psi )
	image_1d_1 = np.sum(original_image, axis=1)
	ax1.plot(image_1d_1)
	ax1.set_title('Original Image')
	ax1.axis('off')
	
	# Plot corrected image
	corrected_fft = fftpack.fftshift(fftpack.fft2(corrected_image))
	corrected_fft = np.abs(corrected_fft)
	corrected_image = gaussian_filter(corrected_image, 4)
	corrected_image = rotate_image( corrected_image, psi )
	image_1d = -np.sum(corrected_image, axis=1)
	ax2.plot(image_1d)
	ax2.set_title('Corrected Image')
	ax2.axis('off')
	
	# Plot CTF if provided
	if ctf is not None:
		ax3.imshow(np.sign(ctf), cmap='gray')
		ax3.set_title('CTF')
		ax3.axis('off')
	
	plt.tight_layout()
	plt.show()

# Example usage
if __name__ == "__main__":
	# Example parameters
	
	pixel_size = 1.548  # Angstroms
	voltage = 300  # kV

	starfile = StarFile("J83_particles_exported.star")
	particles_df = starfile.particles_df
	slice = 85
	psi = float(particles_df["_rlnAnglePsi"][slice])
	image_path = particles_df["_rlnImageName"][slice].split("@")
	defocus_u = particles_df["_rlnDefocusU"][slice]
	defocus_v = particles_df["_rlnDefocusV"][slice]
	ast_angle = particles_df["_rlnDefocusAngle"][slice]

	
	# Correct the image
	corrected_image = correct_image(
		image_path,
		pixel_size,
		voltage,
		defocus_u,
		defocus_v,
		ast_angle,
		output_path="corrected_image.mrc"
	)
	
	# Calculate CTF for visualization
	ctf = calculate_ctf(
		corrected_image.shape,
		pixel_size,
		voltage,
		defocus_u,
		defocus_v,
		ast_angle
	)
	
	# Visualize results
	original_image = readslice(image_path)
	
	visualize_correction(original_image, corrected_image, ctf, psi)