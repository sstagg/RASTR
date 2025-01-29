#! /usr/bin/env python3

import cupy as cp
from cupyx.scipy.ndimage import rotate, shift
from cupyx.scipy.signal import correlate
from matplotlib import pyplot as plt
import numpy as np
from starparse import StarFile
import sys
import mrcfile
import time

def get_fft_real(image, psi=0.0):
	image_array = readslice(image)
	boxsize = image_array.shape[0]
	
	# Normalize input image
	image_array = (image_array - cp.mean(image_array)) / cp.std(image_array)
	
	# pad image to 2^n with zeros
	image_array = cp.pad(image_array, pad_width=boxsize, mode='constant')
	image_array = rotate_image(image_array, psi=psi)

	# Calculate FFT
	fft = cp.fft.fft2(image_array)
	fft = cp.fft.fftshift(fft)
	power_spectrum = cp.abs(fft)**2  # Calculate power spectrum
	#power_spectrum = fft

	
	# Zero out central pixel to remove DC component
	center_y, center_x = fft.shape[0]//2, fft.shape[1]//2
	power_spectrum[center_y, center_x] = 0.0



	return power_spectrum


def rotate_image(image, psi=0, x=0, y=0, order=1):
	if x != 0 or y != 0:
		image = shift(image, (y,x), mode='wrap')
	if psi != 0:
		image = rotate(image, -psi, axes=(0,1), mode='constant', reshape=False, order=order)
	return image

def readslice(image):
	with mrcfile.mmap(image[1], mode='r') as imagestack:
		zslice = int(image[0])-1
		image_array = imagestack.data[zslice]
		image_array = cp.asarray(image_array)
		
		if image_array.ndim < 2:
			image_array = cp.asarray(imagestack.data)
	return image_array

def average_fft(particle_df):
	for index, row in particle_df.iterrows():
		image = row['_rlnImageName'].split('@')
		psi = float(row['_rlnAnglePsi'])
		power_spectrum = get_fft_real(image, psi)
		
		if index == 0:
			average_power = power_spectrum
		else:
			average_power += power_spectrum
			
	average_power = average_power / len(particle_df)
	average_power = rotate_image(average_power, psi=90)
	return average_power


def main():
	# print total time taken
	start_time = time.time()
	starfile = StarFile(sys.argv[1])
	particle_df = starfile.particles_df
	optics_df = starfile.optics_df
	global pixel_size
	pixel_size = float(optics_df.loc[0,'_rlnImagePixelSize'])
	
	# Calculate average power spectrum
	average_power = average_fft(particle_df)
	
	# Save outputs
	#mrcfile.write('average_power.mrc', average_power.get(), overwrite=True)
	with mrcfile.new('average_power.mrc', overwrite=True) as mrc:
		mrc.set_data(average_power.get().astype('float32'))
		mrc.voxel_size = pixel_size
	print("Total time taken: %s seconds" % (time.time() - start_time))

if __name__ == '__main__':
	main()
