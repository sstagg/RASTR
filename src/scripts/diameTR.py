#! /usr/bin/env python

import logging
import cupy as cp
from cupyx.scipy.ndimage import gaussian_filter
from cupyx.scipy.signal import correlate
import numpy as np
import time
import argparse
import mrcfile
from scipy.signal import find_peaks
from matplotlib import pyplot
import math
from src.common.starparse import StarFile
from src.common.mrc_utils import readslice, get_average_power, correct_image
from src.common.volume_utils import rotate_image, low_pass_filter, image_bin, rotate_image_3dto2d
import pandas as pd


def find_psi( starfile, psi_optimizer, model_2d=None, mode='full' ):
	length = psi_optimizer['length']
	pad = psi_optimizer['pad']
	fft_size = psi_optimizer['fft_size']
	bin_factor = psi_optimizer['bin_factor']
	sigma = psi_optimizer['sigma']
	width = psi_optimizer['width']
	strategy = psi_optimizer['strategy']
	
	particles_df = starfile.particles_df
	particle_number = particles_df.shape[0]
	start_time = time.time()

	angle_maxs = cp.zeros(particle_number)
	logger.info("psi iteration 1 at:")
	for line_number, particle in particles_df.iterrows():
		
		print(f'{line_number + 1} /{particle_number}', end='\r')

		image = particle['_rlnImageName'].split('@')
		image_array = readslice(image)
		if sigma != 0:
			image_array = gaussian_filter(image_array, sigma=sigma)
		image_array = cp.pad(image_array, pad_width=pad, mode='constant')
		if bin_factor != 1:
			image_array = image_bin( image_array, bin_factor)

	

		if strategy == 'fft':
			angle_max = find_initial_psi_s1( image_array, angle_step=2.0, width=width, length=length, fft_size=fft_size, bin_factor=bin_factor )			
			angle_max = find_initial_psi_s1( image_array, center=angle_max-90, angle_range=10, angle_step=0.5, width=width, length=length, fft_size=fft_size, bin_factor=bin_factor)
		elif strategy == 'real':
			angle_max = find_initial_psi_s2( image_array, width=width, length=length )
			angle_max = find_initial_psi_s2( image_array, center=angle_max, angle_range=10, angle_step=0.5, width=width, length=length)
		angle_maxs[line_number] = angle_max
		#if abs( angle_max - angle_within180(particles_df.loc[line_number, '_rlnAnglePsi']) ) > 5:
			#pyplot.plot()
			#print(f"{line_number+1}, angle_max: {angle_max}, star: {particles_df.loc[line_number, '_rlnAnglePsi']}")
	print(f'{particle_number} /{particle_number}')
	logger.info("psi determination done        ")
	logger.info(f"time for psi finding: {time.time() - start_time}")

	return angle_maxs.get()



### radon of the center projection of auto correlation
def find_initial_psi_s2( image_array, center=90.0, angle_range=180.0, angle_step=5, width=4, length=200 ):

	correlation_array = correlate ( image_array, image_array, mode='full', method='fft')

	### auto correlation center peak always too high, removing to increase contrast
	center_x, center_y = correlation_array.shape[1] //2, correlation_array.shape[0] //2
	correlation_array[center_y-1:center_y+1,center_x-1:center_x+1] = 0.0

	### get a list of angles to test
	angle_min = center - 0.5 * angle_range
	angle_max = center + 0.5 * angle_range
	thetas = cp.arange(start=angle_min, stop=angle_max, step=angle_step)

	initial_angle = customized_radon( correlation_array, thetas, length, width)

	return angle_within180(initial_angle)


### radon of the center projection of fft 
def find_initial_psi_s1( image_array, center=90.0, angle_range=180, angle_step=5, width=4, length=100, fft_size=240, bin_factor=1 ):

	correlation_array = correlate ( image_array, image_array, mode='full', method='fft')
	center_x, center_y = correlation_array.shape[1] //2, correlation_array.shape[0] //2
	correlation_array[center_y-1:center_y+1,center_x-1:center_x+1] = 0.0

	### the fft of auto correlation == square of fft of raw image. The speed to get the same final fft is similar.
	fft = cp.fft.fft2( correlation_array )
	fft_real = cp.abs( fft )
	real_shift = cp.fft.fftshift( fft_real )

	### clip only the center part.
	center_x, center_y = real_shift.shape[1]//2, real_shift.shape[0]//2
	### somehow making center constant fft pixel zero can make it faster, and since center is used in all angles, changing it to any value won't affect psi determination
	real_shift[center_y, center_x] = 0.0
	half_box = fft_size / 2
	real_shift = real_shift[center_y-half_box : center_y+half_box , center_x-half_box : center_x+half_box]

	### get a list of angles to test
	angle_min = center - 0.5 * angle_range
	angle_max = center + 0.5 * angle_range
	thetas = cp.arange(start=angle_min, stop=angle_max, step=angle_step)

	initial_angle = customized_radon( real_shift, thetas, length, width)
	return angle_within180(initial_angle + 90.0)


### definition: theta is zero when parralel to x axis.
def customized_radon(image, thetas, length=100, width=4):
	thetas = cp.array(thetas)
	center_y, center_x = cp.array(image.shape) / 2 
	y,x = cp.indices(image.shape)

	x_shifted = x - center_x
	y_shifted = y - center_y

	angle_rads = cp.deg2rad(thetas)[:, None, None] # Add dimensions for broadcasting

	# width define the distance of pixel away from box axis
	width_rotated = x_shifted * cp.sin(angle_rads) + y_shifted * cp.cos(angle_rads)
	# length define the distance of the projection of pixel on axis to center of box, 
	length_rotated = -x_shifted * cp.cos(angle_rads) + y_shifted * cp.sin(angle_rads)

	# Define the mask for the sub-array by creating boolen array
	masks = (cp.abs(width_rotated) <= width // 2) & (cp.abs(length_rotated) <= length // 2)

	# Create an empty sub-array with the same shape as the image
	sub_arrays = cp.zeros_like(masks, dtype=image.dtype)

	# Assign the values from the image to the sub-array using the mask
	cp.copyto(sub_arrays, image, where=masks)

	sums = cp.sum(sub_arrays, axis=(1,2))

	return thetas[cp.argmax(sums)]




	

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
	
def find_peak( array_1d, min_gap=80):
	
	peak_one, peak_two = find_peak_strategy(array_1d, min_gap)
	
	return peak_one, peak_two



# maximum difference between positive and negative peak, best for protein decorated tubules
def find_peak_strategy( array_1d, min_gap=80):
	array_1d_np = cp.asnumpy(array_1d)
	center = len(array_1d_np) // 2
	
	array_1d_np = array_1d_np - array_1d_np.min()
	positive_peaks = find_peaks(array_1d_np)[0]
	
	array_1d_np_invert = array_1d_np.max() - array_1d_np
	negative_peaks = find_peaks(array_1d_np_invert)[0]
	
	diff_list = []
	for pos_peak in positive_peaks:
		if pos_peak > center:
			valid_neg_peaks = negative_peaks[negative_peaks > pos_peak]
		else:
			valid_neg_peaks = negative_peaks[negative_peaks < pos_peak]
		
		if len(valid_neg_peaks) > 0:
			closest_neg_peak = valid_neg_peaks[np.argmin(np.abs(valid_neg_peaks - pos_peak))]
			diff = np.abs(array_1d_np[pos_peak] - array_1d_np[closest_neg_peak])

			# find the positive peak
			#diff_list.append((diff, pos_peak))

			# find the negative particle
			#diff_list.append((diff, closest_neg_peak))				

			# find the middle of negative and positive peak
			diff_list.append((diff, 0.5*(pos_peak + closest_neg_peak)))

	diff_list.sort(reverse=True)
	left_peak, right_peak = None, None
	max_index = len(array_1d_np)-1
	for diff, index in diff_list[:6]:
		if left_peak is not None and right_peak is not None:
			break
		if abs(index - center) < min_gap//2:
			continue
		if index < center and index > 0 and left_peak is None:
			left_peak = index
		elif index > center and index < max_index and right_peak is None:
			right_peak = index
	if left_peak is None:
		left_peak = 0
	if right_peak is None:
		right_peak = 0
	return left_peak, right_peak
	#return sorted([index for diff, index in diff_list[:2]])
	

def find_diameter( particles_df, sigma, min_gap, resolution=0):
	
	diameters={'_rlndiameTR':[]}
	particle_number = particles_df.shape[0]
	print("finding diameter: ")
	for line_number in range(particle_number):
		print(f'{line_number + 1} /{particles_df.shape[0]}', end='\r')
		particle_line = particles_df.loc[line_number]

		image_array_rotated = prepare_find_diameter_image(particle_line, sigma, resolution)

		image_1d = cp.sum(image_array_rotated, axis=1)
		peak_one, peak_two = find_peak(image_1d, min_gap)
		

		diameter = abs((peak_one - peak_two) * pixel_size)
		
		diameters['_rlndiameTR'].append(diameter)
	print("done           ")
	return diameters

# sub function from find_diameter
def prepare_find_diameter_image(particle_line, sigma=0, resolution=0):
	psi = angle_within180( particle_line['_rlnAnglePsi'] )
	image = particle_line['_rlnImageName'].split('@')
	
	defocus_u = particle_line['_rlnDefocusU']
	defocus_v = particle_line['_rlnDefocusV']
	astig_angle = particle_line['_rlnDefocusAngle']
	
	image_array = readslice( image )
	#image_array = correct_image( image_array, pixel_size, 300.0, defocus_u, defocus_v, astig_angle)
	if sigma > 0:
		image_array = gaussian_filter(image_array , sigma)
	if resolution > 0:
		image_array = low_pass_filter(image_array, resolution, pixel_size)
	image_array_rotated = rotate_image( image_array, psi=psi)
	return image_array_rotated







def minimize_shift( particle, sigma, min_gap ):
	image = particle['_rlnImageName'].split('@')
	psi = particle['_rlnAnglePsi']

	image_array = readslice(image)
	image_array = gaussian_filter( image_array, sigma)
	image_array_rotated = rotate_image( image_array, psi=psi, x=0, y=0, order=3)
	image_projection = cp.sum(image_array_rotated, axis=1)
	peak_one, peak_two = find_peak( image_projection, min_gap )

	distance = (peak_one + peak_two)/2 - image_array_rotated.shape[1]//2
	x = -distance * math.sin( psi/180.0*math.pi )
	y = -distance * math.cos( psi/180.0*math.pi )

	return x*pixel_size , y*pixel_size


def get_average( optics_df, particles_df ):

	boxsize = int( optics_df.loc[ 0, '_rlnImageSize' ] )
	projection = cp.zeros(( boxsize, boxsize ))

	for line_number, particle in particles_df.iterrows():
		image = particle['_rlnImageName'].split('@')
		image_array = readslice(image)


		psi = float( particle['_rlnAnglePsi'])
		x = float( particle['_rlnOriginXAngst']) / pixel_size
		y = float( particle['_rlnOriginYAngst']) / pixel_size

		image_array_rotated = rotate_image( image_array, psi=psi, x=x, y=y, order=3)

		projection += image_array_rotated
		#pyplot.imshow(projection.get(), origin='lower', cmap='gray')
		#pyplot.show()
		print(f'Processed {line_number + 1} /{particles_df.shape[0]}', end='\r')
	projection /= particles_df.shape[0]

	return projection





def normalize_tube( model_2d ):
	### diameter is manually changed here
	diameter = 412
	radius = diameter//2
	center_y = model_2d.shape[0]//2
	boxsize = model_2d.shape[0]
	model_background = cp.copy( model_2d )
	model_background[center_y-radius:center_y+radius, :] = 0
	background_mean = cp.mean( model_2d )
	#pyplot.imshow(model_background.get(), origin='lower', cmap='gray')
	#pyplot.show()
	#pyplot.imshow(model_2d.get(), origin='lower', cmap='gray')
	#pyplot.show()
	model_2d = model_2d - background_mean
	#radius = 420
	model_2d[0:center_y-radius, :] = 0
	model_2d[center_y+radius:, :] = 0
	model_2d_mean = cp.mean( model_2d[:, int(0.1*boxsize):int(0.9*boxsize)], axis=1 )
	model_2d = cp.array( model_2d_mean*cp.ones((boxsize, boxsize)) )
	model_2d = cp.swapaxes(model_2d, 0, 1)
	#pyplot.imshow(model_2d.get(), origin='lower', cmap='gray')
	#pyplot.show()
	return model_2d


def tube_subtraction( starfile, model_2d, outputfilename, separate=False, diameter_step=10):
	current_starfile = starfile.copy()
	optics_df, particles_df = current_starfile.optics_df, current_starfile.particles_df
	if separate:
		# Keep only particles with valid diameters
		original_order = particles_df.index.tolist()
		
		particles_df = particles_df.loc[ particles_df['_rlndiameTR'] >= 0 ]
		particles_df = particles_df.loc[ particles_df['_rlndiameTR'] <= 450 ]
		
		# Compute the bins for historgram-like separation
		min_diameter = particles_df['_rlndiameTR'].min()
		max_diameter = particles_df['_rlndiameTR'].max()
		bins = list(range(int(min_diameter), int(max_diameter), diameter_step))

		# Split dataframe based on diamter and call the function recursively
		processed_subsets = []
		for i in range(len(bins)-1):
			lower_bound = bins[i]
			upper_bound = bins[i+1]
			starfile_subset = starfile.copy()
			particles_df_subset = particles_df.loc[ (particles_df['_rlndiameTR'] >= lower_bound) & (particles_df['_rlndiameTR'] < upper_bound) ]
			starfile_subset.particles_df = particles_df_subset
			starfile_subset.optics_df = optics_df
			if not particles_df_subset.empty:
				model_2d_subset = get_average( optics_df, particles_df_subset )
				processed_subset = tube_subtraction( starfile_subset, model_2d_subset, outputfilename.replace('.mrc', f'_{lower_bound}_{upper_bound}.mrc'))
				processed_subsets.append(processed_subset.particles_df)
				#print (len(processed_subsets))
				#print(processed_subset.particles_df.index.tolist())
		combined_df = pd.concat(processed_subsets, ignore_index=False)
		combined_df = combined_df.sort_index()
		current_starfile.particles_df = combined_df
	if not separate:
		
		if model_2d is None:
			model_2d = get_average( optics_df, particles_df )
		model_2d = normalize_tube( model_2d )
		#pyplot.imshow(model_2d.get(), origin='lower', cmap='gray')
		#pyplot.show()
		particle_number = particles_df.shape[0]

		scale_factors = []

		pool_size = 50
		scales_pool = cp.linspace(-1, 2, num=pool_size)
		scales_pools = scales_pool.reshape(pool_size,1,1)
	
		with mrcfile.new_mmap(outputfilename, shape=(particle_number, *model_2d.shape)
							  , mrc_mode=2, overwrite=True) as mrc:
			print(f'Writing {outputfilename} ...')
			image_number = 0
			for line_number, particle in particles_df.iterrows():
				psi = float( particle['_rlnAnglePsi'])
				x = float( particle['_rlnOriginXAngst']) / pixel_size
				y = float( particle['_rlnOriginYAngst']) / pixel_size
				model_2d_rotated = rotate_image_3dto2d( model_2d, psi=psi, x=x, y=y, order=3, rotate_mode='nearest')
				#pyplot.imshow(model_2d_rotated.get(), origin='lower', cmap='gray')
				#pyplot.show()
				model_2d_rotates = cp.array([model_2d_rotated] * pool_size)
				model_2d_scaled = model_2d_rotates * scales_pools
				image = particle['_rlnImageName'].split('@')
				image_array = readslice(image)

				image_arrays = cp.array([image_array] * pool_size)
				image_arrays = image_arrays - model_2d_scaled
				image_arrays = image_arrays * model_2d_rotates
				#for z in range(20):
				#	pyplot.imshow(image_arrays[z].get(), origin='lower', cmap='gray')
				#	pyplot.show()
				image_arrays = cp.abs(image_arrays)
				image_arrays = cp.sum(image_arrays, axis=(1,2))
				#image_arrays = cp.abs(image_arrays)
				#pyplot.plot(scales_pool.get(), image_arrays.get())
				#pyplot.show()
				scale_factor = scales_pool[cp.argmin(image_arrays)]
				scale_factors.append(scale_factor.item())
				#print(scale_factor.item())
				mrc.data[image_number] = ( image_array - scale_factor * model_2d_rotated ).get()
				image_number += 1
				#print(f'Processed {line_number + 1} /{particles_df.shape[0]}', end='\r')
				
				#update image name with new output filename
				particles_df.loc[line_number, '_rlnImageName'] =  str(image_number) + '@' + outputfilename
		#print(particles_df.shape)
		#print(scale_factors)
		#print(pd.DataFrame(scale_factors))
		particles_df['_rlnScaleFactorByRASTR'] = pd.Series(scale_factors).values
		#print(particles_df)
		#particles_df = pd.concat([particles_df, pd.DataFrame(scale_factors)], axis=1)

		#print(particles_df.shape)
		current_starfile.particles_df = particles_df
	return current_starfile

		

	

def parseOptions():
	parser = argparse.ArgumentParser()

	parser.add_argument('-i','--i', action='store', dest='input_star',
			help=' initial star file with correct psi angles')

	parser.add_argument('-o','--o', action='store', dest='output_rootname', default='RASTRtubeoutput',
			help=' output rootname')

	parser.add_argument('-m','--model', action='store', dest='model', default=None,
			help=' azimuthal average model as reference')

	parser.add_argument('-d','--diameter', action='store_true', dest='find_diameter', default=False,
			help=' option to classify images based on diameters')

	parser.add_argument('-p','--psi', action='store_true', dest='find_psi', default=False,
			help=' option to find the correct psi angle')

	parser.add_argument('-s','--shift', action='store_true', dest='minimize_shift', default=False,
			help=' option to minimize x, y shift based on correct psi angle')

	parser.add_argument('--particle_number', action='store', dest='particle_number', default=None, type=int,
			help=' option to process a small number of particles first')

	parser.add_argument('--classify', action='store_true', dest='justclassify', default=False, 
			help=' when you just want to classify based on diameters calculated already')

	parser.add_argument('--showaverage', action='store_true', dest='showaverage', default=False,
			help=' check diameter and shifts by calculating an average particle')
	
	parser.add_argument('--subtraction', action='store_true', dest='subtraction', default=False,
			help=' option to subtract the average particle from all particles')
	parser.add_argument('--average_power_spectrum', action='store_true', dest='average_power_spectrum', default=False,
			help=' option to calculate the average power spectrum')

	results=parser.parse_args()

	#make sure required arguments are present
	if results.input_star == None:
		pass

	return (results)


def setup_logger(output_rootname):
	"""Set up the logger for the script."""
	log_file = f"{output_rootname}_log.txt"
	
	# Create logger
	logger = logging.getLogger('RASTR')
	logger.setLevel(logging.INFO)
	
	# Create file handler
	file_handler = logging.FileHandler(log_file)
	file_handler.setLevel(logging.INFO)
	
	# Create console handler
	console_handler = logging.StreamHandler()
	console_handler.setLevel(logging.INFO)
	
	# Create formatter and add it to the handlers
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	file_handler.setFormatter(formatter)
	console_handler.setFormatter(formatter)
	
	# Add handlers to logger
	logger.addHandler(file_handler)
	logger.addHandler(console_handler)
	
	return logger

def log_parameters(logger, args_dict, stage="input", **kwargs):
	"""Log input arguments or optimized parameters."""
	if stage == "input":
		logger.info("Input arguments:")
		for arg, value in args_dict.items():
			logger.info(f"  {arg}: {value}")
	else:
		logger.info(f"Optimized parameters for {stage}:")
		for param, value in kwargs.items():
			logger.info(f"  {param}: {value}")


def main():
	# Parse Options
	results = parseOptions()

	# Set up logger
	global logger
	logger = setup_logger(results.output_rootname)

	# Log input arguments
	args_dict = vars(results)
	log_parameters(logger, args_dict, stage="input")

	# Start timing the execution
	start_time = time.time()
	logger.info("Starting RASTR...")


	starfile = StarFile( results.input_star )
	optics_df = starfile.optics_df

	if results.model is not None:
		logger.info(f"Loading model from {results.model}")
		model = mrcfile.read(results.model)
		model = cp.asarray(model)
		if model.ndim == 3:
			from src.common.volume_utils import project_volume
			model_3d = model
			model_2d = project_volume ( volume=model_3d, rot=0, tilt=90, psi=0, order=3 )	
			logger.info("Model is 3D, projecting to 2D")
		elif model.ndim == 2:
			model_2d = model
			logger.info("Loaded 2D model")
	else:
		model_2d = None
		logger.info("No model provided")


	if results.particle_number is not None:
		logger.info(f"Sampling {results.particle_number} particles from the star file")
		starfile.particles_df = starfile.particles_df.sample(n=int(results.particle_number)).reset_index(drop=True)
		
	particles_df = starfile.particles_df

	global pixel_size
	pixel_size = float(optics_df.loc[0,'_rlnImagePixelSize'])
	logger.info(f"Pixel size: {pixel_size} A/pixel")

	logger.info('Parsing finished, starting calculation')

	

	if results.justclassify:
		logger.info("Starting threshold window classification")
		from src.common.windows_utils import choose_threshold_window
		threshold_window = choose_threshold_window( starfile.particles_df )
		starfile.particles_df = threshold_window.particles_df
		logger.info("Classification finished")

	if results.find_psi:
		logger.info("Starting psi optimization")
		from src.common.windows_utils import optimize_angle_finding
		optimizer_psi = optimize_angle_finding( starfile )
		optimizer_psi.start_window()
		args = optimizer_psi.report()
		log_parameters(logger, {}, stage='PSI angle optimization', **args)
		
	if results.find_diameter:
		logger.info("Starting diameter optimization")
		from src.common.windows_utils import optimize_diameter_parameter
		optimized = False
		while not optimized:
			optimizer = optimize_diameter_parameter( particles_df, pixel_size )
			optimizer.startwindow()
			sigma, min_gap, lowpass, optimized = optimizer.report()
			log_parameters(logger, {}, stage='Diameter optimization', sigma=sigma, min_gap=min_gap, lowpass=lowpass)
			pyplot.close()
		logger.info("Diameter optimization finished")

	if results.find_psi:
		logger.info("Starting psi angle finding")
		
		updated_psi_values = find_psi( starfile, args, model_2d, mode='fast' )
	
		starfile.particles_df['_rlnAnglePsi'] = updated_psi_values
		particles_df = starfile.particles_df

		


	if results.find_diameter:
		logger.info("Starting diameter calculation")
		time_start = time.time()
		diameters = find_diameter( particles_df, sigma, min_gap, lowpass)
		logger.info(f'diameter calculation finished in {time.time() - time_start} seconds')
		diameters_df = pd.DataFrame(diameters)
		# If the diameter column already exists, replace it
		if '_rlndiameTR' in starfile.particles_df.columns:
			starfile.particles_df['_rlndiameTR'] = diameters_df['_rlndiameTR']
		else:
			starfile.particles_df = pd.concat([particles_df, diameters_df], axis=1)

		# Log diameter statistics
		if '_rlndiameTR' in starfile.particles_df.columns:
			diameter_stats = {
				'mean': starfile.particles_df['_rlndiameTR'].mean(),
				'std': starfile.particles_df['_rlndiameTR'].std(),
				'min': starfile.particles_df['_rlndiameTR'].min(),
				'max': starfile.particles_df['_rlndiameTR'].max()
			}
			log_parameters(logger, {}, stage="diameter statistics", **diameter_stats)


	if results.minimize_shift:
		logger.info("Starting shift minimization")
		shifts = [ minimize_shift( starfile.particles_df.iloc[line_number], sigma, min_gap) for line_number in range(starfile.particles_df.shape[0])]
		xshifts, yshifts = zip(*shifts)
		#print(xshifts)
		#print(yshifts)
		starfile.particles_df['_rlnOriginXAngst'] = xshifts
		starfile.particles_df['_rlnOriginYAngst'] = yshifts
		particles_df = starfile.particles_df
		logger.info("Shift minimization finished")


	if results.showaverage:
		logger.info("Calculating and saving average image")
		average = get_average( optics_df, particles_df )
		average = rotate_image( average, psi=-90, x=0, y=0, order=3)
		model_2d = average
		average_name = results.output_rootname + '_average.mrc'
		pyplot.clf()
		with mrcfile.new(average_name, overwrite=True) as f:
			f.set_data( average.get().astype('float32'))
			f.voxel_size = pixel_size
		logger.info(f"Average image saved to {average_name}")
		pyplot.imshow(average.get(), origin='lower', cmap='gray')
		pyplot.show()

		pyplot.plot(cp.sum(average, axis=0).get())
		pyplot.show()

	if results.subtraction:
		starfile = tube_subtraction( starfile, model_2d, results.output_rootname+'_subtracted.mrcs', separate=False, diameter_step=10)
		#starfile = tube_subtraction_projection( starfile, results.output_rootname+'_subtracted_projection.mrcs')		

	if results.average_power_spectrum:
		logger.info("Calculating average power spectrum")
		time_start = time.time()
		average_power = get_average_power( particles_df )
		time_end = time.time()
		output_file = results.output_rootname + '_average_power.mrc'
		pyplot.imshow(average_power.get(), origin='lower', cmap='gray')
		pyplot.show()
		with mrcfile.new(output_file, overwrite=True) as f:
			f.set_data( average_power.get().astype('float32'))
			f.voxel_size = pixel_size
		logger.info(f'Average power spectrum calculated in {time_end - time_start} seconds')
		logger.info(f"Average power spectrum saved to {output_file}")
	starfile.write( results.output_rootname + '.star')
	logger.info(f"Output star file saved to {results.output_rootname}.star")
	logger.info("diameTR finished")

if __name__ == '__main__':
	main()
