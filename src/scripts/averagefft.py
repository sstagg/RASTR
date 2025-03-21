#! /usr/bin/env python3

from src.common.starparse import StarFile
from src.common.mrc_utils import  get_average_power
import sys
import mrcfile
import time


def main():
	# print total time taken
	start_time = time.time()
	starfile = StarFile(sys.argv[1])
	particle_df = starfile.particles_df
	optics_df = starfile.optics_df
	
	pixel_size = float(optics_df.loc[0,'_rlnImagePixelSize'])
	
	# Calculate average power spectrum
	average_power = get_average_power(particle_df)
	
	# Save outputs
	#mrcfile.write('average_power.mrc', average_power.get(), overwrite=True)
	with mrcfile.new('average_power.mrc', overwrite=True) as mrc:
		mrc.set_data(average_power.get().astype('float32'))
		mrc.voxel_size = pixel_size
	print("Total time taken: %s seconds" % (time.time() - start_time))

if __name__ == '__main__':
	main()
