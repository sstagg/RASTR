#! /usr/bin/env python

import cupy as cp
import numpy as np
from common.starparse import StarFile
from cupyx.scipy.signal import correlate
import sys
from matplotlib import pyplot as plt
import mrcfile
import pandas as pd
# Function to get the cross-correlation profile of particles


def get_particle_correlations(particles_df):
    num_particles = particles_df.shape[0]
    # Initialize a DataFrame to store the maximum correlation and corresponding particle index
    particle_correlations = pd.DataFrame(columns=['MaxCorrelation', 'WithParticleIndex'], index=range(num_particles))

    particle_all_correlations = []
    print('Processing particles...')

    # read all slice into gpu memory
    slices = readall(particles_df)

    for i, particle_row in particles_df.iterrows():
        # print out the progress
        print(f'Processing particle {i} of {num_particles}...', end='\r')
        image_info = particle_row['_rlnImageName'].split('@')
        x = float( particle['_rlnOriginXAngst']) / pixel_size
        y = float( particle['_rlnOriginYAngst']) / pixel_size
        image_array_1 = slices[int(image_info[0])-1]
        max_correlation = -np.inf
        max_j = -1

        #a list to store the cross-correlation values
        cross_correlation_values = []
        #normalized_cross_correlation_values = []


        for j in range(num_particles):
            if i == j:  # Skip comparing the particle with itself
                pass
            other_particle_info = particles_df.iloc[j]['_rlnImageName'].split('@')
            image_array_2 = slices[int(other_particle_info[0])-1]
            correlation = correlate(image_array_1, image_array_2, mode='full', method='fft')
            #auto_correlation_1 = correlate(image_array_1, image_array_1, mode='full', method='fft')
            #auto_correlation_2 = correlate(image_array_2, image_array_2, mode='full', method='fft')
            #normalized_correlation = correlation / (cp.sqrt(auto_correlation_1 * auto_correlation_2))

            current_max = cp.asnumpy(correlation).max()
            
            if current_max > max_correlation:
                max_correlation = current_max
                max_j = j
            cross_correlation_values.append(current_max)
        
        #print(len(cross_correlation_values))
        particle_all_correlations.append(cross_correlation_values)
        #plt.plot(cross_correlation_values)
        #plt.show()


        particle_correlations.at[i, 'MaxCorrelation'] = max_correlation
        particle_correlations.at[i, 'WithParticleIndex'] = max_j + 1  # Adjusted for 1-based indexing if necessary

    #output particle_all_correlations to a csv file
    particle_all_correlations_df = pd.DataFrame(particle_all_correlations)
    particle_all_correlations_df.to_csv('refinej17bin2particle_all_correlations.csv')
    return particle_correlations

### Input image should be a two element list with slice number first, image file name second.
### eg. [1, 'image.mrcs']
def readslice( image ):
	with mrcfile.mmap(image[1], mode='r') as imagestack:
		zslice = int(image[0])-1
		image_array = imagestack.data[zslice]
		image_array = cp.asarray(image_array)
	return image_array

def readall( image ):
	with mrcfile.mmap('refinej17_refinej17bin2.mrcs', mode='r') as imagestack:
		image_array = imagestack.data
		image_array = cp.asarray(image_array)
	return image_array



def main():
    # Read in the star file
    star_file = StarFile(sys.argv[1])
    particles_df = star_file.particles_df

    # Get particles cross-correlation profile
    particle_correlations = get_particle_correlations(particles_df)

    # order correlations from highest to lowest
    particle_correlations = particle_correlations.sort_values('MaxCorrelation', ascending=False)

    plt.plot(particle_correlations)
    plt.show()

    # output particle_correlations to a csv file
    particle_correlations.to_csv('particle_correlations.csv')




if __name__ == '__main__':
    main()