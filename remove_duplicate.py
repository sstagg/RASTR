#! /usr/bin/env python

import cupy as cp
import numpy as np
from starparse import StarFile
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
    for i, particle_row in particles_df.iterrows():
        # print out the progress
        print(f'Processing particle {i} of {num_particles}...', end='\r')
        image_info = particle_row['_rlnImageName'].split('@')
        image_array_1 = readslice(image_info)
        max_correlation = -np.inf
        max_j = -1

        #a list to store the cross-correlation values
        cross_correlation_values = []



        for j in range(num_particles):
            if i == j:  # Skip comparing the particle with itself
                pass
            other_particle_info = particles_df.iloc[j]['_rlnImageName'].split('@')
            image_array_2 = readslice(other_particle_info)
            correlation = correlate(image_array_1, image_array_2, mode='full', method='fft')
            current_max = cp.asnumpy(correlation).max()
            
            if current_max > max_correlation:
                max_correlation = current_max
                max_j = j
            cross_correlation_values.append(current_max)
        
        particle_all_correlations.append(cross_correlation_values)
        #plt.plot(cross_correlation_values)
        #plt.show()


        particle_correlations.at[i, 'MaxCorrelation'] = max_correlation
        particle_correlations.at[i, 'WithParticleIndex'] = max_j + 1  # Adjusted for 1-based indexing if necessary

    #output particle_all_correlations to a csv file
    particle_all_correlations_df = pd.DataFrame(particle_all_correlations)
    particle_all_correlations_df.to_csv('particle_all_correlations.csv')
    return particle_correlations

### Input image should be a two element list with slice number first, image file name second.
### eg. [1, 'image.mrcs']
def readslice( image ):
	with mrcfile.mmap(image[1], mode='r') as imagestack:
		zslice = int(image[0])-1
		image_array = imagestack.data[zslice]
		image_array = cp.asarray(image_array)
	return image_array


# Function to do cross-correlation
def cross_correlation(a, b):
    # Compute the cross-correlation of two arrays
    # a and b are 1D arrays
    # Returns a 1D array of the cross-correlation of a and b
    return np.correlate(a, b, mode='full')



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