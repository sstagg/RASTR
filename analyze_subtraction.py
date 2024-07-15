#! /usr/bin/env python
import mrcfile
import cupy as cp
from cupyx.scipy.ndimage import rotate, shift
import numpy as np
from starparse import StarFile
import pandas as pd
from matplotlib import pyplot as plt
import sys


# from raw to volume
def rotate_image( image, psi=0, x=0, y=0, order=1 ):
	if x != 0 or y != 0:
		image = shift ( image, (y,x), mode='wrap' )
	if psi != 0:
		image = rotate ( image, -psi, axes=(0,1), mode='constant', reshape=False, order=order )
	return image






def main():
    starfile = StarFile( sys.argv[1] )
    mrcsfile = sys.argv[2]
    particles_df = starfile.particles_df
    optics_df = starfile.optics_df
    pixel_size = float(optics_df.loc[0,'_rlnImagePixelSize'])
    with mrcfile.mmap(mrcsfile, mode='r') as imagestack:
	    slice_number = imagestack.data.shape[0]
	
    
    if slice_number != len(particles_df):
        print('The number of particles in star does not match the number of slices in the mrcs file')
        #sys.exit(1)

    rmsds = 0
    rmsd_list = []
    for line_number, particle in particles_df.iterrows():
        image = particle['_rlnImageName'].split('@')
        psi = float(particle['_rlnAnglePsi'])
        x = float(particle['_rlnOriginXAngst']) / pixel_size
        y = float(particle['_rlnOriginYAngst']) / pixel_size
        with mrcfile.mmap(mrcsfile, mode='r') as imagestack:
            image_array = imagestack.data[line_number]
        image_array = cp.asarray(image_array)
        image_array = rotate_image(image_array, psi, x, y)
        image_projection = cp.sum(image_array, axis=1)
        #plt.plot(image_projection.get())
        #plt.show()
        image_mean = cp.mean(image_projection)
        rmsd = cp.sqrt(cp.sum(cp.square(image_projection - image_mean)))
        rmsds += rmsd
        rmsd_list.append(rmsd.get())
    rmsds /= len(particles_df)
    print('The number of particles is: ', len(particles_df))
    print('The average rmsd is: ', rmsds)
    plt.hist(rmsd_list, bins=100)
    plt.show()
    


if __name__ == '__main__':
    main()