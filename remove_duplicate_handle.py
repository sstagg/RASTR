#! /usr/bin/env python

# handle correlation files
import cupy as cp
from cupyx.scipy.signal import correlate
from cupyx.scipy.ndimage import rotate, shift
import pandas as pd
from matplotlib import pyplot as plt
from starparse import StarFile
import mrcfile

# Function to do cross-correlation
def cross_correlation(a, b):
    # Compute the cross-correlation of two arrays
    # a and b are 1D arrays
    # Returns a 1D array of the cross-correlation of a and b
    return np.correlate(a, b, mode='full')

def readslice( image ):
	with mrcfile.mmap(image[1], mode='r') as imagestack:
		zslice = int(image[0])-1
		image_array = imagestack.data[zslice]
		image_array = cp.asarray(image_array)
	return image_array

### for rotating raw particle images
def rotate_image( image, psi=0, x=0, y=0, order=1 ):
	if x != 0 or y != 0:
		image = shift ( image, (y,x), mode='wrap' )
	if psi != 0:
		image = rotate ( image, -psi, axes=(0,1), mode='constant', reshape=False, order=order )
	return image


def main():
    # read pandas dataframe from csv file

    reader = pd.read_csv('refinej17bin2particle_all_correlations.csv', iterator=True, chunksize=1)

    #read star file
    star_file = StarFile('refinej17.star')
    particles_df = star_file.particles_df
    optics_df = star_file.optics_df
    pixel_size = float(optics_df.loc[0,'_rlnImagePixelSize'])

        
    #correlation_df = pd.read_csv('particle_all_correlations.csv', nrows=6) #nrows temporarily set to 5 for testing

    #index_temp = 5
    #load the first row into a cupy array
    #first_row = cp.array(correlation_df.iloc[index_temp])
    repeats = []
    for row in reader:
        first_row = cp.array(row)[0]
        image_number = int(first_row[0])+1
        if image_number in repeats:
            continue
        
    #print(rmsd)
        plt.plot(first_row.get())
        plt.show()

        mean_first_row = first_row.mean()
        rmsd = cp.sqrt(cp.mean(cp.square(first_row - mean_first_row)))
    #print index of first row when its value is larger than mean + 10*rmsd
        indexes = cp.where(first_row > mean_first_row + 10*rmsd)[0].get()
        #print(first_row.get())
        #print(indexes)
        
        #print(indexes)
        arrays = []
        for index in indexes:
            if index == 0:
                continue
        
            index = int(index)
		
            #print(index, particles_df.loc[index-1, '_rlnAnglePsi'], particles_df.loc[index-1, '_rlnOriginXAngst'], particles_df.loc[index-1, '_rlnOriginYAngst'])
            x = particles_df.loc[index-1, '_rlnOriginXAngst'] / pixel_size
            y = particles_df.loc[index-1, '_rlnOriginYAngst'] / pixel_size
            psi = particles_df.loc[index-1, '_rlnAnglePsi']
		
            image = readslice(particles_df.loc[index-1, '_rlnImageName'].split('@'))
            image = rotate_image(image, psi, x, y, 3)
            if index == int(first_row[0]) + 1:
                raw_image = image
            arrays.append(image)
        arrays = cp.array(arrays)	
        arrays = arrays * raw_image
        arrays = cp.sum(arrays, axis=(1,2))
        #print (arrays.get())
    ### print the first row's value at indexes
        #print(first_row[indexes].get())
        repeat = indexes[cp.where(arrays > mean_first_row + 10*rmsd)[0].get()]
        print(repeat)

        ### append index where its value in arrays is larger than mean + 10*rmsd
        for i in repeat:
            if int(i) == int(first_row[0]) + 1:
                #print('skip')
                continue
            if int(i) < image_number:
                continue
            if int(i) in repeats:
                continue
            repeats.append(int(i))
        ### print current rows' linenumber in dataframe
        print(image_number, end='\r')
    ### remove duplicate values in repeats

    
    pd.DataFrame(repeats).to_csv('repeatsj17.csv', index=False)
    print(len(repeats))
    
		



    

if __name__ == '__main__':
    main()
