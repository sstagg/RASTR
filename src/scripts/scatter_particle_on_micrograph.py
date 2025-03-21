#! /usr/bin/env python

from common.starparse import StarFile
import sys
import mrcfile

import cupy as cp
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from cupyx.scipy.ndimage import gaussian_filter
def extract_micrograph_path(full_path):
    """Extract just the filename from a full path."""
    return os.path.basename(full_path)

def readslice(image):
	with mrcfile.mmap(image, mode='r') as imagestack:
		image_array = imagestack.data
		image_array = cp.asarray(image_array)
		
		if image_array.ndim < 2:
			image_array = cp.asarray(imagestack.data)
	image_array = gaussian_filter( image_array, sigma=2)
	image_array = cp.flip( image_array, axis=0)
	return image_array

def plot_coordinates(df1, df2, output_dir="micrograph_plots"):
    """
    Load two RELION .star files, and for each micrograph, plot the coordinates
    from each file in different colors.
    """
    
    
    # Convert columns to appropriate types
    for df in [df1, df2]:
        if '_rlnCoordinateX' in df.columns:
            df['_rlnCoordinateX'] = df['_rlnCoordinateX'].astype(float)
        if 'CoordinateY' in df.columns:
            df['_rlnCoordinateY'] = df['_rlnCoordinateY'].astype(float)
    
    # Get unique micrographs
    micrographs1 = df1['_rlnMicrographName'].unique()
    micrographs2 = df2['_rlnMicrographName'].unique()
    all_micrographs = set(micrographs1) | set(micrographs2)
    
    print(f"Found {len(all_micrographs)} unique micrographs")
    
    # Process each micrograph
    for micrograph_path in all_micrographs:
        micrograph_name = extract_micrograph_path(micrograph_path)
        print(f"Processing {micrograph_name}...")
        
        # Get coordinates for this micrograph
        coords1 = df1[df1['_rlnMicrographName'] == micrograph_path][['_rlnCoordinateX', '_rlnCoordinateY']]
        coords2 = df2[df2['_rlnMicrographName'] == micrograph_path][['_rlnCoordinateX', '_rlnCoordinateY']]
        
        # Load micrograph
        img = readslice(micrograph_path)
        
        # Plot
        plt.figure(figsize=(10, 10))
        
        # Display micrograph in grayscale
        plt.imshow(img.get(), origin='lower', cmap='gray')
        
        # Plot coordinates
        if not coords1.empty:
            plt.scatter(coords1['_rlnCoordinateX'], coords1['_rlnCoordinateY'], 
                       color='red', s=30, alpha=0.7, label=f'File 1 ({len(coords1)} particles)')
        
        if not coords2.empty:
            plt.scatter(coords2['_rlnCoordinateX'], coords2['_rlnCoordinateY'], 
                       color='blue', s=30, alpha=0.7, label=f'File 2 ({len(coords2)} particles)')
        
        plt.title(f"Micrograph: {micrograph_name}")
        plt.legend()
        plt.axis('off')
        plt.show()        



def main():
	starfile_1 = StarFile( sys.argv[1])
	starfile_2 = StarFile( sys.argv[2])
	particles_df_1 = starfile_1.particles_df
	particles_df_2 = starfile_2.particles_df

	plot_coordinates( particles_df_1, particles_df_2 )	

if __name__ == '__main__':
	main()
