#!/usr/bin/env python
import cupy as cp
import numpy as np
from cupyx.scipy.signal import correlate
import sys
import argparse
import mrcfile
import pandas as pd
from matplotlib import pyplot as plt
from common.starparse import StarFile
from cupyx.scipy.ndimage import rotate, shift, zoom, gaussian_filter
import os
import csv
from concurrent.futures import ThreadPoolExecutor
import multiprocessing

def get_correlations(images, particles_df, output_file, start_from=0, batch_size=100):
    images = cp.asarray(images)
    num_particles = particles_df.shape[0]
    print("len of images", len(images))
    print("len of pardf", num_particles)
    def process_batch(start, end):
        correlations_batch = []
        for i in range(start, end):
            image_array = images[i]
            correlations = cp.zeros(num_particles, dtype=cp.float32)
            correlations = cp.sum(image_array * images, axis=(1,2))
            correlations[i] = 0  # Set self-correlation to 0
            correlations_batch.append(correlations.get())
        return correlations_batch

    with open(output_file, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        if start_from == 0:
            header = [''] + [str(i) for i in range(num_particles)]
            writer.writerow(header)
        
        for batch_start in range(start_from, num_particles, batch_size):
            batch_end = min(batch_start + batch_size, num_particles)
            correlations_batch = process_batch(batch_start, batch_end)
            
            for i, correlations in enumerate(correlations_batch, start=batch_start):
                row = [str(i)] + [f"{corr:.7f}" for corr in correlations]
				writer.writerow(row)
            csvfile.flush()
            
            print(f"Processed particles {batch_start} to {batch_end-1} of {num_particles}")
            
            with open('progress.txt', 'w') as f:
                f.write(str(batch_end))

def rotate_image(image, psi=0, x=0, y=0, order=1):
    image = cp.asarray(image)
    if x != 0 or y != 0:
        image = shift(image, (y,x), mode='wrap')
    if psi != 0:
        image = rotate(image, -psi, axes=(0,1), mode='constant', reshape=False, order=order)
    return image.get()

def readall(image):
    with mrcfile.mmap(image, mode='r') as mrc:
        images = mrc.data
        images = np.asarray(images)
    return images

def shift_image(image, particles_df, pixel_size):
    images_shifted = np.zeros_like(image)
    for line_number, particle in particles_df.iterrows():
        image_info = particle['_rlnImageName'].split('@')
        slice = int(image_info[0])-1
        x = float(particle['_rlnOriginXAngst']) / pixel_size
        y = float(particle['_rlnOriginYAngst']) / pixel_size
        images_shifted[slice] = rotate_image(image[slice], x=x, y=y, order=3)
    return images_shifted

def main():
    parser = argparse.ArgumentParser(description='Calculate correlations with real-time output and restart capability.')
    parser.add_argument('star_file', help='Input STAR file')
    parser.add_argument('image_file', help='Input image file')
    parser.add_argument('output_file', help='Output CSV file')
    parser.add_argument('--restart', action='store_true', help='Restart from last progress')
    parser.add_argument('--batch-size', type=int, default=100, help='Batch size for processing')
    args = parser.parse_args()

    star_file = StarFile(args.star_file)
    particles_df = star_file.particles_df
    optics_df = star_file.optics_df
    pixel_size = float(optics_df.loc[0,'_rlnImagePixelSize'])

    images = readall(args.image_file)
    images_shifted = shift_image(images, particles_df, pixel_size)
    print('shift finished')

    start_from = 0
    if args.restart and os.path.exists('progress.txt'):
        with open('progress.txt', 'r') as f:
            start_from = int(f.read().strip())
        print(f"Restarting from particle {start_from}")

    get_correlations(images_shifted, particles_df, args.output_file, start_from, args.batch_size)

if __name__ == "__main__":
    main()
