#!/usr/bin/env python
### Script to export particles from cryosparc to relion
### The output is star file with correct imagename, you should do relion_stack_create on your own.
### Usage, cd to export path with file eg. J185_particles_selected_exported.cs first.
### Usage, then run ./csexport.py J185_particles_selected_exported.cs J185_particles_selected_exported.star

import os
import subprocess
import sys
from src.common.starparse import StarFile
import time

def main():
	# Define temporary paths
	temp_star = "csexport_tempt.star"
	
	# Change .cs to .star
	csparcfile = sys.argv[1]
	outputfile = sys.argv[2]
	
	subprocess.run([f"csparc2star.py {csparcfile} {temp_star}"], shell=True)

	# Execute Python script to change image path
	change_image_path(temp_star, outputfile)


	# Remove temporary star file
	os.remove(temp_star)

	print("Remember to run relion_stack_create xxx to get new stack file!")


### In 4.2.1 version, the output star file has *rlnImageName at first column, *rlnMicrographName at second column.
def change_image_path(file_1, file_2):
	starfile = StarFile(file_1)
	particles_df = starfile.particles_df
	current_path = os.path.join(os.getenv('PWD'), '')
	processed_path = set()  # Using a set instead of a list for faster lookups

	# Vectorized operations instead of row-by-row processing
	# Replace '>' at the beginning of image names
	if '_rlnImageName' in particles_df.columns:
		particles_df['_rlnImageName'] = particles_df['_rlnImageName'].str.replace('>', '')
		particles_df['_rlnImageName'] = particles_df['_rlnImageName'].str.replace('@', '@' + current_path)
		
		# Find all unique MRC files to process once
		mrc_mask = particles_df['_rlnImageName'].str.endswith('.mrc')
		mrc_files = particles_df.loc[mrc_mask, '_rlnImageName'].unique()
		
		# Process each unique MRC file only once
		for image in mrc_files:
			original_image = image.split('@')[1]
			image_path = original_image + 's'
			
			if image_path not in processed_path and not os.path.exists(image_path):
				os.symlink(original_image, image_path)
				processed_path.add(image_path)
		
		# Update all MRC paths at once
		particles_df.loc[mrc_mask, '_rlnImageName'] = particles_df.loc[mrc_mask, '_rlnImageName'] + 's'
	
	# Process micrograph paths
	if '_rlnMicrographName' in particles_df.columns:
		particles_df['_rlnMicrographName'] = particles_df['_rlnMicrographName'].str.replace('>', '')
		particles_df['_rlnMicrographName'] = current_path + particles_df['_rlnMicrographName']
	
	# Update the dataframe
	starfile.particles_df = particles_df
	# Write new file
	starfile.write(file_2)

if __name__ == "__main__":
	time_start = time.time()
	main()
	time_end = time.time()
	print(f"Time taken: {time_end - time_start:.2f} seconds")
