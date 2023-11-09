#! /usr/bin/env python

### Script to export particles from cryosparc to relion
### The output is star file with correct imagename, you should do relion_stack_create on your own.
### Usage, cd to export path with file eg. J185_particles_selected_exported.cs first.
### Usage, then run ./csexport.py J185_particles_selected_exported.cs J185_particles_selected_exported.star


import os
import subprocess
import sys
from starparse import StarFile

def main():
	# Define temporary paths
	dir_name = "links"
	temp_star = "tempt.star"

	# get cryosparc project path. Current version 4.2.1
	# Please be aware the path pattern 'CS-' might be changed in future version
	current_path = os.getcwd()
	path_parts = current_path.split(os.sep)
	cs_project_path = os.sep.join((path_parts)[:next(i for i,part in enumerate(path_parts) if 'CS-' in part)+1])

	# Change .cs to .star
	csparcfile = sys.argv[1]
	outputfile = sys.argv[2]
	subprocess.run([f"~/pyem/csparc2star.py {csparcfile} {temp_star}"], shell=True)

	# Get particle stack file path
	starfile = StarFile ( temp_star )
	particles_df = starfile.particles_df


	# The output star file by csparc2star.py has image name: 
	#eg. 000001@>J172/imported/014723505993541221887_23jan20ajob006_x312r128.mrcs
	image = particles_df.loc[ 0, '_rlnImageName'].split('>')[-1]

	# Check if the particle files suffix
	# Suffix mrc cannot be directly used by relion. These are usually extractions by cryosparc itself
	if image.split('.')[-1] == 'mrc':

		image_dir = os.path.dirname( image )

		# Create mrcs file by softlinks
		# Create a directory named "links"
		os.makedirs(dir_name, exist_ok=True)

		# Create symbolic links for .mrc files
		subprocess.run(["ln -s "+cs_project_path+'/'+image_dir+"/*.mrc "+dir_name+"/."], shell=True)

		# Rename link files from .mrc to .mrcs
		subprocess.run(["rename.ul .mrc .mrcs "+dir_name+"/*.mrc"], shell=True)

		# Move .mrcs files to the particle directory
		subprocess.run(["mv "+dir_name+"/*.mrcs "+cs_project_path+'/'+image_dir+"/."], shell=True)

		# Remove the directory "links"
		os.rmdir(dir_name)

	# Execute Python script to change image path
	change_image_path(temp_star, outputfile, cs_project_path)

	# Remove temporary star file
	os.remove(temp_star)
	print ("Remember to run relion_stack_create xxx to get new stack file!")

### In 4.2.1 version, the output star file has _rlnImageName at first column, _rlnMicrographName at second column.
def change_image_path(file_1, file_2, job_prefix):

	starfile = StarFile( file_1 )
	particles_df = starfile.particles_df
	
	job_prefix = os.path.join(job_prefix, '')

	for line_number in range(particles_df.shape[0]):
		image = particles_df.loc[ line_number, '_rlnImageName']

		if image[-4:] == '.mrc':
			particles_df.loc[ line_number, '_rlnImageName'] = job_prefix.join(image.split('>')) + 's'
		elif image[-5:] == '.mrcs':
			particles_df.loc[ line_number, '_rlnImageName'] = job_prefix.join(image.split('>'))
		
		if '_rlnMicrographName' in particles_df.columns:
			micrograph = particles_df.loc[ line_number, '_rlnMicrographName']
			particles_df.loc[ line_number, '_rlnMicrographName'] = job_prefix.join(micrograph.split('>'))

	# Update the dataframe
	starfile.particles_df = particles_df

	# Write new file
	starfile.write( file_2 )

if __name__ == "__main__":
	main()
