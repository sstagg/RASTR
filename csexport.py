#!/usr/bin/env python
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
    
    
    # Change .cs to .star
    csparcfile = sys.argv[1]
    outputfile = sys.argv[2]
    
    subprocess.run([f"csparc2star.py {csparcfile} {temp_star}"], shell=True)

    # Get particle stack file path
    starfile = StarFile(temp_star)
    particles_df = starfile.particles_df
    

    # Execute Python script to change image path
    change_image_path(temp_star, outputfile)
    
    # Remove temporary star file
    os.remove(temp_star)
   
    # change export job links with +'s'
    for j_dir in os.listdir():
        if not os.path.isdir(j_dir) or not j_dir.startswith('J'):
            continue
        extract_path = os.path.join(j_dir, 'extract')
        if not os.path.isdir(extract_path):
            continue

        for file in os.listdir(extract_path):
            if os.path.isfile(os.path.join(extract_path, file)):
                old_path = os.path.join(extract_path, file)
                if old_path.endswith('mrc'):
                    new_path = os.path.join(extract_path, file + 's')
                elif old_path.endswith('mrcs'):
                    new_path = old_path
                os.rename(old_path, new_path)
    print("Remember to run relion_stack_create xxx to get new stack file!")

### In 4.2.1 version, the output star file has *rlnImageName at first column, *rlnMicrographName at second column.
def change_image_path(file_1, file_2):
    starfile = StarFile(file_1)
    particles_df = starfile.particles_df
    current_path = os.path.join(os.getenv('PWD'), '')
    for line_number in range(particles_df.shape[0]):
        image = particles_df.loc[line_number, '_rlnImageName']
        if image.endswith('.mrc'):
            image = image + 's'
        
        image = image.replace('>', '')
        image = image.replace('@', '@'+current_path) 
        particles_df.loc[line_number, '_rlnImageName'] = image
		
		

        
        if '_rlnMicrographName' in particles_df.columns:
            micrograph = particles_df.loc[line_number, '_rlnMicrographName']
            micrograph = micrograph.replace('>', '')
            micrograph = os.path.join(current_path, micrograph)
            particles_df.loc[line_number, '_rlnMicrographName'] = micrograph

    
    # Update the dataframe
    starfile.particles_df = particles_df
    # Write new file
    starfile.write(file_2)

if __name__ == "__main__":
    main()
