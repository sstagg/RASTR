#!/usr/bin/env python
### Script to export particles from cryosparc to relion
### The output is star file with correct imagename, you should do relion_stack_create on your own.
### Usage, cd to export path with file eg. J185_particles_selected_exported.cs first.
### Usage, then run ./csexport.py J185_particles_selected_exported.cs J185_particles_selected_exported.star

import os
import subprocess
import sys
from src.common.starparse import StarFile

def main():
    # Define temporary paths
    temp_star = "tempt.star"
    
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
    for line_number in range(particles_df.shape[0]):
        image = particles_df.loc[line_number, '_rlnImageName']
    
            
        # old csparc2star output has '>' at the beginning of the image name
        image = image.replace('>', '')
        image = image.replace('@', '@'+current_path) 

        if image.endswith('.mrc'):
            original_image = image.split('@')[1]
            image = image + 's'
            image_path = original_image + 's'

            particles_df.loc[line_number, '_rlnImageName'] = image
            # check if image_path exist
            if not os.path.exists(image_path):
                # create symbolic link to the original image
                os.symlink(original_image, image_path)
            else:
                continue
                

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
