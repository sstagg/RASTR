#! /usr/bin/env python

from common.starparse import StarFile
import pandas as pd
from matplotlib import pyplot as plt
import sys


def main():
    reader = pd.read_csv(sys.argv[1])
    repeat_list = [int(x) for x in reader.iloc[:,0].tolist()]
    #print (repeat_list)
    
    #read star file
    star_file = StarFile(sys.argv[2])
    particles_df = star_file.particles_df
    optics_df = star_file.optics_df
    pixel_size = float(optics_df.loc[0,'_rlnImagePixelSize'])
    #print(len(particles_df))
    index_to_delete = []
    for line_number, particle in particles_df.iterrows():
        image = particle['_rlnImageName'].split('@')
        if int(image[0]) in repeat_list:
            index_to_delete.append(line_number)
    
    particles_df = particles_df.drop(index_to_delete)

    #print (len(particles_df))
    #print (len(index_to_delete))

    star_file.particles_df = particles_df
    star_file.write(sys.argv[3])
    
    pass


if __name__ == '__main__':
    main()
