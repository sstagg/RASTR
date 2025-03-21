#! /usr/bin/env python
from common.starparse import StarFile
import sys
import numpy as np
import time
import pandas as pd
import argparse

def within180(lis):
    return (lis + 360)%180

def parse_arguments():
    parser = argparse.ArgumentParser(description='Compare psi angles between two star files')
    parser.add_argument('file1', help='First star file')
    parser.add_argument('file2', help='Second star file')
    parser.add_argument('--output', '-o', default='large_psi_diff.star', dest='output', action='store',
                       help='Output star file for particles with psi difference > 30 (default: large_psi_diff.star)')
    parser.add_argument('--write', '-w', action='store_true', dest='write',
                       help='Write output star file')
    return parser.parse_args()

def main():
    args = parse_arguments()
    star_file_1 = StarFile(args.file1)
    star_file_2 = StarFile(args.file2)


    particle_df_1 = star_file_1.particles_df
    particle_df_2 = star_file_2.particles_df
    if len(particle_df_1) != len(particle_df_2):
        print('star files not matching, using image name to match particles which may be slow')
    #print( particle_df_1)
    dif = []
    dif_all = []
    dif_large = []
    
    for linenumber, row in particle_df_1.iterrows():
        if row['_rlnImageName'] != particle_df_2['_rlnImageName'][linenumber]:
            print('star file not matching')

        psi_1 = float(row['_rlnAnglePsi'])
        psi_2 = float(particle_df_2['_rlnAnglePsi'][linenumber])

        psi_1 = within180(psi_1)
        psi_2 = within180(psi_2)
        dif_o = abs(psi_1 - psi_2)
        
        if dif_o > 90:
            dif_o = 180-dif_o
        
        dif_all.append(dif_o)
        
        if dif_o <= 30:
            dif.append(dif_o)
        else:
            dif_large.append(row)
    


    if args.write:
        particle_df_3 = pd.DataFrame(dif_large, columns=particle_df_1.columns)
        star_file_1.particles_df = particle_df_3
        star_file_1.write(args.output)

    print(f"Particles with psi difference ≤ 30°: {len(dif)}")
    dif = np.array(dif)
    dif_all = np.array(dif_all)

    rmsd = np.sqrt(((dif) ** 2).mean())
    print(f"RMSD (≤ 30°): {rmsd}")

    print(f"Total particles: {len(dif_all)}")
    rmsd = np.sqrt(((dif_all)**2).mean())
    print(f"RMSD (all): {rmsd}")

if __name__ == "__main__":
    main()




