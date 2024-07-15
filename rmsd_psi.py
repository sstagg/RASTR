#! /usr/bin/env python
from starparse import StarFile

import sys
import numpy as np
import time
import pandas as pd
def within180(lis):
	lis =  (lis + 360)%180
	return lis
file_1 = sys.argv[1]
file_2 = sys.argv[2]

star_file_1 = StarFile( file_1 )
star_file_2 = StarFile( file_2 )

particle_df_1 = star_file_1.particles_df
particle_df_2 = star_file_2.particles_df


#create an empty pandas dataframe
particle_df_3 = pd.DataFrame()
dif = []

# iterate through the dataframe
for linenumber, row in particle_df_1.iterrows():
	
	psi_1 = float(row['_rlnAnglePsi'])
	psi_2 = float(particle_df_2['_rlnAnglePsi'][linenumber])


	psi_1 = within180(psi_1)
	psi_2 = within180(psi_2)
	dif_o  = abs( psi_1 - psi_2)
	if dif_o > 90:
		dif_o = 180-dif_o
	
	
	
	if dif_o <= 50:
		dif.append(dif_o)
	#else:
		#append row to the dataframe
		#particle_df_3 = particle_df_3.append(row, ignore_index=True)
	else:
		print(linenumber)


#star_file_1.particles_df = particle_df_3
#star_file_1.write('output.star')

print(len(dif))
dif = np.array(dif)
	


rmsd = np.sqrt(((dif) ** 2).mean())
print(rmsd)





