#! /usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import random
from starparse import StarFile
import sys
def parserInput(args):
	parser = argparse.ArgumentParser(description='Modify Relion star file values.')
	parser.add_argument('-i', '--i', action='store', default=None, dest='input')
	parser.add_argument('-o', '--o', action='store', default=None, dest='output')
	parser.add_argument('-mic', '--micrograph', action='store', default=None, dest='micrograph')
	parser.add_argument('-img', '--image', action='store', default=None, dest='image')
	parser.add_argument('-mag', '--magnification', action='store', default=None, dest='mag')
	parser.add_argument('-rot','--anglerot',action='store',default=None,dest='Rot')
	parser.add_argument('-tilt','--angletilt',action='store',default=None,dest='Tilt')
	parser.add_argument('-psi','--anglepsi',action='store',default=None,dest='Psi')
	parser.add_argument('-x','--shiftx',action='store',default=None,dest='X')
	parser.add_argument('-y','--shifty',action='store',default=None, dest='Y')
	parser.add_argument('-n','--norm',action='store',default=None, dest='norm')
	parser.add_argument('--remove', action='store', default=None, dest='remove', help="comma separated list of columns to remove (_rlnAnglePsi)")
	parser.add_argument('--optic', action='store', default=1, dest='opticgroup')
	parser.add_argument('--tilt-range', nargs=2, default=False, type=float, metavar=('MIN', 'MAX'), help='Filter tilt within the specified range')
	parser.add_argument('--sort', action='store_true', default=False, dest='sort', help='Sort particles based on slice number, only work when there is just one particle stack file')
	return parser.parse_args(args)





def filter_tilt_range(star_file, tilt_min, tilt_max):
	star_file.particles_df  = star_file.particles_df[
        (star_file.particles_df["_rlnAngleTilt"] >= tilt_min) &
        (star_file.particles_df["_rlnAngleTilt"] <= tilt_max)
    ]
	
	return star_file




# mapping from command line arguments to star file column names


def change_star_file_values(args):
	arg_to_column = {
		"micrograph": "_rlnMicrographName",
		"image": "_rlnImageName",
		"mag": "_rlnMagnification",
		"Rot": "_rlnAngleRot",
		"Tilt": "_rlnAngleTilt",
		"Psi": "_rlnAnglePsi",
		"X": "_rlnOriginXAngst",
		"Y": "_rlnOriginYAngst",
		"class" : "_rlnClassNumber",
		"opticgroup" : "_rlnOpticsGroup",
		"norm" : "_rlnNormCorrection",
		}
	parsed_args = parserInput(args)
	star_file = StarFile(parsed_args.input)

	if parsed_args.remove is not None:
		columns_to_remove = parsed_args.remove.split(',')
		star_file.particles_df = star_file.particles_df.drop(columns=columns_to_remove)

	for arg, value in vars(parsed_args).items():
		if value is not None and arg in arg_to_column:
			column = arg_to_column[arg]

			if arg == 'image':
				new_image_name = star_file.particles_df[column].str.replace(r'@.*$', f'@{value}', regex=True)
				star_file.particles_df[column] = new_image_name

			elif arg == 'micrograph':
				star_file.particles_df[column] = value

			elif arg == 'opticgroup':
				star_file.particles_df[column] = value

			elif arg in ['Rot', 'Tilt', 'Psi', 'X', 'Y', 'norm']:
				if value.startswith('r'):
					if arg in ['Rot', 'Tilt', 'Psi']:
						star_file.particles_df[column] = np.random.uniform(0, float(value[1:]), star_file.particles_df.shape[0])
					elif arg in ['X', 'Y']:
						star_file.particles_df[column] = np.random.uniform(-float(value[1:]), float(value[1:]), star_file.particles_df.shape[0])
				elif value.startswith('+'):
					star_file.particles_df[column] += float(value[1:])
				elif value.startswith('-'):
					star_file.particles_df[column] -= float(value[1:])
				elif value.startswith('*'):
					star_file.particles_df[column] *= float(value[1:])
				elif value.startswith('/'):
					star_file.particles_df[column] /= float(value[1:])
				else:
					star_file.particles_df[column] = float(value)

	if parsed_args.tilt_range:
		tilt_min, tilt_max = parsed_args.tilt_range
		star_file = filter_tilt_range(star_file, tilt_min, tilt_max)

	if parsed_args.sort:
		## function not done, detect if more than one particles stack
		if check_only_one_stack():
			print("two stack file detected, cannot sort, exiting")
			sys.exit()


		particles_df = star_file.particles_df
		particles_df = particles_df.sort_values("_rlnImageName", key=lambda x: x.apply(extract_image_int))
		star_file.particles_df = particles_df		

	star_file.write(parsed_args.output)

def check_only_one_stack():
	return False

def extract_image_int(values):
	return int(values.split('@')[0])

# run the script
if __name__ == "__main__":
	import sys
	change_star_file_values(sys.argv[1:])

