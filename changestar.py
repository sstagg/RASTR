#! /usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import random
from starparse import StarFile
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
	parser.add_argument('--remove', action='store', default=None, dest='remove')
	return parser.parse_args(args)

# mapping from command line arguments to star file column names
arg_to_column = {
	"micrograph": "_rlnMicrographName",
	"image": "_rlnImageName",
	"mag": "_rlnMagnification",
	"Rot": "_rlnAngleRot",
	"Tilt": "_rlnAngleTilt",
	"Psi": "_rlnAnglePsi",
	"X": "_rlnOriginXAngst",
	"Y": "_rlnOriginYAngst",
}

def change_star_file_values(args):
	parsed_args = parserInput(args)
	star_file = StarFile(parsed_args.input)

	if parsed_args.remove is not None:
		columns_to_remove = [star_file.particles_df.columns[int(i)-1] for i in parsed_args.remove.split(',')]
		star_file.particles_df = star_file.particles_df.drop(columns=columns_to_remove)

	for arg, value in vars(parsed_args).items():
		if value is not None and arg in arg_to_column:
			column = arg_to_column[arg]

			if arg == 'image':
				new_image_name = star_file.particles_df[column].str.replace(r'@.*$', f'@{value}', regex=True)
				star_file.particles_df[column] = new_image_name

			elif arg == 'micrograph':
				star_file.particles_df[column] = value

			elif arg in ['Rot', 'Tilt', 'Psi', 'X', 'Y']:
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


	star_file.write(parsed_args.output)

# run the script
if __name__ == "__main__":
	import sys
	change_star_file_values(sys.argv[1:])

