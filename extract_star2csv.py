#! /usr/bin/env python
from starparse import StarFile
import pandas as pd
import sys
import argparse

def parserInput():
	parser = argparse.ArgumentParser(description='Extract Relion star file values.')
	parser.add_argument('-i', '--i', action='store', default=None, dest='input')
	parser.add_argument('-o', '--o', action='store', default=None, dest='output')
	parser.add_argument('-mic', '--micrograph', action='store_true', default=None, dest='micrograph')
	parser.add_argument('-img', '--image', action='store_true', default=None, dest='image')
	parser.add_argument('-mag', '--magnification', action='store_true', default=None, dest='mag')
	parser.add_argument('-rot','--anglerot',action='store_true',default=None,dest='Rot')
	parser.add_argument('-tilt','--angletilt',action='store_true',default=None,dest='Tilt')
	parser.add_argument('-psi','--anglepsi',action='store_true',default=None,dest='Psi')
	parser.add_argument('-x','--shiftx',action='store_true',default=None,dest='X')
	parser.add_argument('-y','--shifty',action='store_true',default=None, dest='Y')
	parser.add_argument('--optic', action='store_true', default=None, dest='opticgroup')
	parser.add_argument('--loglikeli', action='store_true', default=None, dest='loglikeli')
	return parser.parse_args()




def main():
	parsed_args = parserInput()
	input_file = parsed_args.input
	starfile = StarFile(input_file)
	
	#loop all args
	columns_to_keep = []
	for arg in parsed_args.__dict__:
		if parsed_args.__dict__[arg] is not None and arg in arg_to_column:
			columns_to_keep.append(arg_to_column[arg])
	
	#drop dataframes
	starfile.particles_df = starfile.particles_df[columns_to_keep]

	starfile.write(parsed_args.output)

if __name__=="__main__":
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
		"loglikeli" : "_rlnLogLikeliContribution"
		
		}
	main()