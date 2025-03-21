#!/usr/bin/env python
"""
RELION Analysis Script

Analyzes RELION particle data to show resolution trends and class distributions across iterations.
Usage: Navigate to the job directory and run ./analyze_relion.py
"""

import os
import argparse
from matplotlib import pyplot as plt


def count_iterations(files):
	"""Count the total number of iterations based on model.star files."""
	return sum(1 for file in files if file.endswith('model.star'))


def get_project_info(files):
	"""Extract initial project information including root name and class count."""
	project_info = {}
	
	# Find root name from optimiser file
	for star_file in files:
		if star_file.endswith('optimiser.star'):
			with open(star_file, 'r') as file:
				for line in file:
					if line.strip().startswith('_rlnOutputRootName'):
						root_name = line.split()[1]
						project_info['rootname'] = root_name
						break
			break
	
	# Get class count from initial model file
	model_file = f"{root_name.split('/')[-1]}_it000_model.star"
	with open(model_file, 'r') as file:
		for line in file:
			if line.strip().startswith('_rlnNrClasses'):
				project_info['class_number'] = int(line.split()[1])
				break
				
	return project_info


def analyze_column_positions(filename):
	"""Determine column positions for resolution and distribution data."""
	position = {}
	
	with open(filename, 'r') as file:
		for line in file:
			if not line.strip():
				continue
				
			if line.split()[0] == '_rlnClassDistribution':
				position['percent'] = int(line.split()[1][1:]) - 1
			elif line.split()[0] == '_rlnEstimatedResolution':
				position['resolution'] = int(line.split()[1][1:]) - 1
			
			if line.strip() == 'data_model_class_1':
				break
				
	return position


def plot_results(iteration_count, resolutions, distributions):
	"""Generate and display plots for resolution and distribution data."""
	# Plot resolutions
	plt.figure(figsize=(10, 6))
	for class_name, res_values in resolutions.items():
		plt.plot(range(iteration_count), res_values)
	
	plt.legend(list(resolutions.keys()), loc='upper left')
	plt.xlabel('Iteration')
	plt.ylabel('Resolution (Å)')
	plt.title('Resolution Progression by Class')
	plt.grid(True, linestyle='--', alpha=0.7)
	plt.tight_layout()
	plt.show()
	
	# Plot distributions
	plt.figure(figsize=(10, 6))
	for class_name, dist_values in distributions.items():
		plt.plot(range(iteration_count), dist_values)
	
	plt.legend(list(distributions.keys()), loc='upper left')
	plt.xlabel('Iteration')
	plt.ylabel('Class Distribution (%)')
	plt.title('Particle Distribution by Class')
	plt.grid(True, linestyle='--', alpha=0.7)
	plt.tight_layout()
	plt.show()


def analyze_data():
	"""Main analysis function to extract and process RELION data."""
	files = os.listdir('.')
	iteration_count = count_iterations(files)
	project_info = get_project_info(files)
	
	root_name = project_info['rootname']
	class_count = project_info['class_number']
	
	# Create dictionary for each class
	class_names = [f'class{i+1:03d}' for i in range(class_count)]
	resolutions = {name: [] for name in class_names}
	distributions = {name: [] for name in class_names}
	
	# Get column positions for data extraction
	model_base = root_name.split('/')[-1]
	positions = analyze_column_positions(f"{model_base}_it000_model.star")
	
	# Process each iteration
	for i in range(iteration_count):
		iteration = f"it{i:03d}"
		model_file = f"{model_base}_{iteration}_model.star"
		
		with open(model_file, 'r') as file:
			for line in file:
				if not line.strip():
					continue
					
				tokens = line.split()
				if tokens and tokens[0].startswith(root_name):
					# Extract class name from filename
					class_name = tokens[0][(len(root_name)+7):-4]
					
					# Store resolution and distribution values
					resolutions[class_name].append(float(tokens[positions['resolution']]))
					distributions[class_name].append(float(tokens[positions['percent']]))
				
				if line.strip() == 'data_model_class_1':
					break
	
	return iteration_count, resolutions, distributions


def main():
	"""Main function to run the analysis and display results."""
	# Parse command line arguments
	parser = argparse.ArgumentParser(description='Analyze RELION particle data.')
	parser.add_argument('--save', action='store_true', help='Save plots instead of displaying them')
	args = parser.parse_args()
	
	# Get analysis data
	iteration_count, resolutions, distributions = analyze_data()
	
	# Plot results
	plot_results(iteration_count, resolutions, distributions)


if __name__ == '__main__':
	main()