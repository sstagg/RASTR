#! /usr/bin/env python

### Analyse function to show particles, resolution trends among iterations.
### Usage, cd to the job directory, then ./analyzerelion.py

import os
from matplotlib import pyplot

### get total iterations
def count_its(files):
	return sum(1 for file in files if file.endswith('model.star'))


### get information for project
def initial(files):
	result = {}
	for sfile in files:
		### get into any optimiser file
		if sfile.endswith('optimiser.star'):
			with open(sfile,'r') as a:
				lines = a.readlines()
			for line in lines:
				if len(line) > 18 and line[:18] == '_rlnOutputRootName':
					### get output root name for relion project, typically 'run'
					name = line.split()[1]
					result['rootname'] = name
					break
			break
	# Name typically looks as Class3D/job001/run
	with open(name.split('/')[-1] + '_it000_model.star','r') as a:
		lines = a.readlines()

	### get total classes
	for line in lines:
		if len(line)>13 and line[:13] == '_rlnNrClasses':
				classes = line.split()[1]
				result['class_number'] = int(classes)
				break
	return result

		


### pyplot show function, show resolution first, then distribution
def show(itnumber, res=None, percent=None):
	for i in res:
		pyplot.plot( range(itnumber), res[i] )
	pyplot.legend([i for i in res], loc='upper left')
	pyplot.xlabel('iteration')
	pyplot.ylabel('resolution')
	pyplot.show()
	for i in res:
		pyplot.plot( range(itnumber), percent[i] )
	pyplot.legend( [i for i in res], loc='upper left')
	pyplot.xlabel('iteration')
	pyplot.ylabel('distribution')
	pyplot.show()



### open it000_model file to get position of resolution and distribution
def analyze_position( filename ):
	with open(filename,'r') as a:
		lines = a.readlines()
	position = {}
	for line in lines:
		if len(line) > 14:
			if line.split()[0] == '_rlnClassDistribution':
				position['percent'] = int(line.split()[1][1:])-1
			elif line.split()[0] == '_rlnEstimatedResolution':
				position['resolution'] = int(line.split()[1][1:])-1
			
			if line == 'data_model_class_1\n':
				break
	return position

def analyze():
	files = os.listdir('.')
	itnumber = count_its(files)
	initials = initial(files)
	rootname = initials['rootname']
	classes = initials['class_number']	

	### get a list of strings of classname
	classnames = ['class' + str('%03i' %(i+1)) for i in range(classes)]
	resolutions = {i: [] for i in classnames}
	percent = {i: [] for i in classnames}
	
	position = analyze_position (rootname.split('/')[2]+'_it000_model.star')

	### append model file stepwise
	for i in range(itnumber):
		iteration = 'it'+str('%03i' %i)
		with open(rootname.split('/')[2]+'_'+iteration+'_model.star','r') as a:
			lines = a.readlines()
		
		for line in lines:
			if len(line) > len(rootname):
				words = line.split()
				if words[0][:len(rootname)] == rootname:
					resolutions[words[0][(len(rootname)+7):-4]].append(float(words[position['resolution']]))
					percent[words[0][(len(rootname)+7):-4]].append(float(words[position['percent']]))
				if line == 'data_model_class_1\n':
					break
	return itnumber, resolutions, percent


def main():
	# Get lists of resolutions, percentages
	itnumber, res, percent = analyze()

	# Plot
	show( itnumber, res, percent)


if __name__ == '__main__':
	main()
