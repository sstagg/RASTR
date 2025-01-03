#! /usr/bin/env python
#### script for RASTR by Ruizi
#### this script takes a pre-made mask. 
#### a typical run command: ./RASTR_02.py --star_in run_it005_data.star  --model run_it005_class001azavg.mrc  --angpix 2.02  -k -o rootname -ma spheremask.mrc  -al 0,90,180,270  --pad 3


import numpy as np
import sys
import os, errno
import subprocess
import argparse
import logging
import copy
from multiprocessing import Process
import mrcfile
import cupy as cp
from cupyx.scipy.ndimage import gaussian_filter, rotate
from starparse import StarFile
##Functions

#Makes a directory
def make_tmp(directory):
	try:
		os.makedirs(directory)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

#run a command, argument is a string
def run_command(command):
	print (command)
	process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
	output, error = process.communicate()

#flips the values in an array, only values between 0 and 1
def flip(maskvolume):
	if maskvolume.max() <= 1 and maskvolume.min() >= 0:
		return 1 - maskvolume
	else:
		print ('not binary mask')
	
### using scipy to rotate a 3d matrix with rot, tilt, psi angle, rotation order verified using cisTEM generate.
def volumerotation( volume, rot=0, tilt=0, psi=0, order=1):
	newvolume = copy.deepcopy(volume)
	if rot != 0:
		newvolume = rotate( newvolume, rot, axes=[1,2], order=order, reshape=False)
	if tilt != 0:	
		newvolume = rotate( newvolume, -tilt, axes=[0,2], order=order, reshape=False)
	if psi != 0:
		newvolume = rotate( newvolume, psi, axes=[1,2], order=order, reshape=False)
	return newvolume

### parse the top part of star file, return a dictionary with different column and number. eg paracolumn['rot']=8, then the 9th value in a row represent rot angle.
def column(filename):
	if filename[-4:] == 'star':
		paracolumn = {}
		fobj = open(filename,'r')
		lines = fobj.readlines()
		fobj.close()
		datastart = False
		for line in lines:
			words = line.split()
			if line == 'data_particles\n' or line == 'data_\n':
				datastart = True
			if len(words) == 2 and datastart: 
				if words[0] == '_rlnMicrographName':
					paracolumn['mic'] = int(words[1][1:])-1
				if words[0] == '_rlnImageName':
					paracolumn['image'] = int(words[1][1:])-1
				if words[0] == '_rlnMagnification':
					paracolumn['mag'] = int(words[1][1:])-1
				if words[0] == '_rlnAngleRot':
					paracolumn['rot'] = int(words[1][1:])-1
				if words[0] == '_rlnAngleTilt':
					paracolumn['tilt'] = int(words[1][1:])-1
				if words[0] == '_rlnAnglePsi':
					paracolumn['psi'] = int(words[1][1:])-1
				if words[0] == '_rlnOriginX':
					paracolumn['shx'] = int(words[1][1:])-1
				if words[0] == '_rlnOriginY':
					paracolumn['shy'] = int(words[1][1:])-1
### this is copied from other scirpt for modify meta data, par may not be needed here.
	elif filename[-3:] == 'par':
		paracolumn = {}
		paracolumn['mag'] = 6
		paracolumn['psi'] = 1
		paracolumn['tilt'] = 2
		paracolumn['rot'] = 3
		paracolumn['shx'] = 4
		paracolumn['shy'] = 5
	else:
		print ("unknown file type")
	return paracolumn

### funcition for subtration, utilizing relion_project --subtract_exp
### mrcs_subtracted is a global dictionary in main, angle is its key
def subtraction(subtractionmodel, tempangle, angpix, starfilename, outputfilename, pad):
	relion_subtraction = 'relion_project --i '+subtractionmodel+'  --o '+outputfilename+'  --angpix '+angpix+' --ang '+starfilename+' --ctf --pad '+str(pad)+' --subtract_exp'
	logging.info('angle '+str(tempangle)+': '+relion_subtraction)
	run_command(relion_subtraction)
	logging.info('angle '+str(tempangle)+': done!')


### logger function, inherited from Peter
def setupLogger(results):

	logging.basicConfig(level=logging.DEBUG,
	                    format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
	                    datefmt  = '%m-%d %H:%M',
	                    filename = results.output_rootname+'.log',
	                    filemode = 'w')
	#define a handler which writes INFO messages or higher to the sys.stderr
	console = logging.StreamHandler()
	console.setLevel(logging.INFO)
	#set format for which is simpler for console use
	formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
	#tell the handler to use this format
	console.setFormatter(formatter)
	# add the handler to the root logger
	logging.getLogger('').addHandler(console)
	
	
	#seperate loggers for different section of program
	
	logger1 = logging.getLogger('Making models         ')
	logger2 = logging.getLogger('Subtraction           ')
	logger3 = logging.getLogger('Extraction            ')
	logger4 = logging.getLogger('Creating final files  ')
	#output the initial values

	commandline=' '
	logging.info(commandline.join(sys.argv))
	logging.info('star in          = '+ results.star_in)
	logging.info('model to delete  = '+ results.model)
	logging.info('angpix           = '+ results.angpix)
	logging.info('output rootname  = '+ results.output_rootname)
	logging.info('keep scratch     = '+ str(results.keep_scratch))
	logging.info('scratch folder   = '+ results.scratch)	

	return logger1, logger2, logger3, logger4

def readslice( image ):
	with mrcfile.mmap(image[1], mode='r') as imagestack:
		zslice = int(image[0]) - 1
		imagearray = imagestack.data[zslice]
	imagearray_new = cp.asarray(imagearray)
	return imagearray_new

def relionmask(angle, whitemaskfile, maskprojectionfile, angpix, gauss, initial_star_in, mrcs_masked):
	logging.info('starting extraction for angle: '+ str(angle))
	logging.info('projecting masks')
	maskproj_run='relion_project --i '+whitemaskfile+'  --o '+maskprojectionfile+'  --angpix '+angpix+' --ang '+initial_star_in+' --pad 3'
	run_command(maskproj_run)
	logging.info('done projection masks for angle: '+str(angle))

	newfile=mrcs_masked
	
	mask_starfile = StarFile( maskprojectionfile+'.star' )
	mask_particles_df = mask_starfile.particles_df
	particle_starfile = StarFile( initial_star_in )
	particles_df = particle_starfile.particles_df

	with mrcfile.open(maskprojectionfile+'.mrcs') as mrc:
		mrcsshape = mrc.data.shape
	print (mrcsshape)
	with mrcfile.new_mmap(newfile, shape=mrcsshape, mrc_mode=2, overwrite=True) as mrcobj:
		for linenumber, particle in particles_df.iterrows():
			image = particle['_rlnImageName'].split('@')
			
			mask_image = mask_particles_df.iloc[linenumber]['_rlnImageName'].split('@')
			i = int(mask_image[0])
			singleslice = readslice(image)
			maskslice = readslice(mask_image)

			maskslice = gaussian_filter(maskslice,sigma=gauss)
			threshold = 8
			maskslice[maskslice<=threshold] = 0.0
			maskslice[maskslice>threshold] = 1.0

			maskedslice = maskslice*singleslice
			mrcobj.data[i-1]=maskedslice.get()
				
	logging.info('done extraction for angle: '+str(angle))



def parseOptions():

	parser = argparse.ArgumentParser()
	
	parser.add_argument('-s','--star_in', action='store', dest='star_in',
						help='      initial star file for the aligned model')
	parser.add_argument('-m','--model', action='store', dest='model',
						help='      aligned model to subtract')
	parser.add_argument('-a','--angpix', action='store', dest='angpix',
						help='      angstroms per pixel')
	parser.add_argument('-r','--radius', type=int, action='store', dest='radius',
						help='      Size of sphere to reconstruct, pixels, default is 3/16 of box size')
	parser.add_argument('-x','--x_start', type=int, action='store', dest='x_start',
						help='      center of sphere to mask in x at phi=0, pixels, IMOD coordinates, default is 3/4 of box size')
	parser.add_argument('-n','--n_spheres', type=int, action='store', dest='n_spheres', default=4,
						help='      Number of spheres to mask around axis, convenient if 360 is divisable by n_spheres, defualt is 9')
	parser.add_argument('-t','--tube_radius', type=int, action='store', dest='trad',
						help='      radius of the membrane to be subtracted, this will remove the tube out to a certain radius. Used for decorations')
	parser.add_argument('-c','--center', action='store_true', default=False, dest='center', 
						help='      center the masked area in a new smaller box')
	parser.add_argument('-o', '--output', action='store', dest='output_rootname', default='RASTR_particles',
						help='     output rootname for star/mrcs, default is RASTR_particles')
	parser.add_argument('-k', '--keep_scratch', action='store_true', default=False, dest='keep_scratch', 
						help='      keep scratch folder, default is False, only final star, mrcs, and reconstruction output on default')
	parser.add_argument('-pad','--pad', action='store',dest='pad',default=2,
					 	help='	  padding for relion_project, default is 2')
	parser.add_argument('-or','--order',action='store',dest='order',default=3,
					 	help='	  order for interpolation, default is 3')
	parser.add_argument('-f','--test', action='store_true', default=False, dest='test_TF', 
						help='      option for testing, dont use, was needed during setup')
	parser.add_argument('-b','--both',action='store_true', default=False, dest='both',
						help='	do both masking and center. End with two stacks, one with section masked keeping same box size as original and the second with it centered in a smaller box')	
	parser.add_argument('-g','--gauss', type=int, default=5, dest='gauss', 
						help='	Sigma for gaussian edge to add to sphere. Default = 5')	
	parser.add_argument('-al','--anglelist',type=str,action='store',default='0,90,180,270', dest='anglelist',
						help=' input list of angles, separate by comma')
	parser.add_argument('-ma','--mask', type=str, action='store',default=None,dest='mask',
						help='input mask for refinement, make the region to be kept 1')
	parser.add_argument('-sin','--single',action='store_true',default=False,dest='single', 
						help='single thread mode, default is multi-thread mode')
	parser.add_argument('--phi', type=float, action='store', default=None, dest='phi',
						help='starting phi angle for rotation, default is using star file angles')
	results = parser.parse_args()
	
	#make sure all the required arguments are here

	if results.star_in == None or results.model == None or results.angpix == None:
		print ('Required arguments: --star_in <star file> --model <model file> -- angpix <angpix>; use -h or --help to see all options')
		sys.exit()
	
	return results
		
#Main program
def main():
	#Input arguments
	results = parseOptions()
	
	#Create scratch folder
	scratch = 'scratch'
	scratch_num = 2
	if os.path.isdir(scratch) is True:
		while os.path.isdir(scratch) is True:
			scratch = 'scratch'+str(scratch_num)
			scratch_num+=1
	results.scratch = scratch
	make_tmp(scratch)


	#get box size first because some defaults depend on it:
	input_model = mrcfile.read(results.model)
	input_model = cp.asarray(input_model.data)
	box_size = int(input_model.shape[0])

	#set defaults
	pad = results.pad
	order = results.order
	output_rootname = results.output_rootname
	angles=results.anglelist.split(',')
	for i in range(len(angles)):
		angles[i]=float(angles[i])

	subtractionmodels={}
	blackmasks={}
	whitemasks={}
	mrcs_subtracted={}
	maskprojection={}
	mrcs_masked={}
	for angle in angles:
		subtractionmodels[angle] = scratch+'/'+'subtractionmodel'+str(angle)+'.mrc'
		blackmasks[angle] = scratch + '/blackmask' + str(angle) + '.mrc'
		whitemasks[angle] = scratch + '/whitemask' + str(angle) + '.mrc'
		mrcs_subtracted[angle] = scratch + '/' + results.output_rootname + '_' + str(angle) + '_subtracted'
		maskprojection[angle]=scratch+'/maskprojection_'+str(angle)
		mrcs_masked[angle]=mrcs_subtracted[angle]+'_masked.mrcs'
	

	#create a log file
	logger1,logger2,logger3,logger4=setupLogger(results)

	#####             #####
	##### START RASTR #####
	#####             #####
	logging.info('Starting RASTR')

	#Default values to be printed in log file	
	if results.phi == None:
		initial_star_in = results.star_in
	else:
		logging.info('phi angle provided, changing star file')
		initial_star_in = output_rootname+'_tmp_phi.star'
		from changestar import parserInput, change_star_file_values
		change_star_file_values(['--i', str(results.star_in), '--o', initial_star_in, '-rot', str(results.phi)])
		logging.info('Done editing star file')

	whitemask = mrcfile.read(results.mask)
	whitemask = cp.asarray(whitemask.data)

	### making subraction models
	logging.info('making subtraction models')
	for angle in angles:
		subtractionfile = subtractionmodels[angle]
		b_maskfile = blackmasks[angle]
		w_maskfile = whitemasks[angle]
		
		whitemaskvolume = volumerotation(whitemask, rot=angle, tilt=0, psi=0, order=order)
		whitemaskvolume = gaussian_filter(whitemaskvolume,sigma=3)
		whitemaskvolume[whitemaskvolume>1] = 1.0
		whitemaskvolume[whitemaskvolume<0] = 0.0
		blackmaskvolume = flip(whitemaskvolume)

		mrcfile.write(w_maskfile, data=whitemaskvolume.get())
		mrcfile.write(b_maskfile, data=blackmaskvolume.get())

		subtractionvolume = blackmaskvolume*input_model
		mrcfile.write(subtractionfile, data=subtractionvolume.get())
		logger1.info(subtractionfile)
	logging.info('Done making subtraction models')



	####subtraction
	if results.single:
		for angle in angles:
			outputfilename = mrcs_subtracted[angle] 
			subtraction(subtractionmodels[angle], angle, results.angpix, initial_star_in, outputfilename, pad)
	else:
		processes = []
		for angle in angles:
			outputfilename = mrcs_subtracted[angle]
			p = Process(target=subtraction, args=(subtractionmodels[angle], angle, results.angpix, initial_star_in, outputfilename, pad))
			p.start()
			processes.append(p)

		for p in processes:
			p.join()
	logging.info('Done with subtraction')
		
	#### mask target region
	####  speed of this step is mainly affected by order number. 
	#### It is within range 0-5, the larger the longer for scipy.ndimage.rotate to run
	logging.info('')
	logger3.info('Extracting particles')


	paracolumn=column(initial_star_in)
	
#angle, whitemaskfile, maskprojectionfile, angpix, initial_star_in, mrcs_masked
	if results.single:
		for angle in angles:
			relionmask(angle, whitemasks[angle], maskprojection[angle], results.angpix, results.gauss, initial_star_in, mrcs_masked[angle])
	else:
		processes = []
		for angle in angles:
			p=Process(target=relionmask,args=(angle, whitemasks[angle], maskprojection[angle], results.angpix, results.gauss, initial_star_in, mrcs_masked[angle],))
			p.start()
			processes.append(p)
		for p in processes:
			p.join()
	logging.info('Done with masking')
	

	### join star file together
	outputstar=open(output_rootname+'_tmp.star','w')
	firstanglestar=open(mrcs_subtracted[angles[0]]+'.star','r')
	lines=firstanglestar.readlines()
	firstanglestar.close()
	### this loop for writing header parts


	opticstart=False
	particlesstart=False
	for line in lines:
		words=line.split()
		if line=='data_optics\n':
			opticstart=True
		if line=='data_particles\n':
			particlesstart=True
	#### when in optics line, write them directly into new star file
		if opticstart and not particlesstart:
			outputstar.write(line)
	### when in particle column part, write column part
		elif particlesstart and len(words)<=2:
			outputstar.write(line)
	### this loop for writing data lines
	for angle in angles:
		fileobj=open(mrcs_subtracted[angle]+'.star','r')
		lines=fileobj.readlines()
		fileobj.close()
		opticstart=False
		particlesstart=False
		for line in lines:
			try:
				if line[0]=='#':
					continue
			except:
				pass
			words=line.split()
			if line=='data_optics\n':
				opticstart=True
			if line=='data_particles\n':
				particlesstart=True
			### skip optic part and particle column part, write only data lines
			if opticstart and particlesstart and len(words)>2:
			#print words
				words[paracolumn['rot']]=str((float(words[paracolumn['rot']])+angle)%360.0)
				words[paracolumn['image']]=words[paracolumn['image']][:-5]+'_masked.mrcs'
				newline=' '.join(words)+'\n'
				outputstar.write(newline)	
	outputstar.close()
	### join stacks together
	relion_preprocess_run='relion_stack_create --i '+output_rootname+'_tmp.star  --o '+output_rootname + '  --one_by_one'
	run_command(relion_preprocess_run)


	run_command('relion_reconstruct --i '+output_rootname+'.star '+'--o '+output_rootname+'.mrc '+' --ctf')
		
if __name__=="__main__":
	main()