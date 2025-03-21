#! /usr/bin/env python
#### script for RASTR by Ruizi
#### this script takes a pre-made mask. 
#### a typical run command: ./RASTR_02.py --star_in run_it005_data.star  --model run_it005_class001azavg.mrc  --angpix 2.02  -k -o rootname -ma spheremask.mrc  -al 0,90,180,270  --pad 3

import csv
import numpy as np
import sys
import math
import os, errno
import subprocess
import argparse
import time
import logging
import copy
from pyami import mrc
from pyami import convolver
from scipy.ndimage import rotate,shift,gaussian_filter
from multiprocessing import Process
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
	
	process=subprocess.Popen(command.split(), stdout=subprocess.PIPE)
	output, error = process.communicate()

#flips the values in an array, only values between 0 and 1
def flip(maskvolume):
	if maskvolume.max()<=1 and maskvolume.min()>=0:
		return 1-maskvolume
	else:
		print ('not binary mask')
	
### using scipy to rotate a 3d matrix with rot, tilt, psi angle, rotation order verified using cisTEM generate.
def volumerotation(volume,rot=0,tilt=0,psi=0,order=1):
	newvolume=copy.deepcopy(volume)
	if rot!=0:
		newvolume=rotate(newvolume,rot,axes=[1,2],order=order,reshape=False)
	if tilt!=0:	
		newvolume=rotate(newvolume,-tilt,axes=[0,2],order=order,reshape=False)
	if psi!=0:
		newvolume=rotate(newvolume,psi,axes=[1,2],order=order,reshape=False)
	return newvolume



### project 3d matrix to xy plane
def volumeprojection(volume,rot=0,tilt=0,psi=0,x=0,y=0,order=0):
	volume_rotated=volumerotation(volume,rot,tilt,0,order=order)
	slice_prj=np.sum(volume_rotated,axis=0)
	slice_prj=rotate(slice_prj,psi,axes=(0,1),reshape=False,order=order)
	slice_prj=shift(slice_prj,(-y,-x),mode='wrap')
	return slice_prj
### projection of a mask, max 1, min 0

#### optimization of maskprojection function is not finished!!!!!!!!!!!!!!
def mask_projection(volume,rot=0,tilt=0,psi=0,x=0,y=0,order=0,sigma=0):
	slice_mask=volumeprojection(volume,rot,tilt,psi,x,y,order)
	slice_mask[slice_mask<8]=0.0
	slice_mask[slice_mask>=8]=1.0
	#slice_mask=gaussian_filter(slice_mask,5)
	#slice_mask[slice_mask<0]=0.0
	#slice_mask[slice_mask>1]=1.0
	return slice_mask

### parse the top part of star file, return a dictionary with different column and number. eg paracolumn['rot']=8, then the 9th value in a row represent rot angle.
def column(filename):
	if filename[-4:]=='star':
		paracolumn={}
		fobj=open(filename,'r')
		lines=fobj.readlines()
		fobj.close()
		datastart=False
		for line in lines:
			words=line.split()
			if line=='data_particles\n' or line=='data_\n':
				datastart=True
			if len(words)==2 and datastart: 
				if words[0]=='_rlnMicrographName':
					paracolumn['mic']=int(words[1][1:])-1
				if words[0]=='_rlnImageName':
					paracolumn['image']=int(words[1][1:])-1
				if words[0]=='_rlnMagnification':
					paracolumn['mag']=int(words[1][1:])-1
				if words[0]=='_rlnAngleRot':
					paracolumn['rot']=int(words[1][1:])-1
				if words[0]=='_rlnAngleTilt':
					paracolumn['tilt']=int(words[1][1:])-1
				if words[0]=='_rlnAnglePsi':
					paracolumn['psi']=int(words[1][1:])-1
				if words[0]=='_rlnOriginX':
					paracolumn['shx']=int(words[1][1:])-1
				if words[0]=='_rlnOriginY':
					paracolumn['shy']=int(words[1][1:])-1
### this is copied from other scirpt for modify meta data, par may not be needed here.
	elif filename[-3:]=='par':
		paracolumn={}
		paracolumn['mag']=6
		paracolumn['psi']=1
		paracolumn['tilt']=2
		paracolumn['rot']=3
		paracolumn['shx']=4
		paracolumn['shy']=5
	else:
		print ("unknown file type")
	return paracolumn

### funcition for subtration, utilizing relion_project --subtract_exp
### mrcs_subtracted is a global dictionary in main, angle is its key
def subtraction(tempangle,starfilename,pad):
	mrcs_subtracted[tempangle]=scratch+'/'+results.output_rootname+'_'+str(tempangle)+'_subtracted'
	relion_subtraction='relion_project --i '+subtractionmodels[tempangle]+'  --o '+mrcs_subtracted[tempangle]+'  --angpix '+angpix+' --ang '+starfilename+' --ctf --pad '+str(pad)+' --subtract_exp'
	logger2.info('angle '+str(tempangle)+': '+relion_subtraction)
	run_command(relion_subtraction)
	logger2.info('angle '+str(tempangle)+': done!')

### using numpy to bin 3d volume by 2
def bin_3d(volume):
	boxsize=volume.shape[0]
	if boxsize%2 !=0:
		sys.exit()
	newboxsize=boxsize/2
	volume=volume.reshape(boxsize/2,2,boxsize/2,2,boxsize/2,2)
	newvolume=np.mean(volume,axis=5)
	newvolume=np.mean(newvolume,axis=3)
	newvolume=np.mean(newvolume,axis=1)
	return newvolume

### using numpy to unbin 2d array. eg 3x3 to 6x6
def unbin_2d(singleslice):
	singleslice=np.repeat(singleslice,2,axis=0)
	singleslice=np.repeat(singleslice,2,axis=1)
	return singleslice
### logger function, inherited from Peter
def setupLogger():

	logging.basicConfig(level=logging.DEBUG,
						format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
						datefmt='%m-%d %H:%M',
						filename=output_rootname+'.log',
						filemode='w')
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
	
	logger1 = logging.getLogger('Making models              ')
	logger2 = logging.getLogger('Subtraction                        ')
	logger3 = logging.getLogger('Extraction                 ')
	logger4 = logging.getLogger('Creating final files       ')
	#output the initial values

	commandline=' '
	logging.info(commandline.join(sys.argv))
	logging.info('star in          = '+ results.star_in)
	logging.info('model to delete  = '+ results.model)
	logging.info('angpix           = '+ results.angpix)
	logging.info('output rootname  = '+ output_rootname)
	logging.info('keep scratch     = '+ str(results.keep_scratch))
	logging.info('scratch folder   = '+ scratch)	

	return logger1, logger2, logger3, logger4


def relionmask(angle):
	logger3.info('starting extraction for angle: '+ str(angle))
	logger3.info('projecting masks')
	maskproj_run='relion_project --i '+whitemasksbin[angle]+'  --o '+maskprojection[angle]+'  --angpix '+angpix+' --ang '+initial_star_in+' --pad 3'
	run_command(maskproj_run)
	logger3.info('done projection masks for angle: '+str(angle))

	mrcsfileobj=open(mrcs_subtracted[angle]+'.mrcs','rb')
	masksfileobj=open(maskprojection[angle]+'.mrcs','rb')
	
	maskheaderbytes=masksfileobj.read(1024)
	maskheaderdict=mrc.parseHeader(maskheaderbytes)

	newfile=mrcs_masked[angle]
	outputobj=open(newfile,'wb')
	### writing the same header
	outputobj.write(mrc.makeHeaderData(headerdict))

	f=open(maskprojection[angle]+'.star','r')
	lines=f.readlines()
	f.close()
	for line in lines:
		words=line.split()
		if len(words)>10:
			i=int(words[paracolumn['image']].split('@')[0])-1

			singleslice=mrc.readDataFromFile(mrcsfileobj,headerdict,zslice=i)
			maskslice=mrc.readDataFromFile(masksfileobj,maskheaderdict,zslice=i)

			maskslice[maskslice<=8]=0.0
			maskslice[maskslice>8]=1.0
			maskslice=unbin_2d(maskslice)
			maskslice=convolver.Convolver().convolve(maskslice,convolver.gaussian_kernel(5))

			maskedslice=maskslice*singleslice
			mrc.appendArray(maskedslice,outputobj)
	mrcsfileobj.close()
	masksfileobj.close()
	outputobj.close()
	logger3.info('done extraction for angle: '+str(angle))


### main function for masking mrcs.
def mask2d(angle):
	logger3.info('starting extraction for angle: '+str(angle))
	mrcsfileobj=open(mrcs_subtracted[angle]+'.mrcs','rb')
	newfile=mrcs_subtracted[angle]+'_masked.mrcs'
	outputobj=open(newfile,'wb')
	### writing the same header 
	outputobj.write(mrc.makeHeaderData(headerdict))
	
	f=open(initial_star_in,'r')
	lines=f.readlines()
	f.close()
	### reading different white mask for each angle, using the same star file orientations, this way, white mask match the mrcs
	maskvolume=mrc.read(whitemasks[angle])

	### bin volume to run faster
	binrepeat=0
	while maskvolume.shape[0]>150:
		maskvolume=bin_3d(maskvolume)
		binrepeat+=1
	scale=2**binrepeat
	### extraction one by one
	for line in lines:
		words=line.split()
		if len(words)>2:
			rot=float(words[paracolumn['rot']])
			tilt=float(words[paracolumn['tilt']])
			psi=float(words[paracolumn['psi']])
			### dividing by scale because volume is binned to reduce running time
			x=float(words[paracolumn['shx']])/scale
			y=float(words[paracolumn['shy']])/scale
			i=int(words[paracolumn['image']].split('@')[0])-1
			singleslice=mrc.readDataFromFile(mrcsfileobj,headerdict,zslice=i)
			maskproj=mask_projection(maskvolume,rot,tilt,psi,x,y,2,3)
			for j in range(binrepeat):
				maskproj=unbin_2d(maskproj)
			maskproj=convolver.Convolver().convolve(maskproj,convolver.gaussian_kernel(5))
			maskedslice=maskproj*singleslice
			mrc.appendArray(maskedslice,outputobj)
	outputobj.close()
	logger3.info('done extration for angle: '+str(angle))
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
	parser.add_argument('-c','--center', action='store_true', default=False,
						dest='center',
						help='      center the masked area in a new smaller box')
	parser.add_argument('-o', '--output', action='store', dest='output_rootname',
						 help='     output rootname for star/mrcs, default is RASTR_particles')
	
	parser.add_argument('-k', '--keep_scratch', action='store_true', default=False,
						dest='keep_scratch',
						help='      keep scratch folder, default is False, only final star, mrcs, and reconstruction output on default')
	parser.add_argument('-pad','--pad',action='store',dest='pad',default=2)
	parser.add_argument('-or','--order',action='store',dest='order',default=3)
	parser.add_argument('-f','--test', action='store_true', default=False,
						dest='test_TF',
						help='      option for testing, dont use, was needed during setup')
	parser.add_argument('-b','--both',action='store_true', default=False,
							dest='both',
						help='	do both masking and center. End with two stacks, one with section masked keeping same box size as original and the second with it centered in a smaller box')	

	parser.add_argument('-g','--gauss',type=int, default=5, dest='gauss', help='	Sigma for gaussian edge to add to sphere. Default = 5')	
	
	parser.add_argument('-al','--anglelist',type=str,action='store',default='0,90,180,270', dest='anglelist',
				help=' input list of angles, separate by comma')
	parser.add_argument('-ma','--mask',type=str,action='store',default=None,dest='mask',help='input mask for refinement, make the region to be kept 1')


	parser.add_argument('-sin','--single',action='store_true',default=False,dest='single')

	results = parser.parse_args()
	
	#make sure all the required arguments are here

	if results.star_in == None:
		print ('Required arguments: --star_in <star file> --model <model file> -- angpix <angpix>; use -h or --help to see all options')
		sys.exit()
	
	if results.model == None:
		print ('Required arguments: --star_in <star file> --model <model file> -- angpix <angpix>; use -h or --help to see all options')
		sys.exit()
	
	
	if results.angpix == None:
		print ('Required arguments: --star_in <star file> --model <model file> -- angpix <angpix>; use -h or --help to see all options')
		sys.exit()
	return (results)
		
#Main program
if __name__=="__main__":
	#Input arguments
	results=parseOptions()
	
	#get box size first because some defaults depend on it:

	input_model=mrc.read(results.model)
	box_size = int(input_model.shape[0])

	#set defaults
	pad=results.pad
	order=results.order
	
	if results.output_rootname == None:
			output_rootname = 'RASTR_particles'
	else:
			output_rootname = results.output_rootname

	### initial dictionaries of different files,key is angle, item is filename.. eg blackmasks['0.0']='blackmask0.0.mrc'
	angles=[]
	subtractionmodels={}
	blackmasks={}
	whitemasks={}
	whitemasksbin={}
	maskprojection={}
	mrcs_subtracted={}
	mrcs_masked={}

	#Default values to be printed in log file	

	initial_star_in = results.star_in
	angpix = str(results.angpix)
	

	#create a scratch directory

	scratch = 'scratch'
	scratch_num = 2
	if os.path.isdir(scratch) is True:
			while os.path.isdir(scratch) is True:
					scratch = 'scratch'+str(scratch_num)
					scratch_num+=1
	make_tmp(scratch)

	#create a log file
	logger1,logger2,logger3,logger4=setupLogger()

	#####             #####
	##### START RASTR #####
	#####             #####
	

	### parse a list of angles for subtraction
	angles=results.anglelist.split(',')
	for i in range(len(angles)):
		angles[i]=float(angles[i])

	whitemask=mrc.read(results.mask)


	### making subraction models
	for angle in angles:
		subtractionfile=scratch+'/'+'subtractionmodel'+str(angle)+'.mrc'
		b_maskfile=scratch+'/blackmask'+str(angle)+'.mrc'
		w_maskfile=scratch+'/whitemask'+str(angle)+'.mrc'
		subtractionmodels[angle]=subtractionfile
		whitemasks[angle]=w_maskfile
		blackmasks[angle]=b_maskfile
		whitemasksbin[angle]=scratch+'/whitemask'+str(angle)+'bin2.mrc'

		whitemaskvolume=volumerotation(whitemask,angle,0,0,order=5)
		whitemaskvolume=gaussian_filter(whitemaskvolume,sigma=3)
		whitemaskvolume[whitemaskvolume>1]=1.0
		whitemaskvolume[whitemaskvolume<0]=0.0
		blackmaskvolume=flip(whitemaskvolume)

		mrc.write(whitemaskvolume,w_maskfile)
		mrc.write(blackmaskvolume,b_maskfile)

		whitemaskvolume_bin2=bin_3d(whitemaskvolume)
		mrc.write(whitemaskvolume_bin2,whitemasksbin[angle])

		subtractionvolume=blackmaskvolume*input_model
		mrc.write(subtractionvolume,subtractionfile)
		logger1.info(subtractionfile)
	logging.info('Done making subtraction models')



	####subtraction
	if results.single:
		for angle in angles:
			subtraction(angle,initial_star_in,pad)
	else:
		processes=[]
		for angle in angles:
			p=Process(target=subtraction,args=(angle,initial_star_in,pad))
			p.start()
			processes.append(p)

		for p in processes:
			p.join()
	logging.info('Done with subtraction')

	for angle in angles:
		mrcs_subtracted[angle]=scratch+'/'+results.output_rootname+'_'+str(angle)+'_subtracted'


	#### mask target region
	####  speed of this step is mainly affected by order number. 
	#### It is within range 0-5, the larger the longer for scipy.ndimage.rotate to run
	logging.info('')
	logger3.info('Extracting particles')
	processes=[]
	paracolumn=column(initial_star_in)
	inputobj=open(mrcs_subtracted[angles[0]]+'.mrcs','rb')
	headerbytes=inputobj.read(1024)
	inputobj.close()
	headerdict=mrc.parseHeader(headerbytes)

	if results.single:
		for angle in angles:
			maskprojection[angle]=scratch+'/maskprojection_'+str(angle)
			mrcs_masked[angle]=mrcs_subtracted[angle]+'_masked.mrcs'
			relionmask(angle)
	else:
		for angle in angles:
			maskprojection[angle]=scratch+'/maskprojection_'+str(angle)
			mrcs_masked[angle]=mrcs_subtracted[angle]+'_masked.mrcs'
		for angle in angles:
			p=Process(target=relionmask,args=(angle,))
			p.start()
			processes.append(p)
		for p in processes:
			p.join()
	logging.info('Done with masking')
	

	### join star file together
	outputstar=open('tmp.star','w')
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
	relion_preprocess_run='relion_stack_create --i tmp.star  --o '+output_rootname + '  --one_by_one'
	run_command(relion_preprocess_run)


	run_command('relion_reconstruct --i '+output_rootname+'.star '+'--o '+output_rootname+'.mrc '+' --ctf')
		
