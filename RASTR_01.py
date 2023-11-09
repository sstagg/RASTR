#!/usr/bin/env python

import csv
import numpy as np
import sys
import math
import os, errno
import subprocess
import argparse
import fileinput
import shutil
import logging
from pyami import imagefun
from pyami import convolver
from pyami import mrc


##Functions
#makes a sphere in a box
def make_sphere(box, radius, xyz_center_pos):
	vol=np.ones((box, box, box), dtype='float32')
	for x in range(box):
		for y in range(box):
			for z in range(box):
				dx=x-xyz_center_pos[0]
				dy=y-xyz_center_pos[1]
				dz=z-xyz_center_pos[2]
				d=math.sqrt(dx*dx+dy*dy+dz*dz)
				if d < radius:
					vol[z,y,x]=0
	return vol

#Makes a directory
def make_tmp(directory):
	try:
		os.makedirs(directory)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise



#phi rotation of box to get center coord takes IMOD coordinates
def z_rotate(coord_0, angle, box_0):
	c0=np.array([coord_0]).reshape(3,1)-(box_0/2)
	a = np.deg2rad(float(angle))
        rot_matrix=np.array([\
        np.cos(a), -np.sin(a), 0, \
        np.sin(a), np.cos(a), 0, \
        0, 0, 1]).reshape((3,3))
        c1=np.dot(rot_matrix, c0)
	c1=np.round(c1)
	return c1

#generating rotation matrix 
def full_rotate(phi, theta, psi1):
        r=np.deg2rad(float(phi))
        t=np.deg2rad(float(theta))
        p=np.deg2rad(float(psi1))
        rm=np.array([\
        np.cos(p)*np.cos(t)*np.cos(r)-np.sin(p)*np.sin(r), -np.cos(r)*np.sin(p)-np.cos(p)*np.cos(t)*np.sin(r), np.cos(p)*np.sin(t), \
        np.cos(p)*np.sin(r)+np.cos(t)*np.cos(r)*np.sin(p), np.cos(p)*np.cos(r)-np.cos(t)*np.sin(p)*np.sin(r), np.sin(p)*np.sin(t), \
        -np.cos(r)*np.sin(t), np.sin(t)*np.sin(r), np.cos(t)]).reshape((3, 3))
        return rm

#rounds up
def roundup(x):
        return int(math.ceil(x / 10.0)) * 10

#counts the number of models in a  directory - ended up not using
def models_count(directory):
	a = len([f for f in os.listdir(directory)
		if f.endswith('.mrc') and os.path.isfile(os.path.join(directory, f))])
	return a

#makes a list of all the files in a directory with the extension extens (given as a string, no period), ended up not using
def files_list(directory, extens):
	files_in_dir = []
	for root, dirs, files in os.walk(directory):
		for file in files:
			if file.endswith('.'+ extens ):
				files_in_dir.append(file)
	return files_in_dir
#Makes a list of files in a directory that end with _extens, ended up not using

def files_list_dash(directory, extens):
        files_in_dir = []
        for root, dirs, files in os.walk(directory):
                for file in files:
                        if file.endswith('_'+ extens ):
                                files_in_dir.append(file)
        return files_in_dir


#takes in a star file and outputs a list with ((column name, column number -1))
def parse_star(file):
        ssvin = open(file, 'rb')
        reader=csv.reader(ssvin, delimiter=' ')
        number_of_arguments = 0
        argument_list = []
        for row in reader:
                row = filter(None, row)
                if len(row) == 2:
                        number_of_arguments += 1
                        num = row[1]
                        num =int(num.replace("#",""))
                        argument_list.append((row[0], num-1))
        ssvin.close()
        return argument_list


#makes a zero box in a larger box of ones, coordinates from 3dmod (basically corner as origin)
#should be 'box mask=np.ones' and 'start[2]+k]=0' if opposite is because testing
def make_box(box, x_y_z_center_point, cut_out):
	if results.test_TF is True:
		x_y_z_starting_point = x_y_z_center_point-(cut_out/2)
		box_mask=np.zeros((box,box,box),dtype='float32')
		start=np.array((x_y_z_starting_point[2],x_y_z_starting_point[1],x_y_z_starting_point[0]),dtype='int')
		size=np.array((cut_out,cut_out,cut_out),dtype='int')
		for i in range(size[0]+1):
			for j in range(size[1]+1):
				for k in range(size[2]+1):
					box_mask[start[0]+i, start[1]+j, start[2]+k]=1
	else:
		x_y_z_starting_point = x_y_z_center_point-(cut_out/2)
		box_mask=np.ones((box,box,box),dtype='float32')
		start=np.array((x_y_z_starting_point[2],x_y_z_starting_point[1],x_y_z_starting_point[0]),dtype='int')
		size=np.array((cut_out,cut_out,cut_out),dtype='int')
		for i in range(size[0]+1):
			for j in range(size[1]+1):
				for k in range(size[2]+1):
					box_mask[start[0]+i, start[1]+j, start[2]+k]=0

	return box_mask

#run a command, argument is a string
def run_command(command):
	print command
	process=subprocess.Popen(command.split(), stdout=subprocess.PIPE)
	output, error = process.communicate()

##Like make box but masks and takes coord with origin in center
def mask_box(mrcs_i, mrcs_o, box, x_y_z_center_point, cut_out, x_shift, y_shift, particle_number):
	x_y_z_starting_point=x_y_z_center_point+box/2-(cut_out/2)
	start=np.array((int(round(float(x_y_z_starting_point[1])-y_shift)), int(round(float(x_y_z_starting_point[0])-x_shift))), dtype='int')
	size=np.array((int(cut_out), int(cut_out)), dtype='int')
	for i in range(size[0]+1):
		if start[0]+i <= (box-1):
			for j in range(size[1]+1):
				if start[1]+j <= (box-1):
					mrcs_o[particle_number, start[0]+i, start[1]+j]=mrcs_i[particle_number,start[0]+i,start[1]+j]


#Like mask_box, but instead masks all but a sphere	
def mask_circle(mrcs_i, mrcs_o, box, x_y_z_center_point, radius, x_shift, y_shift, particle_number, mean):
	xi = int(x_y_z_center_point[0]-x_shift+box/2)
	yi = int(x_y_z_center_point[1]-y_shift+box/2)	
	for x in range(box):
		for y in range(box):
			dx=x-xi
			dy=y-yi
			d=math.sqrt(dx*dx+dy*dy)
			if d < radius:
				mrcs_o[particle_number, y, x]=mrcs_i[particle_number, y, x]
			else:
				mrcs_o[particle_number, y, x]=mean


#instead of masking a square, this extracts a square and adds to a stack with the same size 
def extract_square(mrcs_i, mrcs_o, box, x_y_z_center_point, radius, x_shift, y_shift, particle_number, mean):
	xi = int(x_y_z_center_point[0]+box/2-radius-x_shift)
	yi = int(x_y_z_center_point[1]+box/2-radius-y_shift)
	for x in range(2*radius):
		if xi+x <= box-1:
			for y in range(2*radius):
				if yi+y <= box-1:
					mrcs_o[particle_number, y, x]=mrcs_i[particle_number, yi+y, xi+x]
				else:
					mrcs_o[particle_number, y, x]=mean
		else:
			mrcs_o[particle_number, y, x]=mean
			

#Makes a cylinder in a box 
def make_cylinder(box, radius):
        vol=np.ones((box, box, box), dtype='float32')
        for z in range(box):
                for x in range(box):
                        for y in range(box):
                                dx=x-(box/2)
                                dy=y-(box/2)
                                d=math.sqrt(dx*dx+dy*dy)
                                if d < radius:
                                        vol[z,y,x]=0
        return vol

#flips the values in an array, only values between 0 and 1
def flip(mrcs_in):
        for x in range(mrcs_in.shape[2]):
                for y in range(mrcs_in.shape[1]):
                        for z in range(mrcs_in.shape[0]):
                                mrcs_in[z,y,x]=1-mrcs_in[z,y,x]


def mask_circle2(mrcs_i, box, x_y_z_center_point, radius, x_shift, y_shift, particle_number, gauss_sigma=20):
	mask1 = imagefun.filled_circle((box, box), radius, center=(-x_y_z_center_point[1]-y_shift+float(box)/2, x_y_z_center_point[0]-x_shift+float(box)/2))
	mask1 = 1 - mask1
	con = convolver.gaussian_kernel(gauss_sigma)
	conv = convolver.Convolver()
	cmask = conv.convolve(mask1, con)
	mrcs_i[particle_number] = np.multiply(mrcs_i[particle_number], cmask)


def make_sphere2(box, radius, xys_center_pos, gauss):
	mask1 = imagefun.filled_sphere((box,box,box), center=(xyz_center_pos[2], xyz_center_pos[1], xyz_center_pos[0]))
	mask1 = 1 - mask1
	con = convolver.gaussian_kernel(gauss)
	conv = convolver.Convolver()
	cmask = conv.convolve(mask1, con)
	return cmask


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
	logging.info('x_start          = '+ str(x_start))
	logging.info('radius           = '+ str(radius))
	logging.info('n boxes to cut   = '+ str(n_spheres))
	logging.info('tube radius      = '+ str(results.trad))
	logging.info('center           = '+ str(results.center))
	logging.info('gauss		= '+ str(results.gauss))
	logging.info('output rootname  = '+ output_rootname)
	logging.info('keep scratch     = '+ str(results.keep_scratch))
	logging.info('test             = '+ str(results.test_TF))
	logging.info('scratch folder   = '+ scratch)	
	logging.info('both             = '+ str(results.both))
	return logger1, logger2, logger3, logger4	

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
	
	
	parser.add_argument('-f','--test', action='store_true', default=False,
	                    dest='test_TF',
	                    help='      option for testing, dont use, was needed during setup')
	parser.add_argument('-b','--both',action='store_true', default=False,
                            dest='both',
	                    help='	do both masking and center. End with two stacks, one with section masked keeping same box size as original and the second with it centered in a smaller box')	

	parser.add_argument('-g','--gauss',type=int, default=5, dest='gauss', help='	Sigma for gaussian edge to add to sphere. Default = 5')	
	
	results = parser.parse_args()
	
	#make sure all the required arguments are here

	if results.star_in == None:
	        print 'Required arguments: --star_in <star file> --model <model file> -- angpix <angpix>; use -h or --help to see all options'
	        sys.exit()
	
	if results.model == None:
	        print 'Required arguments: --star_in <star file> --model <model file> -- angpix <angpix>; use -h or --help to see all options'
	        sys.exit()
	
	
	if results.angpix == None:
	        print 'Required arguments: --star_in <star file> --model <model file> -- angpix <angpix>; use -h or --help to see all options'
	        sys.exit()
	return (results)
		
#Main program
if __name__=="__main__":

	#Input arguments
	results=parseOptions()
	
	#get box size first because some defaults depend on it:

	input_mrc=mrc.read(results.model)
	box_size = int(input_mrc.shape[0])

	#set defaults
	
	if results.x_start == None:
	        x_start = int(round((float(box_size)*3)/4))
	else:
	        x_start = results.x_start
	
	if results.n_spheres == None:
	        n_spheres = 9
	else:
	        n_spheres = results.n_spheres
	
	if results.radius == None:
	        radius = int(round((float(box_size)*3)/16))
	else:
	        radius = results.radius
	if results.output_rootname == None:
	        output_rootname = 'RASTR_particles'
	else:
	        output_rootname = results.output_rootname


	#Default values to be printed in log file	

        initial_star_in = results.star_in
        angpix = str(results.angpix)
        angle_to_add = float(360)/float(n_spheres)
	
	#create a scratch directory

	scratch = 'scratch'
	scratch_num = 2
	if os.path.isdir(scratch) is True:
        	while os.path.isdir(scratch) is True:
                	if scratch_num == 2:
                        	scratch = scratch+str(scratch_num)
                	else:
                        	scratch = scratch.replace((str(int(scratch_num)-1)),str(scratch_num))
                	scratch_num+=1


        make_tmp(scratch)

	#create a log file
	logger1,logger2,logger3,logger4=setupLogger()

	#####             #####
	##### START RASTR #####
	#####             #####



	#FIND CENTER OF AREA TO EXTRACT/MASK
	xyz_center = [int(x_start), box_size/2, box_size/2]
	#CREATE A MODELS DICTIONARY TO STORE VALUES OF EACH MODEL FOR LATER USE
	models_a = {}
	#START AT THE SECTION CLOSEST TO EDGE
	angle = 0
	#CREATING MODELS FOR SUBTRACTION
#	mask_a=make_sphere(box_size, radius, xyz_center)
	##CREATE MODEL IF LEAVING TUBE
#	if results.trad == None:
#		temp_sphere=scratch+'/gauss_sphere_p000.mrc'
#		mrc.write(mask_a, temp_sphere)
#		low_pass_command='proc3d '+temp_sphere+' '+temp_sphere+' apix=1 lp='+str(results.gauss)
#        	run_command(low_pass_command)
#		mask_b=mrc.read(temp_sphere)
#	##CREATE MODEL IF REMOVING TUBE
#	else:
#		mask_a=make_sphere(box_size, radius, xyz_center)
#		cyl=make_cylinder(box_size, results.trad)
#		flip(mask_a)
#		bite=np.multiply(cyl, mask_a)
#		mrc.write(bite, scratch+'/temp_mask_a.mrc')
#		low_pass_command='proc3d '+scratch+'/temp_mask_a.mrc '+scratch+'/temp_mask_b.mrc apix=1 lp='+str(results.gauss)
#		run_command(low_pass_command)
#		mask_b=mrc.read(scratch+'/temp_mask_b.mrc')
#		flip(mask_b)
	##MASKING ON THE INPUT MODEL
#        final_array=np.multiply(input_mrc, mask_b)

        ##WRITE OUT MODEL TO DISC
#        mrc_to_write=scratch+'/z_0_phi_000_.mrc'
#	if results.test_TF == True:
#		flip(final_array)
#	mrc.write(final_array,mrc_to_write)
#	logger1.info(mrc_to_write) 
	##ROTATE MODEL AROUND Z TO CREATE ADDITIONAL SECTIONS
	###ADD INITIAL MODEL TO DICTIONARY
#	models_a[mrc_to_write] = xyz_center[0:4], angle
	###FIND ANGLE TO STEP
#	angle=angle+angle_to_add
	###STEP AROUND Z AXIS USING PROC3D
        while abs(angle) < 360:
                if len(str(int(abs(angle)))) == 3:
			mrc2=scratch+'/z_0_phi_'+ str(int(abs(angle)))+'_.mrc'
                elif len(str(int(abs(angle)))) == 2:
                        mrc2=scratch+'/z_0_phi_0'+ str(int(abs(angle)))+'_.mrc'
		elif len(str(int(abs(angle)))) == 1:
                        mrc2=scratch+'/z_0_phi_00'+ str(int(abs(angle)))+'_.mrc'
		nc_array = z_rotate(xyz_center, angle, box_size)
		new_center = [float(nc_array[0])+(box_size/2), float(nc_array[1])+(box_size/2), float(nc_array[2])+box_size/2]
		mask_a = make_sphere(box_size, radius, new_center)
		temp_sphere=scratch+'/gauss_sphere_'+str(int(abs(angle)))+'.mrc'
		mrc.write(mask_a, temp_sphere)
		low_pass_command='proc3d '+temp_sphere+' '+temp_sphere+' apix=1 lp='+str(results.gauss)
		run_command(low_pass_command)
		mask_b=mrc.read(temp_sphere)
		final_array= np.multiply(input_mrc, mask_b)
		mrc.write(final_array, mrc2)
		logger1.info(mrc2)
		models_a[mrc2]=xyz_center[0:4], angle
		angle=angle+angle_to_add
	logging.info('Done making models to subtract!')
	logging.info('')
	#SUBTRACTION OF PROJECTED MODEL FROM INITIAL PARTICLES

	logger2.info('Starting relion subtractions')	

        
	model_number = int(0)
        particles_a={}
	for model in models_a:
		mrcs_out=model.split('.')[0]+'particles'
		particles_a[mrcs_out]=models_a[model]
		#when done add in the --ctf and subtract command
		if results.test_TF == True:
			relion_subtract_run = 'relion_project --i '+model+' --o '+mrcs_out+' --angpix '+angpix+' --ang '+initial_star_in
		else:
			relion_subtract_run = 'relion_project --i '+model+' --o '+mrcs_out+' --angpix '+angpix+' --ang '+initial_star_in+' --ctf --subtract_exp'
		logger2.info('Job '+str(model_number+1)+': '+relion_subtract_run)
		run_command(relion_subtract_run)
		logger2.info('Job '+str(model_number+1)+' done!')
		model_number+=1
	logging.info('Done with subtractions')


	
	##MASK OR EXTRACT PREVIOUSLY MASKED AREAS
	 
	logging.info('')
	logger3.info('Extracting particles')
	stack_number=0
	if results.both == True or results.center == True:
		final_list_cent = []
	if results.both == True or results.center == False:
		final_list = []
	for particles in particles_a:
		rn_out=particles.split('.')[0]+'_masked'
		if results.both == True or results.center == True:
			star_out_cent = rn_out+'_cent.star'
			mrcs_out_cent = rn_out+'_cent.mrcs'
			tsvout_2 = open(star_out_cent, 'wb')
			writer2=csv.writer(tsvout_2, delimiter=' ')
		if results.both == True or results.center == False:
			star_out=rn_out+'.star'
			mrcs_out=rn_out+'.mrcs'
			tsvout = open(star_out, 'wb')
			writer = csv.writer(tsvout, delimiter=' ')
		star_in=particles+'.star'
		mrcs_in=particles+'.mrcs'
		tsvin=open(star_in, 'rb')
		reader=csv.reader(tsvin, delimiter=' ')
		mrcs_to_read=mrc.read(mrcs_in)
		
		#VALUES NEEDED FOR MASKING/EXTRACTION

		mean_list=np.mean(mrcs_to_read, axis=(1,2))
		particle_number=int(mrcs_to_read.shape[0])
		center_circle=z_rotate(particles_a[particles][0],360-particles_a[particles][1], box_size)
	
		

		#PARSE THE STAR FILE FOR NEEDED VALUES
		star_arguments=parse_star(star_in)
		for col in star_arguments:
			if col[0] == '_rlnOriginX':
				OX = col[1]
		for col in star_arguments:
			if col[0] == '_rlnOriginY':
				OY = col[1]
		for col in star_arguments:
			if col[0] == '_rlnImageName':
				IName = col[1]
                for col in star_arguments:
                        if col[0] == '_rlnAngleRot':
                                AROT = col[1]
                for col in star_arguments:
                        if col[0] == '_rlnAngleTilt':
                                ATILT = col[1]
                for col in star_arguments:
                        if col[0] == '_rlnAnglePsi':
                                APSI = col[1]

		row_number=0
		###IF EXTRACTING TO A SMALLER BOX SIZE
		if results.center == True or results.both == True:
			temp_mrcs_array_cent=np.zeros((mrcs_to_read.shape[0], 2*radius, 2*radius), dtype='float32')
		###IF LEAVING BOX SIZE THE SAME
		if results.center == False or results.both == True:	
			temp_mrcs_array=np.zeros((mrcs_to_read.shape),dtype='float32')
		
		logger3.info('Stack '+str(stack_number+1)+' of '+str(len(particles_a)))
		for row in reader:
			#REMOVING EMPTY COLUMNS
			row = filter(None, row)
			#SELECTING ROWS WITH INFORMATION (NON-HEADER ROWS)
			if len(row) == len(star_arguments):
			#GENERATE ROTATION MATRIX
				rot_matrix=full_rotate(float(row[AROT]), float(row[ATILT]), float(row[APSI]))
				###CENTER OF MASKED SECTION
				end_position=np.dot(rot_matrix, center_circle)
				### IF EXTRACTING TO SMALLER BOX
				if results.center == True or results.both == True:
                                        extract_square(mrcs_to_read, temp_mrcs_array_cent, box_size, end_position, radius, float(row[OX]), float(row[OY]), row_number, mean_list[row_number])
				if results.center == False or results.both == True:
					mask_circle2(mrcs_to_read, box_size, end_position, radius, float(row[OX]), float(row[OY]), row_number, gauss_sigma=results.gauss)
					mask_rowOY = float(row[OY])-float(center_circle[2])
				### IF LEAVING IN BOX
			#CREATE NEW STAR FILE WITH THE PROPER SHIFTS FOR ALL THE PARTICLES TO ALIGN PROPERLY
				row[IName]= str(row_number+1)+'@'+mrcs_out
				row[AROT] = float(row[AROT]) - float(particles_a[particles][1])
				if results.center == False or results.both == True:
					row[OY] = mask_rowOY
					writer.writerow(row)
				if results.center == True or results.both == True:
					row[IName]=str(row_number+1)+'@'+mrcs_out_cent
					row[OY] = 0
					row[OX] = 0	
					writer2.writerow(row)
#				logger4.info('Particle '+str(row_number+1)+' of '+str(particle_number))
				row_number+=1
			#WRITE HEADER FOR STAR FILE
			else:
				if results.center == False or results.both == True:
					writer.writerow(row)
				if results.center == True or results.both == True:
					writer2.writerow(row)
		if results.center == False or results.both == True:
			mrc.write(mrcs_to_read, mrcs_out)
			normalize_cmd = 'e2proc2d.py '+mrcs_out+' '+rn_out+'_norm.mrcs --process=normalize'
			logger3.info('Normalizing new stack')
			logger3.info(normalize_cmd)
			run_command(normalize_cmd)
			logger3.info('Finished normalization, replacing old stack')
			mv_cmd = 'mv '+rn_out+'_norm.mrcs '+mrcs_out
			os.system(mv_cmd)
			logger3.info('Complete')
			tsvout.close()
			final_list.append(rn_out)
		if results.center == True or results.both == True:
			mrc.write(temp_mrcs_array_cent, mrcs_out_cent)
			normalize_cmd = 'e2proc2d.py '+mrcs_out_cent+' '+rn_out+'_cent_norm.mrcs --process=normalize'
			logger3.info('Normalizing new stack')
			logger3.info(normalize_cmd)
			run_command(normalize_cmd)
			logger3.info('Finished normalization, replacing old stack')
			mv_cmd = 'mv '+rn_out+'_cent_norm.mrcs '+mrcs_out_cent
			os.system(mv_cmd)
			logger3.info('Complete')
			tsvout_2.close()
			final_list_cent.append(rn_out+'_cent')
		##NORMALIZE THE STACK JUST CREATED (LESS LOAD THEN DOING ON THE LARGE STACK AT END
		tsvin.close()
		stack_number+=1
	
	## ZYZ Rotation matrix based of Euler angles 3 for initial Z(rot), 2 for Y(tilt), and 1 for final Z(psi), courtesy of wikipedia
	#Z(1)Y(2)Z(3)=	[	c1c2c3-s1s3,	-c3s1-c1c2s3,	c1s2	]
	#		[	c1s3+c2c3s1,	c1c3-c2s1s3,	s1s2	]
	#		[	-c3s2,		s2s3,		c2	]
	#	else:
	#		writer.writerow(row)


	
	#CREATE ONE STACK/STAR FILE FOR ALL
	logger4.info('Making single output star/stack')
	logger4.info('To combine:')
	if results.both == True or results.center == False:
		tmp_star_out=scratch+'/tmp_particles.star'
		starout=open(tmp_star_out, 'wb')
		writer3=csv.writer(starout, delimiter=' ')
		star_number = 0
		for rn in final_list:
			logger4.info(rn)		
	
			star=rn+'.star'
			star_columns=parse_star(star)
			readFile = open(star)
			lines = readFile.readlines()
			readFile.close()
			w = open(star, 'w')
			#WAS GETTING GAP IN PARTICLES BECAUSE OF EXTRA BLANK LINE, THIS REMOVES THAT
			w.writelines([item for item in lines[:-1]])
			w.close()
			starin=open(star, 'rb')
			reader2=csv.reader(starin, delimiter=' ')
			##ONLY 1 HEADER IN STAR FILE
			if star_number == 0:
				for row in reader2:
					row = filter(None, row)
					writer3.writerow(row)
			else:
				for row in reader2:
					row = filter(None, row)
					if len(row) == len(star_columns):
						writer3.writerow(row)
			starin.close()
			star_number+=1
			logger4.info('Stack '+str(star_number)+' of '+str(len(final_list)))
		starout.close()

		relion_preprocess_run='relion_preprocess --operate_on '+scratch+'/tmp_particles.star --operate_out '+output_rootname
		logger4.info('Running relion_preprocess job: '+relion_preprocess_run)
		run_command(relion_preprocess_run)
		os.system('mv '+output_rootname+'.mrcs.mrcs '+output_rootname+'.mrcs')
		logging.info('Preprocess finished, final outputs written to:')
		logging.info('star:    '+output_rootname+'.star')
		logging.info('stack:   '+output_rootname+'.mrcs')
		logging.info('Creating reconstruction of masked particles:')
		#add back ctf and maxres?? maybe for maxres, after test
		##RELION RECONSTRUCT TO USE AS MODEL OR JUST TO CHECK HOW WELL WORKED
		if results.test_TF == True:
		        relion_reconstruct_command='relion_reconstruct --i '+output_rootname+'.star --angpix '+angpix+' --o '+output_rootname+'.mrc'
		else:
		        relion_reconstruct_command='relion_reconstruct --i '+output_rootname+'.star --angpix '+angpix+' --o '+output_rootname+'.mrc --ctf'
		run_command(relion_reconstruct_command)
		logging.info('Created output model: '+output_rootname+'.mrc')


	if results.both == True or results.center == True:
                tmp_star_out_cent=scratch+'/tmp_particles_cent.star'
                starout=open(tmp_star_out_cent, 'wb')
                writer4=csv.writer(starout, delimiter=' ')
                star_number = 0
                for rn in final_list_cent:
                        logger4.info(rn)

                        star=rn+'.star'
                        star_columns=parse_star(star)
                        readFile = open(star)
                        lines = readFile.readlines()
                        readFile.close()
                        w = open(star, 'w')
                        #WAS GETTING GAP IN PARTICLES BECAUSE OF EXTRA BLANK LINE, THIS REMOVES THAT
                        w.writelines([item for item in lines[:-1]])
                        w.close()
                        starin=open(star, 'rb')
                        reader3=csv.reader(starin, delimiter=' ')
                        ##ONLY 1 HEADER IN STAR FILE
                        if star_number == 0:
                                for row in reader3:
                                        row = filter(None, row)
                                        writer4.writerow(row)
                        else:
                                for row in reader3:
                                        row = filter(None, row)
                                        if len(row) == len(star_columns):
                                                writer4.writerow(row)
                        starin.close()
                        star_number+=1
                        logger4.info('Stack '+str(star_number)+' of '+str(len(final_list_cent)))
                starout.close()

		##USE PREPROCESS TO COMBINE STACKS
		relion_preprocess_run='relion_preprocess --operate_on '+scratch+'/tmp_particles_cent.star --operate_out '+output_rootname+'_cent'
		logger4.info('Running relion_preprocess job: '+relion_preprocess_run)
		run_command(relion_preprocess_run)
		os.system('mv '+output_rootname+'_cent.mrcs.mrcs '+output_rootname+'_cent.mrcs')
		logging.info('Preprocess finished, final outputs written to:')
		logging.info('star:    '+output_rootname+'_cent.star')
		logging.info('stack:   '+output_rootname+'_cent.mrcs')
		logging.info('Creating reconstruction of masked particles:')
		#add back ctf and maxres?? maybe for maxres, after test
		##RELION RECONSTRUCT TO USE AS MODEL OR JUST TO CHECK HOW WELL WORKED
		if results.test_TF == True:
			relion_reconstruct_command='relion_reconstruct --i '+output_rootname+'_cent.star --angpix '+angpix+' --o '+output_rootname+'_cent.mrc'
		else:
			relion_reconstruct_command='relion_reconstruct --i '+output_rootname+'_cent.star --angpix '+angpix+' --o '+output_rootname+'_cent.mrc --ctf'
		run_command(relion_reconstruct_command)
		logging.info('Created output model: '+output_rootname+'_cent.mrc')
	###REMOVE SCRATCH IF STATED
	if results.keep_scratch is False:
		shutil.rmtree(scratch)
