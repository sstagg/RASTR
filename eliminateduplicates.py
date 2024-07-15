#!/usr/bin/env python

import sys
from pyami import correlator, peakfinder, convolver, mrc, imagefun
from appionlib.apImage import imagenorm
import numpy
import glob
import scipy
from scipy.stats import pearsonr
from appionlib import starFile
import optparse

def getShift(ref,image,corr='cc', test=False, kernel=1):
	if corr=='cc':
		cc=correlator.cross_correlate(ref, image)
	elif corr=='pc':
		cc=correlator.phase_correlate(ref, image, zero=False)
		conv=convolver.Convolver()
		kernel=convolver.gaussian_kernel(kernel)
		conv.setKernel(kernel)
		cc=conv.convolve(image=cc)		
	if test is True:
		mrc.write(cc, 'cc.mrc')
	peak=peakfinder.findSubpixelPeak(cc)
	shift=correlator.wrap_coord(peak['subpixel peak'],ref.shape)
	#returns y,x
	#print shift
	return shift

def shiftParticle(image,shift, test=False):
	shifted=scipy.ndimage.shift(image, shift, mode='wrap', cval=0.0)
	if test is True:
		mrc.write(shifted,'shifted.mrc')
	return shifted

def writeStar(starname,opticdata,particledata):
	outstar=starFile.StarFile(starname)
	
	#optics block
	opticlabelsout, opticdataout=starDataOut(opticdata)
	outstar.buildLoopFile(opticname,opticlabelsout,opticdataout)
	
	#particles block
	particlelabelsout,particledataout=starDataOut(particledata)
	outstar.buildLoopFile(particlesname,particlelabelsout,particledataout)
	
	outstar.write()

def starDataOut(blockdata):
	datalst=[]
	labels=list(blockdata[0].keys())
	for line in blockdata:
		params = ""
		for p in list(line.values()):
			params += " %s" % p
		datalst.append(params)
	return labels, datalst


def parseOptions():
	parser=optparse.OptionParser()
	parser.add_option('--star', dest='star', type='str', help='input star file')
	parser.add_option('--mrc', dest='mrc', type='str', help='input mrc file')
	parser.add_option('--outprefix', dest='outprefix', type='str', help='output prefix')
	parser.add_option('--thresh', dest='thresh', default=0.2, type='float', help='threshold for picking')
	parser.add_option('--searchdepth', dest='searchdepth', default=25, type='int', help='threshold for picking')

	options, args=parser.parse_args()

	if len(args) != 0 or len(sys.argv) == 1:
		parser.print_help()
		sys.exit()
	
	return options


if __name__=='__main__':
	
	options=parseOptions()

	#parse input star file
	star=starFile.StarFile(options.star)
	star.read()
	opticname='data_optics'
	opticblock=star.getDataBlock(opticname)
	opticlabels=opticblock.getLabelDict()
	opticdata=opticblock.getLoopDict()
	
	particlesname='data_particles'
	particlesblock=star.getDataBlock(particlesname)
	particlelabels=particlesblock.getLabelDict()
	particledata=particlesblock.getLoopDict()
	
	totalparticles=len(particledata)




	duplicates=[]
#	for n in range(totalparticles-options.searchdepth):
	newstar=[]
	for n in range(totalparticles):
		if n in duplicates:
			print ("skipping",n)
			continue

		newstar.append(particledata[n])			
		a=mrc.read(options.mrc,zslice=n)
		end=n+1+options.searchdepth
		if end > totalparticles:
			end=totalparticles
		for p in range(n+1,end):
			b=mrc.read(options.mrc,zslice=p)

			shift=getShift(a, b , corr='cc')
			aligned=shiftParticle(b,shift)
			comparea=imagefun.clipImage(a,5)
			#comparea=(comparea-comparea.mean())/comparea.std()
			comparea=comparea.flatten()

			compareb=imagefun.clipImage(aligned,5)
			#compareb=(compareb-compareb.mean())/compareb.std()
			compareb=compareb.flatten()

			cc=pearsonr(comparea,compareb)
			print (n,p,cc)
			if cc[0] > options.thresh:
				print("skipping",n,p,cc)
				duplicates.append(p)

	print (duplicates)

	#newstar=[]
	#for n in range(25):
	#	if n in duplicates:
	#		continue
	#	newstar.append(particledata[n])
	
	writeStar(options.outprefix+'.star', opticdata, newstar)
			
