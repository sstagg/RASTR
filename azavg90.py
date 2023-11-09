#!/usr/bin/env python

### usage: ./azavg.py inputfile

import numpy
import mrcfile
import sys
import argparse
### mask means do averaging along z axis within the mask.
### keep means keep the voxels outside the mask
def parse(args):
	parser = argparse.ArgumentParser(description='full box average or within mask')
	parser.add_argument('-m','--mask',action='store',default=None,dest='mask')
	parser.add_argument('-k','--keep',action='store_true',default=False,dest='keep')
	return parser.parse_args(args)


invol = sys.argv[1]
outvol = invol.split('.')[0]+'azavg.mrc'
parse = parse(sys.argv[2:])


vol = mrcfile.read(invol)
vol = vol.astype(numpy.float32)
if parse.mask:
	maskvol = mrcfile.read(parse.mask)
	if maskvol.shape != vol.shape:
		sys.exit('mask shape not consistent with volume')
	elif maskvol.any()>1 or maskvol.any()<0:
		sys.exit('mask data not binary')
	else:
		nvol = vol*maskvol
		average = nvol.mean(axis=0)
		maskavg = maskvol.mean(axis=0)
		maskavg[maskavg<1] = 1
		average = average/maskavg
		nvol[:] = average
		nvol = nvol*maskvol
		### if voxel outside mask need to be kept
		if parse.keep:
			reverse_mask = 1-mask
			nvol = nvol+(vol*reverse_mask)
		nvol = nvol.astype(numpy.float32)
		mrcfile.write(outvol, data=nvol)		
else:
	zavg = numpy.average(vol[int(0.1*vol.shape[0]):int(0.9*vol.shape[0])][:][:], axis=0)
	cylavg = numpy.zeros(vol.shape)
	for n in range(vol.shape[0]):
		cylavg[n] = zavg
	cylavg = cylavg.astype(numpy.float32)
	mrcfile.write(outvol, data=cylavg)

