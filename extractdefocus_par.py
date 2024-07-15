#! /usr/bin/env python
###extract defocus value to a txt file for azimuthal_average program
### usage ./extractdefocus_par.py parfile&starfile outputfile
import sys

inputfile=sys.argv[1]
outputfile=sys.argv[2]

fileobj=open(inputfile,'r')
lines=fileobj.readlines()
fileobj.close()

newlines=[]
#newlines.append('C  DF1  DF2\n')



def get_from_par():
	count=1
	for line in lines:
		words=line.split()
		if words[0]=='C':
			continue
		
		newlines.append(words[8]+' '+words[9]+' '+words[10]+'\n')
	
	fileobj=open(outputfile,'w')
	fileobj.writelines(newlines)
	fileobj.close()


def get_from_star():
	column_defocusu=None
	column_defocusv=None

	for line in lines:
		words=line.split()
		if len(words)==2:
			if words[0]=='_rlnDefocusU':
				column_defocusu=int(words[1][1:])-1
			if words[0]=='_rlnDefocusV':
				column_defocusv=int(words[1][1:])-1
			if words[0]=='_rlnDefocusAngle':
				column_defocusangle=int(words[1][1:])-1

		if len(words)>2 and column_defocusu !=None and column_defocusv!=None:
			newlines.append(words[column_defocusu]+' '+words[column_defocusv]+' '+words[column_defocusangle]+'\n')
			

	fileobj=open(outputfile,'w')
	fileobj.writelines(newlines)
	fileobj.close()


input_suffix=inputfile.split('.')[-1]
if input_suffix=='star':
	get_from_star()
elif input_suffix=='par':
	get_from_par()
