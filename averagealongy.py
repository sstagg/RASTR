#! /usr/bin/env python

import mrcfile
import numpy as np
from matplotlib import pyplot
import sys

def averagey(image):
	average = np.sum(image, axis=1)
	pyplot.plot(average)
	pyplot.show()


def main():
	imagename = sys.argv[1]
	image = mrcfile.read(imagename)
	averagey(image)


if __name__ == '__main__':
	main()
