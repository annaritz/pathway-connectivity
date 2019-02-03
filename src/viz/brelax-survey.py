# generate heatmap of b_relaxation survey.
import sys
import os
import matplotlib.pyplot as plt
from math import log

# tutorial from http://blog.nextgenetics.net/?e=43
# and http://code.activestate.com/recipes/578175-hierarchical-clustering-heatmap-python/
import numpy, scipy
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist

## get data
def main(infile,outprefix):
	print('making matrix...')
	dataMatrix =[]
	with open(infile) as fin:
		for line in fin:
			if line[0] == '#': 
				continue
			row =line.strip().split('\t')
			row = [int(i) for i in row[1:]] # skip first element (hid name)
			dataMatrix.append(row)
	print('data is %d by %d' % (len(dataMatrix),len(dataMatrix[0])))

	#get your data into a 2d array where rows are genes, and columns 
	#are conditions
	data = numpy.array(dataMatrix)
	
	coarse_data = dataMatrix
	for i in range(len(dataMatrix)):
		for j in range(len(dataMatrix[i])):
			coarse_data[i][j] = int(coarse_data[i][j]*1000)/1000.0
	cdata = numpy.array(coarse_data)

	ind =numpy.lexsort((cdata[:,0],cdata[:,1]))

	data = data[ind]

	## log transform
	data = numpy.log10(data)

	fig, (ax1,ax2) = plt.subplots(ncols=2, nrows=1, figsize=(8,5))
	cax = ax1.matshow(data, aspect='auto') 
	fig.colorbar(cax)

	plt.tight_layout()
	plt.savefig(outprefix+'.png')
	print('saved to '+outprefix+'.png')
	plt.savefig(outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
	print('saved to '+outprefix+'.pdf')

if __name__ == '__main__':
	
	if len(sys.argv) != 3:
		print('USAGE: python3 brelax-survey.py <BRELAX_FILE> <OUTPREFIX>')
	main(sys.argv[1],sys.argv[2])
