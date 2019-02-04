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

def main(infile,outprefix):
	print('making matrix...')
	dataMatrix =[]
	with open(infile) as fin:
		for line in fin:
			if line[0] == '#': 
				continue
			row =line.strip().split('\t')
			row = [int(i) for i in row[2:]] # skip first element (hid name) and second element (time)
			dataMatrix.append(row)
	print('data is %d by %d' % (len(dataMatrix),len(dataMatrix[0])))

	#get your data into a 2d array where rows are genes, and columns 
	#are conditions
	data = numpy.array(dataMatrix)

	# only get the first five columns
	inds = numpy.lexsort((data[:,0],data[:,5],data[:,10]))

	data = data[inds]

	numpy.seterr(divide='ignore')
	data = numpy.log10(data)
	numpy.seterr(divide='warn')# https://stackoverflow.com/questions/21752989/numpy-efficiently-avoid-0s-when-taking-logmatrix

	fig, (ax1,ax2) = plt.subplots(ncols=2, nrows=1, figsize=(5,8))
	ax1.matshow(data, aspect='auto', origin='lower') 
	ax2.matshow(data, aspect='auto', origin='lower') 
	#plt.colorbar()
	plt.tight_layout()
	plt.savefig(outprefix+'.png')
	print('saved to '+outprefix+'.png')
	plt.savefig(outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
	print('saved to '+outprefix+'.pdf')
	return

def clustering(infile,outprefix):
	print('making matrix...')
	dataMatrix =[]
	with open(infile) as fin:
		for line in fin:
			if line[0] == '#': 
				continue
			row =line.strip().split('\t')
			row = [int(i) for i in row[2:]] # skip first element (hid name)
			dataMatrix.append(row)
	print('data is %d by %d' % (len(dataMatrix),len(dataMatrix[0])))

	# pull off all zeros
	num_zero_rows = sum([dataMatrix[i][-1]==0 for i in range(len(dataMatrix))])
	print('there are %d zeroed-out rows. Ignore.' % (num_zero_rows))
	dataMatrix = [dataMatrix[i] for i in range(len(dataMatrix)) if dataMatrix[i][-1]!=0 ]
	print('data is now %d by %d' % (len(dataMatrix),len(dataMatrix[0])))

	#get your data into a 2d array where rows are genes, and columns 
	#are conditions
	data = numpy.array(dataMatrix)

	# only get the first five columns
	smalldata = data[:100,:10]

	print('calculating distances...')
	#calculate a distance matrix
	distMatrix = dist.pdist(smalldata)

	#convert the distance matrix to square form. The distance matrix 
	#calculated above contains no redundancies, you want a square form 
	#matrix where the data is mirrored across the diagonal.
	distSquareMatrix = dist.squareform(distMatrix)

	print('linkage matrix...')
	#calculate the linkage matrix 
	linkageMatrix = hier.linkage(distSquareMatrix)

	print('dendrogram...')
	dendro = hier.dendrogram(linkageMatrix)

	#get the order of rows according to the dendrogram 
	leaves = dendro['leaves'] 
	print(leaves)

	print('reorder...')
	#reorder the original data according to the order of the 
	#dendrogram. Note that this slice notation is numpy specific.
	#It just means for every row specified in the 'leaves' array,
	#get all the columns. So it effectively remakes the data matrix
	#using the order from the 'leaves' array.
	transformedData = data[leaves,:]

	print('log transform...')
	## log transform
	for i in range(len(dataMatrix)):
		for j in range(len(dataMatrix[i])):
			if dataMatrix[i][j] != 0:
				dataMatrix[i][j] = log(dataMatrix[i][j],10)
	numpy.seterr(divide='ignore')
	transformedData = numpy.log10(transformedData)
	numpy.seterr(divide='warn')# https://stackoverflow.com/questions/21752989/numpy-efficiently-avoid-0s-when-taking-logmatrix

	fig, (ax1,ax2) = plt.subplots(ncols=2, nrows=1, figsize=(5,8))
	ax1.matshow(dataMatrix, aspect='auto', origin='lower') 
	ax2.matshow(transformedData, aspect='auto', origin='lower') 

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
