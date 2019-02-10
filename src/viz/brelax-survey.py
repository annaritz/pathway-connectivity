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

def main(infile,infile_filtered,outprefix):

	data1 = sort_by_col(infile)
	data2 = sort_by_col(infile_filtered)

	fig, (ax1,ax2) = plt.subplots(ncols=2, nrows=1, figsize=(8,4))
	ca1 = ax1.matshow(data1, aspect='auto') 
	ax1.set_title('Nodes with B-relaxation distance $\leq k$\n(Full Hypergraph)')
	ax1.set_xlabel('k')
	ax1.set_ylabel('Source Nodes (%d)' % len(data1))
	fig.colorbar(ca1,ax=ax1)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.set_yticklabels([])
	#fig.colorbar(ax1)
	#cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
	#cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically oriented colorbar

	#ax2.matshow(numpy.log2(data), aspect='auto') 
	ca2 = ax2.matshow(data2, aspect='auto') 
	ax2.set_title('Nodes with B-relaxation distance $\leq k$\n(Blacklisted Nodes Removed)')
	ax2.set_xlabel('k')
	ax2.set_ylabel('Source Nodes (%d)' % len(data2))
	fig.colorbar(ca2,ax=ax2)
	ax2.xaxis.set_ticks_position('bottom')
	ax2.set_yticklabels([])

	#plt.colorbar()
	plt.tight_layout()
	plt.savefig(outprefix+'.png')
	print('saved to '+outprefix+'.png')
	plt.savefig(outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
	print('saved to '+outprefix+'.pdf')
	return

def three_panel(infile,infile_filtered,infile_removed,outprefix):

	data1 = sort_by_col(infile)
	data2 = sort_by_col(infile_filtered)
	data3 = sort_by_col(infile_removed)

	fig, (ax1,ax2,ax3) = plt.subplots(ncols=3, nrows=1, figsize=(12,4))

	ca1 = ax1.matshow(numpy.ma.log10(data1).filled(0), aspect='auto') 
	ax1.set_title('Hypergraph\n$B$-relaxation Distance')
	ax1.set_xlabel('k')
	ax1.set_ylabel('Source Nodes (%d)' % len(data1))
	cbar = fig.colorbar(ca1,ax=ax1)
	cbar.ax.set_title('$\log|B_{\leq k}|$',fontsize=10)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.set_yticklabels([])
	#fig.colorbar(ax1)
	#cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
	#cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically oriented colorbar

	#ax2.matshow(numpy.log2(data), aspect='auto') 
	ca2 = ax2.matshow(numpy.ma.log10(data2).filled(0), aspect='auto') 
	ax2.set_title('Hypergraph $B$-relaxation Distance\n(Blacklisted Nodes Removed)')
	ax2.set_xlabel('k')
	ax2.set_ylabel('Source Nodes (%d)' % len(data2))
	cbar = fig.colorbar(ca2,ax=ax2)
	cbar.ax.set_title('$\log|B_{\leq k}|$',fontsize=10)
	ax2.xaxis.set_ticks_position('bottom')
	ax2.set_yticklabels([])

	ca3 = ax3.matshow(numpy.ma.log10(data3).filled(0), aspect='auto') 
	ax3.set_title('Hypergraph $B$-relaxation Distance\n(Small Molecules Removed)')
	ax3.set_xlabel('k')
	ax3.set_ylabel('Source Nodes (%d)' % len(data3))
	cbar = fig.colorbar(ca3,ax=ax3)
	cbar.ax.set_title('$\log|B_{\leq k}|$',fontsize=10)
	ax3.xaxis.set_ticks_position('bottom')
	ax3.set_yticklabels([])

	#plt.colorbar()
	plt.tight_layout()
	plt.savefig(outprefix+'.png')
	print('saved to '+outprefix+'.png')
	plt.savefig(outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
	print('saved to '+outprefix+'.pdf')
	return

def sort_by_col(infile):
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
	inds = numpy.lexsort((data[:,0],data[:,5],data[:,10],data[:,15],data[:,20],data[:,30],data[:,40]))
	data = data[inds]

	#numpy.seterr(divide='ignore')
	#data = numpy.log2(data)
	#numpy.seterr(divide='warn')# https://stackoverflow.com/questions/21752989/numpy-efficiently-avoid-0s-when-taking-logmatrix

	return data


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
	
	if len(sys.argv) != 4 and len(sys.argv) != 5:
		print('USAGE: python3 brelax-survey.py <BRELAX_FILE> <BRELAX_FILE_FILTERED> <OPTIONAL THIRD FILE> <OUTPREFIX>')

	if len(sys.argv) == 4:
		main(sys.argv[1],sys.argv[2],sys.argv[3])
	else:
		three_panel(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
