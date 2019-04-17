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


def three_panel(graph_file,bipartite_file,compound_file,outprefix):

	data1 = sort_by_col(graph_file)
	data2 = sort_by_col(bipartite_file)
	data3 = sort_by_col(compound_file)

	fig, (ax1,ax2,ax3) = plt.subplots(ncols=3, nrows=1, figsize=(12,4))

	ca1 = ax1.matshow(numpy.ma.log10(data1).filled(0), aspect='auto') 
	ax1.set_title('Graph Connectivity')
	ax1.set_xlabel('Distance')
	ax1.set_ylabel('Source Nodes')
	cbar = fig.colorbar(ca1,ax=ax1)
	cbar.ax.set_title('$\log_{10}$(# of Nodes)',fontsize=10)
	ax1.xaxis.set_ticks_position('bottom')
	#ax1.set_yticklabels([])
	#fig.colorbar(ax1)
	#cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
	#cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically oriented colorbar

	#ax2.matshow(numpy.log2(data), aspect='auto') 
	ca2 = ax2.matshow(numpy.ma.log10(data2).filled(0), aspect='auto') 
	ax2.set_title('Bipartite Graph Connectivity')
	ax2.set_xlabel('Distance')
	ax2.set_ylabel('Source Nodes')
	cbar = fig.colorbar(ca2,ax=ax2)
	cbar.ax.set_title('$\log_{10}$(# of Nodes)',fontsize=10)
	ax2.xaxis.set_ticks_position('bottom')
	#ax2.set_yticklabels([])

	ca3 = ax3.matshow(numpy.ma.log10(data3).filled(0), aspect='auto') 
	ax3.set_title('Compound Graph Connectivity')
	ax3.set_xlabel('Distance')
	ax3.set_ylabel('Source Nodes')
	cbar = fig.colorbar(ca3,ax=ax3)
	cbar.ax.set_title('$\log_{10}$(# of Nodes)',fontsize=10)
	ax3.xaxis.set_ticks_position('bottom')
	#ax3.set_yticklabels([])

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
	max_len = -1
	even = True
	with open(infile) as fin:
		for line in fin:
			if line[0] == '#': 
				continue
			row =line.strip().split('\t')
			row = [int(i) for i in row[1:]] # skip first element (hid name)
			dataMatrix.append(row)
			if max_len == -1:
				max_len = len(row)
			else:
				if max_len != len(row):
					even = False
					max_len = max(max_len,len(row))
	if not even: # pad 0's
		for i in range(len(dataMatrix)):
			if len(dataMatrix[i]) != max_len:
				if len(dataMatrix[i]) == 0:
					dataMatrix[i] = [0]*max_len
				else:
					dataMatrix[i] += [dataMatrix[i][-1]]*(max_len-len(dataMatrix[i]))
				if len(dataMatrix[i]) != max_len:
					print('ERROR: lengths still not even.')

	print('data is %d by %d' % (len(dataMatrix),len(dataMatrix[0])))

	#get your data into a 2d array where rows are genes, and columns 
	#are conditions
	data = numpy.array(dataMatrix)
	jump = max_len//2//5
	print('jump is %d' % (jump))
	inds = numpy.lexsort((data[:,0],data[:,jump],data[:,2*jump],data[:,3*jump],data[:,4*jump]))
	data = data[inds]

	#numpy.seterr(divide='ignore')
	#data = numpy.log2(data)
	#numpy.seterr(divide='warn')# https://stackoverflow.com/questions/21752989/numpy-efficiently-avoid-0s-when-taking-logmatrix

	return data

if __name__ == '__main__':
	
	## hard-coded for now.
	graph_file = '../SIF/outfiles/reactome.txt.parameterized'
	bipartite_file = '../hypergraph_code/output/graph-reactome.txt.parameterized'
	compound_file = '../BioPAXSTREAM/output/reactome_parameterized_filtered.txt'
	outprefix = 'parameterized-connectivity-survey'
	three_panel(graph_file,bipartite_file,compound_file,outprefix)
