# generate heatmap of b_relaxation survey.
import sys
import os
import matplotlib.pyplot as plt
from math import log
from scipy.misc import comb

# tutorial from http://blog.nextgenetics.net/?e=43
# and http://code.activestate.com/recipes/578175-hierarchical-clustering-heatmap-python/
import numpy, scipy
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist


def cumulative_histogram(graph_file,bipartite_file,compound_file,hgraph_file,outprefix):

	#fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(10,8))
	#plot_ax(ax1,graph_file,bipartite_file,compound_file,hgraph_file,False,False,'Connectivity')
	#plot_ax(ax2,graph_file,bipartite_file,compound_file,hgraph_file,True,False,'Connectivity')
	#plot_ax(ax3,graph_file,bipartite_file,compound_file,hgraph_file,False,True,'Cumulative Connectivity')
	#plot_ax(ax4,graph_file,bipartite_file,compound_file,hgraph_file,True,True,'Cumulative Connectivity')

	fig, (ax1,ax2,ax3) = plt.subplots(ncols=3, nrows=1, figsize=(12,4))
	plot_ax(ax1,graph_file,bipartite_file,compound_file,hgraph_file,False,False,'Unreachable s-t Paths',zero=True)
	plot_ax(ax2,graph_file,bipartite_file,compound_file,hgraph_file,False,False,'Reachable s-t Paths')
	plot_ax(ax3,graph_file,bipartite_file,compound_file,hgraph_file,True,True,'Cumulative Reachable s-t Paths')

	plt.tight_layout()
	plt.savefig(outprefix+'.png')
	print('saved to '+outprefix+'.png')
	plt.savefig(outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
	print('saved to '+outprefix+'.pdf')
	return

def plot_ax(ax,graph_file,bipartite_file,compound_file,hgraph_file,norm,cumulative,title,zero=False):

	data1,num1,dist1,zeros1 = get_hist(graph_file,norm=norm, cumulative=cumulative)
	data2,num2,dist2,zeros2 = get_hist(bipartite_file,norm=norm, cumulative=cumulative)
	data3,num3,dist3,zeros3 = get_hist(compound_file,norm=norm, cumulative=cumulative)
	data4,num4,dist4,zeros4 = get_hist(hgraph_file,norm=norm, cumulative=cumulative)

	if zero:
		#print('** ZEROS **')
		#print(zeros1,zeros2,zeros3,zeros4)
		#print(num1,num2,num3,num4)
		#print(num1**2,num2**2,num3**2,num4**2)
		#print(zeros1/num1**2,zeros2/num2**2,zeros3/num3**2,zeros4/num4**2)
		ax.bar([0],[zeros1/float(num1**2)],.6)
		ax.bar([1],[zeros3/float(num3**2)],.6)
		ax.bar([2],[zeros2/float(num2**2)],.6)
		ax.bar([3],[zeros4/float(num4**2)],.6)
		ax.set_title(title)
		ax.set_ylabel('Proportion of Unreachable s-t Paths')
		ax.set_ylim(0,1)
		ax.set_xticks([0,1,2,3])
		ax.set_xticklabels(['Graph\n(%.3f)' % (zeros1/num1**2),'Compound\nGraph\n(%.3f)' % (zeros3/num3**2),'Bipartite\nGraph\n(%.3f)' % (zeros2/num2**2),'Hypergraph\n(%.5f)' % (zeros4/num4**2)])
		return

	# ax.plot(data1,'-o',lw=1,ms=5,label='Graph\n(%d nodes)' % (num1))
	# ax.plot(data3,'-o',lw=1,ms=5,label='Compound Graph\n(%d nodes)' % (num3))
	# ax.plot(data2,'-o',lw=1,ms=5,label='Bipartite Graph\n(%d nodes)' % (num2))
	# ax.plot(data4,'-o',lw=1,ms=5,label='Hypergraph\n(%d nodes)' % (num4))

	ax.plot(data1,'-',lw=3,label='Graph\n(%d nodes)' % (num1))
	ax.plot(data3,'-',lw=3,label='Compound Graph\n(%d nodes)' % (num3))
	ax.plot(data2,'-',lw=3,label='Bipartite Graph\n(%d nodes)' % (num2))
	ax.plot(data4,'-',lw=3,label='Hypergraph\n(%d nodes)' % (num4))
	
	ax.set_title(title)

	ax.set_xlabel('Distance')
	if norm:
		ax.set_ylabel('Proportion of s-t Paths')
		ax.set_ylim(-.01,1.01)
	else:
		ax.set_ylabel('Number of s-t Paths')
	ax.legend()
	return

def get_hist(infile,norm=False,cumulative=False):
	print('making histogram from %s...' % infile)

	
	# first get max distance - this is used later.
	max_dist =0
	num_datapoints = 0
	with open(infile) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			num_datapoints+=1
			max_dist = max(max_dist,len(line.strip().split())-1)

	data = [0]*max_dist
	num_zeros = 0
	with open(infile) as fin:
		for line in fin:
			if line[0] == '#': 
				continue
			
			row =line.strip().split('\t')
			row = [int(i) for i in row[1:]] # skip first element (hid name)

			for i in range(len(row)):  
				if i==0:
					data[i] += row[i]
				else: # row is cumulative -- "undo" that for this.
					data[i]+=row[i]-row[i-1]
			# count number of zeros. Row is cumulative so we can use that.
			if len(row)==0:
				num_zeros+=num_datapoints-1 # ignore THe source node
			else:
				num_zeros += num_datapoints - row[-1]

	## append num_zeros to data
	## data = [num_zeros] + data
			
	print('ORIG:',data[:10])
	if cumulative:
		for i in range(1,len(data)):
			data[i]+=data[i-1]
	print('POST-CUMULATIVE:',data[:10])
	if norm:
		for d in range(len(data)):
			data[d] = data[d]/float(num_datapoints**2)
	print('POST-NORM',data[:10])
	#sys.exit()
	return data, num_datapoints, max_dist, num_zeros

if __name__ == '__main__':
	
	## hard-coded for now.
	graph_file = '../SIF/outfiles/reactome.txt.parameterized'
	bipartite_file = '../hypergraph_code/output/graph-reactome.txt.parameterized'
	compound_file = '../BioPAXSTREAM/output/reactome_parameterized_filtered.txt'
	hgraph_file = '../hypergraph_code/output/reactome-bconn-parameterized.txt'
	outprefix = 'cumulative-histogram'
	cumulative_histogram(graph_file,bipartite_file,compound_file,hgraph_file,outprefix)
