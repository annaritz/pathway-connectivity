# generate histograms of connectivity survey.
import sys
import os
import matplotlib.pyplot as plt

# from https://stackoverflow.com/questions/44731152/matplotlib-create-broken-axis-in-subplot
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_hist(file1,file2,file3,file4,outprefix):
	a = read_file(file1)
	a = [e/len(a) for e in a]
	b = read_file(file2)
	b = [e/len(b) for e in b]
	c = read_file(file3)
	c = [e/len(c) for e in c]
	d = read_file(file4)
	d = [e/len(d) for e in d]

	fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(10,8))

	n1,bins1,patches1,js1,je1 = add_ax_hist(fig,ax1,a,'Graph Connectivity\n%d nodes surveyed' % \
		(len(a)),nbins=500,broken=True, shift=0.02)
	n2,bins2,patches2,js2,je2 = add_ax_hist(fig,ax2,b,'Bipartite Graph Connectivity\n%d nodes surveyed' % \
		(len(b)),nbins=500,broken=True,shift=0.01)
	n3,bins3,patches3,js3,je3 = add_ax_hist(fig,ax3,c,'Compound Graph Connectivity\n%d nodes surveyed' % \
		(len(c)),nbins=50,broken=False)#True,shift=0.01)
	n4,bins4,patches4,js4,je4 = add_ax_hist(fig,ax4,d,'Hyperedge Connectivity\n%d nodes surveyed' % \
		(len(d)),nbins=30,broken=False)

	## overwrite jump start & end for 3 & 4
	js1,je1 = [.01,.80]
	js2,je2 = [.01,.40]
	js3,je3 = [.10,.80]
	js4,je4 = [.001,.001]

	plt.tight_layout()
	plt.savefig(outprefix+'.png')
	print('saved to '+outprefix+'.png')
	plt.savefig(outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
	print('saved to '+outprefix+'.pdf')

	print('computing statistics...')
	print_stats('graph',file1,a,js1,je1)
	print_stats('bipartite graph',file2,b,js2,je2)
	print_stats('compound graph',file3,c,js3,je3)
	print_stats('hypergraph',file4,d,js4,je4)

	return

def print_stats(name,filename,vals,js,je):
	num_vals = len(vals)
	print("FILE %s: %s.  %d total nodes surveyed." % (name,filename,num_vals))
	num_lower_tail = len([v for v in vals if v <= js])
	num_upper_tail = len([v for v in vals if v >= je])
	print(' %d source nodes (%f) <= %f reached nodes' % (num_lower_tail,num_lower_tail/num_vals,js))
	print(' %d source nodes (%f) >= %f reached nodes' % (num_upper_tail,num_upper_tail/num_vals,je))
	return

def add_ax_hist(fig,ax,vals,title,nbins=10,broken=True,ylabel='# of Source Nodes',log=False,shift=0.05, d = 0.01):
	vals = sorted(vals)
	jump_start = None
	jump_end = None
	#nbins = 1000 # number of bins
	#shift = .05 # how far in to make the broken axis from last coodinates
#	d = .01  # how big to make the diagonal lines in axes coordinates

	ax.set_ylabel(ylabel)
	#ax.set_xlabel('# of Reachable Nodes (%d total)' % (len(vals)))
	ax.set_xlabel('Percent of Reachable Nodes')
	ax.set_title(title)
	if broken:
		# get break limits
		i = 0
		for j in range(len(vals)-1):
			if vals[j+1]-vals[j] >= vals[i+1]-vals[i]:
				i = j

		print('biggest jump:',i,vals[i],vals[i+1])
		jump_start = vals[i]
		jump_end = vals[i+1]
		# make them equal distances
		range_to_keep = max(vals[i],max(vals)-vals[i+1])+shift

		divider = make_axes_locatable(ax)
		ax2 = divider.new_horizontal(size="100%", pad=0.2)
		fig.add_axes(ax2)

		n,bins,patches = ax.hist(vals,bins=nbins,color='#4C8DD6',edgecolor='k',log=log)
		ax.set_xlim(-range_to_keep/10,range_to_keep)
		ax.spines['right'].set_visible(False)

		ax2.hist(vals,bins=nbins,color='#4C8DD6',edgecolor='k',log=log)
		ax2.set_xlim(max(vals)-range_to_keep,max(vals)+range_to_keep/10)
		ax2.tick_params(left="off", labelleft='off')
		ax2.spines['left'].set_visible(False)

		# From https://matplotlib.org/examples/pylab_examples/broken_axis.html
		
		# arguments to pass to plot, just so we don't keep repeating them
		kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
		ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
		ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal

		kwargs.update(transform=ax.transAxes)  # switch to the bottom axes
		ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # bottom-left diagonal
		ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagona
	else:
		n,bins,patches = ax.hist(vals,bins=nbins,color='#4C8DD6',edgecolor='k',log=log)

	return n,bins,patches, jump_start, jump_end


def read_file(infile):
	if not os.path.isfile(infile):
		return []
	entries = {}
	with open(infile) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			else:
				row = line.strip().split('\t')
				if len(row) >= 2:
					entries[row[0]] = int(row[1])
	return list(entries.values())

if __name__ == '__main__':
	if len(sys.argv) != 6:
		print('USAGE: python3 connectivity-survey.py <GRAPH> <BIPARTITE_GRAPH> <COMPOUND_GRAPH> <HYPERGRAPH> <OUTPREFIX>')

	plot_hist(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])