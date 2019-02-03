# generate histograms of connectivity survey.
import sys
import os
import matplotlib.pyplot as plt

# from https://stackoverflow.com/questions/44731152/matplotlib-create-broken-axis-in-subplot
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_hist(graph,compound_graph,hypergraph_nodes,hypergraph_edges,outprefix):
	graph_vals = read_file(graph)
	cgraph_vals = read_file(compound_graph)
	hgraph_vals = read_file(hypergraph_nodes)
	hedge_vals = read_file(hypergraph_edges)

	fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(10,8))

	add_ax_hist(fig,ax1,graph_vals,'Graph Connectivity',nbins=1000,broken=True)
	add_ax_hist(fig,ax2,cgraph_vals,'Compound Graph Connectivity',nbins=30,broken=False)
	add_ax_hist(fig,ax3,hgraph_vals,'Hypergraph Connectivity',nbins=30,broken=False)
	add_ax_hist(fig,ax4,hedge_vals,'Hyperedge Connectivity',nbins=30,broken=False,ylabel='# of Hyperedge Tails')

	plt.tight_layout()
	plt.savefig(outprefix+'.png')
	print('saved to '+outprefix+'.png')
	plt.savefig(outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
	print('saved to '+outprefix+'.pdf')
	return

def add_ax_hist(fig,ax,vals,title,nbins=10,broken=True,ylabel='# of Source Nodes'):
	vals = sorted(vals)
	#nbins = 1000 # number of bins
	shift = 100 # how far in to make the broken axis from last coodinates
	d = .01  # how big to make the diagonal lines in axes coordinates

	ax.set_ylabel(ylabel)
	#ax.set_xlabel('# of Reachable Nodes (%d total)' % (len(vals)))
	ax.set_xlabel('# of Reachable Nodes')
	ax.set_title(title)
	if broken:
		# get break limits
		i = 0
		for j in range(len(vals)-1):
			if vals[j+1]-vals[j] >= vals[i+1]-vals[i]:
				i = j

		print('biggest jump:',i,vals[i],vals[i+1])
		# make them equal distances
		range_to_keep = max(vals[i],max(vals)-vals[i+1])+shift

		divider = make_axes_locatable(ax)
		ax2 = divider.new_horizontal(size="100%", pad=0.2)
		fig.add_axes(ax2)

		ax.hist(vals,bins=nbins,color='#4C8DD6',edgecolor='k')
		ax.set_xlim(-10,range_to_keep)
		ax.spines['right'].set_visible(False)

		ax2.hist(vals,bins=nbins,color='#4C8DD6',edgecolor='k')
		ax2.set_xlim(max(vals)-range_to_keep,max(vals))
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
		ax.hist(vals,bins=nbins,color='#4C8DD6',edgecolor='k')

	return 


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
					entries[row[0]]= int(row[1])
	return list(entries.values())

if __name__ == '__main__':
	if len(sys.argv) != 6:
		print('USAGE: python3 connectivity-survey.py <GRAPH> <COMPOUND_GRAPH> <HYPERGRAPH_NODES> <HYPERGRAPH_EDGES> <OUTPREFIX>')

	plot_hist(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])