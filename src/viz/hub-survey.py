# generate histograms of connectivity survey.
import sys
import os
import matplotlib.pyplot as plt

# from https://stackoverflow.com/questions/44731152/matplotlib-create-broken-axis-in-subplot
from mpl_toolkits.axes_grid1 import make_axes_locatable

types = ['Protein','SmallMolecule','Complex']
plural = {'Protein':'Proteins',
	'SmallMolecule':'Small Molecules',
	'Complex':'Complexes',
	'Other':'Other Types'}
thres = 150

def plot_hist(file1,outprefix):
	entries = read_file(file1)
	x = []
	labels = []
	for t in types + ['Other']:
		x.append([e for e in entries[t] if e <= thres])
		num_skipped = len(entries[t])-len(x[-1])
		print('TYPE %s: %d have a threshold > %d' %(t,num_skipped,thres))
		if num_skipped > 0:
			labels.append(plural[t]+' (%d nodes > %d)' % (num_skipped,thres))
		else:
			labels.append(plural[t])

	fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6,4))
	nbins = 25
	#add_ax_hist(fig,ax,x,'Hub Survey',nbins=100,broken=True,log=True)
	ax.hist(x,bins=nbins,histtype='bar', stacked=True,edgecolor='k',log=True)
	ax.set_ylabel('# of Source Nodes')
	ax.set_xlabel('# of Reachable Nodes from the Source\'s Outgoing Hyperedges')
	ax.set_title('Hypergraph Hubs')
	plt.legend(labels)
	plt.tight_layout()
	plt.savefig(outprefix+'.png')
	print('saved to '+outprefix+'.png')
	plt.savefig(outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
	print('saved to '+outprefix+'.pdf')

	return

def add_ax_hist(fig,ax,vals,title,nbins=10,broken=True,ylabel='# of Source Nodes',log=False,shift=0.05, d = 0.01):
	vals = [sorted(v) for v in vals]
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
		flat_vals = sum(vals,[])
		i = 0
		for j in range(len(flat_vals)-1):
			if flat_vals[j+1]-flat_vals[j] >= flat_vals[i+1]-flat_vals[i]:
				i = j

		print('biggest jump:',i,flat_vals[i],flat_vals[i+1])
		jump_start = flat_vals[i]
		jump_end = flat_vals[i+1]
		# make them equal distances
		range_to_keep = max(flat_vals[i],max(flat_vals)-flat_vals[i+1])+shift

		divider = make_axes_locatable(ax)
		ax2 = divider.new_horizontal(size="100%", pad=0.2)
		fig.add_axes(ax2)

		#n,bins,patches = ax.hist(vals,bins=nbins,color='#4C8DD6',edgecolor='k',log=log)
		n,bins,patches = ax.hist(vals,bins=nbins,histtype='bar',edgecolor='k', stacked=True,log=log)
		ax.set_xlim(-range_to_keep/10,range_to_keep)
		ax.spines['right'].set_visible(False)

		#ax2.hist(vals,bins=nbins,color='#4C8DD6',edgecolor='k',log=log)
		ax2.hist(vals,bins=nbins,histtype='bar', stacked=True,edgecolor='k',log=log)
		ax2.set_xlim(max(flat_vals)-range_to_keep,max(flat_vals)+range_to_keep/10)
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
		#n,bins,patches = ax.hist(vals,bins=nbins,color='#4C8DD6',edgecolor='k',log=log)
		n,bins,patches = ax.hist(vals,bins=nbins,histtype='bar', stacked=True,edgecolor='k',log=log)

	return n,bins,patches, jump_start, jump_end


def read_file(infile):
	if not os.path.isfile(infile):
		return []
	entries = {t:[] for t in types}
	entries['Other'] = []
	with open(infile) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			else:
				row = line.strip().split('\t')
				if len(row) >= 2:
					found = False
					for t in types:
						if t in row[0]:
							found = True
							entries[t].append(int(row[2])) #  column 3 is # in Bconn set
					if not found:
						entries['Other'].append(int(row[2])) #  column 3 is # in Bconn set
			
	return entries

if __name__ == '__main__':
	if len(sys.argv) != 6:
		print('USAGE: python3 hub-survey.py <HUB_FILE> <OUTPREFIX>')

	plot_hist(sys.argv[1],sys.argv[2])