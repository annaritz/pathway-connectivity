import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os
#from viz_utils import *

def main_single(scores_file,sigs_file,outprefix):
	print('SCORES FILE',scores_file)
	print('SIGS FILE',sigs_file)

	scores,pathways = read_table(scores_file)
	sigs,ignore = read_table(sigs_file)
	print('SCORES: %d by %d' % (len(scores),len(scores[0])))
	print('SIGS: %d by %d' % (len(sigs),len(sigs[0])))

	## make grid of scatter plots.
	x = []
	y = []
	areas = []
	colors = []
	for i in range(len(pathways)):
		for j in range(len(pathways)):
			x.append(i)
			y.append(j)
			areas.append(get_area(sigs[i][j]))
			colors.append(scores[i][j])


	fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6.1,6))
	plt.gca().invert_yaxis()

	ca = ax.scatter(x,y,s=areas,c=colors,cmap=cm.get_cmap('Blues'),vmin=0,vmax=1.0, edgecolors='k',linewidths=0.1)
	ax.set_xlim(-0.5,len(pathways)-0.5)
	ax.set_ylim(len(pathways)-0.5,-0.5)
	ax.set_xticks(range(len(pathways)))
	ax.set_yticks(range(len(pathways)))
	ax.set_xticklabels([NAMES[p] for p in pathways], rotation=270, fontsize=9)
	ax.set_yticklabels([NAMES[p] for p in pathways], fontsize=9)

	fig.colorbar(ca, ax=ax, fraction=0.1, aspect=30)

	# get title
	k = int(scores_file.split('_')[-1].split('.')[0])
	ax.set_title('Influence Score $s_{%d}$' % (k))

	plt.tight_layout()
	plt.savefig(outprefix+'.png')
	print('saved to '+outprefix+'.png')
	plt.savefig(outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
	print('saved to '+outprefix+'.pdf')

	return

def main_single(scores_file,sigs_file,outprefix):
	print('SCORES FILE',scores_file)
	print('SIGS FILE',sigs_file)

	scores,pathways = read_table(scores_file)
	sigs,ignore = read_table(sigs_file)
	print('SCORES: %d by %d' % (len(scores),len(scores[0])))
	print('SIGS: %d by %d' % (len(sigs),len(sigs[0])))

	## make grid of scatter plots.
	x = []
	y = []
	areas = []
	colors = []
	for i in range(len(pathways)):
		for j in range(len(pathways)):
			x.append(i)
			y.append(j)
			areas.append(get_area(sigs[i][j]))
			colors.append(scores[i][j])


	fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6.1,6))
	plt.gca().invert_yaxis()

	ca = ax.scatter(x,y,s=areas,c=colors,cmap=cm.get_cmap('Blues'),vmin=0,vmax=1.0, edgecolors='k',linewidths=0.1)
	ax.set_xlim(-0.5,len(pathways)-0.5)
	ax.set_ylim(len(pathways)-0.5,-0.5)
	ax.set_xticks(range(len(pathways)))
	ax.set_yticks(range(len(pathways)))
	ax.set_xticklabels([NAMES[p] for p in pathways], rotation=270, fontsize=9)
	ax.set_yticklabels([NAMES[p] for p in pathways], fontsize=9)

	fig.colorbar(ca, ax=ax, fraction=0.1, aspect=30)

	# get title
	k = int(scores_file.split('_')[-1].split('.')[0])
	ax.set_title('Influence Score $s_{%d}$' % (k))

	plt.tight_layout()
	plt.savefig(outprefix+'.png')
	print('saved to '+outprefix+'.png')
	plt.savefig(outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
	print('saved to '+outprefix+'.pdf')

	return

def main_summary(scores_prefix,sigs_prefix,outprefix):
	k_range = [0,1,2,3,4,5,6,7,8,9,10,15,20,30,40]
	nrows = 5
	fig, axes_box = plt.subplots(ncols=int((len(k_range)+1)/nrows), nrows=nrows, figsize=(4,8))
	plt.gca().invert_yaxis()

	axes = []
	for a1 in axes_box:
		for a2 in a1:
			axes.append(a2)
	for i in range(len(k_range)):
		k = k_range[i]
		ax = axes[i]

		scores_file = '%s%02d.txt' % (scores_prefix,k)
		sigs_file = '%s%02d.txt' % (sigs_prefix,k)
		print('SCORES FILE',scores_file)
		print('SIGS FILE',sigs_file)

		scores,pathways = read_table(scores_file)
		sigs,ignore = read_table(sigs_file)

		## make grid of scatter plots.
		x = []
		y = []
		areas = []
		colors = []
		for i in range(len(pathways)):
			for j in range(len(pathways)):
				x.append(i)
				y.append(j)
				areas.append(get_area(sigs[i][j],5))
				colors.append(scores[i][j])
		ca = ax.scatter(x,y,s=areas,c=colors,cmap=cm.get_cmap('Blues'),vmin=0,vmax=1.0, edgecolors='k',linewidths=0.02)
		ax.set_xlim(-0.5,len(pathways)-0.5)
		ax.set_ylim(len(pathways)-0.5,-0.5)
		ax.set_xticks([])
		ax.set_yticks([])
		ax.set_title('$s_{%d}$' % (k),fontsize=14)

	plt.tight_layout()
	plt.savefig(outprefix+'.png')
	print('saved to '+outprefix+'.png')
	plt.savefig(outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
	print('saved to '+outprefix+'.pdf')
	return

def get_area(sig,factor=60):
	return (1-sig)*factor


def read_table(infile, transpose=True):
	pathways = []
	M = []
	with open(infile) as fin:
		for line in fin:
			if len(pathways) == 0:
				pathways = line.strip().split()
				continue
			M.append([float(r) for r in line.strip().split()[1:]])

	# need to transpose M.
	if transpose:
		M2 = []
		for i in range(len(M[0])):
			M2.append([0]*len(M))
			for j in range(len(M)):
				M2[i][j] = M[j][i]
		M = M2
	return M,pathways

if __name__ == '__main__':
	if len(sys.argv) != 5 or (sys.argv[1] != 'single' and sys.argv[1] != 'summary'):
		print('USAGE: python3 significant-pathway-scores.py single <SCORES_FILE> <SIGS_FILE> <OUTPREFIX>\n       python3 significant-pathway-scores.py summary <SCORES_PREFIX> <SIGS_PREFIX> <OUTPREFIX>')
		sys.exit()

	if sys.argv[1] == 'single':
		main_single(sys.argv[2],sys.argv[3],sys.argv[4])
	else:
		main_summary(sys.argv[2],sys.argv[3],sys.argv[4])