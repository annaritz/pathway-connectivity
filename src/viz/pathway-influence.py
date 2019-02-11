# generate heatmap of b_relaxation survey.
import sys
import os
import matplotlib.pyplot as plt
from math import log
import glob

# tutorial from http://blog.nextgenetics.net/?e=43
# and http://code.activestate.com/recipes/578175-hierarchical-clustering-heatmap-python/
import numpy, scipy
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist

from viz_utils import *

def main(inprefix,outprefix):
	pathways = read_files(inprefix)
	num = len(sorted_pathways)
	k_range = [-1,0,1,2,3,4,5,6,7,8,9,10,15,20,30,40]
	
	#sorted_pathways = sorted(pathways.keys())
	all_data = []
	for k in k_range:
		print('k=%d' % (k))
		M = get_data(k,sorted_pathways,pathways,num)
		all_data.append(M)
		#plot_single(M,k,outprefix,sorted_pathways,num)

	make_summary_plot(k_range[1:],all_data[1:],outprefix+'_full')
	#make_summary_plot(k_range[1:],all_data[1:],outprefix+'_full_log',logvals=True)
	return


def make_summary_plot(k_range,all_data,outprefix,logvals=False):
	nrows = 5
	fig, axes_box = plt.subplots(ncols=int((len(k_range)+1)/nrows), nrows=nrows, figsize=(4,8))
	axes = []
	for a1 in axes_box:
		for a2 in a1:
			axes.append(a2)
	max_val = 0

	for i in range(len(all_data)):
		for j in range(len(all_data[i])):
			max_val = max(max_val,max(all_data[i][j]))
			if logvals:
				for k in range(len(all_data[i][j])):
					#if all_data[i][j][k] == 0:
					#	all_data[i][j][k] = 0.00000001 # small epsilon
					all_data[i][j][k] = log(all_data[i][j][k]+1,10)

	#print('MAX VAL IS',max_val)
	for i in range(len(axes)):
		if len(all_data) == i:
			axes[i].axis('off')
			continue
		ax = axes[i]
		k = k_range[i]
		M = all_data[i]
		if logvals:
			ca = ax.matshow(M, aspect='auto', vmin=log(1,10), vmax=log(2,10))
		else:
			ca = ax.matshow(M, aspect='auto', vmin=0, vmax=1,cmap=plt.get_cmap('Blues'))
		#fig.colorbar(ca,ax=ax)
		ax.set_xticks([])
		ax.set_yticks([])
		ax.set_title('$s_{%d}$' % (k),fontsize=14)

	plt.tight_layout(w_pad=0)
	plt.savefig(outprefix+'.png')
	print('saved to '+outprefix+'.png')
	plt.savefig(outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
	print('saved to '+outprefix+'.pdf')
	return

def plot_single(M,k,outprefix,sorted_pathways,num):
	fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6,6))
	ca = ax.matshow(M, aspect='auto', vmin=0.0, vmax=1,cmap=plt.get_cmap('Blues'))
	fig.colorbar(ca,ax=ax)
	if k == -1:
		ax.set_title('Pathway Overlap')
		pathway_outprefix = outprefix+'_init'
	else:
		ax.set_title('Influence Score $s_{%d}$' % (k))
		pathway_outprefix = outprefix+'_k_%s' % (str(k).zfill(2))

	ax.xaxis.set_ticks_position('bottom')
	ax.set_xticks(range(num))
	ax.set_yticks(range(num))
	ax.set_xticklabels([NAMES[p] for p in sorted_pathways], rotation=270, fontsize=9)
	ax.set_yticklabels([NAMES[p] for p in sorted_pathways], fontsize=9)

	plt.tight_layout()
	plt.savefig(pathway_outprefix+'.png')
	print('saved to '+pathway_outprefix+'.png')
	plt.savefig(pathway_outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (pathway_outprefix,pathway_outprefix))
	print('saved to '+pathway_outprefix+'.pdf')

	out = open(pathway_outprefix+'.txt','w')
	out.write('\t'.join(sorted_pathways)+'\n')
	for i in range(num):
		out.write('%s\t%s\n' % (sorted_pathways[i],'\t'.join([str(e) for e in M[i]])))
	out.close()
	print('saved output to %s'  %(pathway_outprefix+'.txt'))
	return

def get_data(k,sorted_pathways,pathways,num):
	M = []
	for i in range(num):
		M.append([0]*num)
		p1 = sorted_pathways[i]
		for j in range(num):
			p2 = sorted_pathways[j]
			
			if k==-1:
				# if 'Leptin' in p1 and 'IGF' in p2:
				# 	print(p1,len(pathways[p1][-1]))
				# 	print(p2,len(pathways[p2][-1]))
				# 	print(len(pathways[p1][-1].intersection(pathways[p2][-1])))
				# 	print(asymmetric_jaccard(pathways[p1],pathways[p2]))
				# 	sys.exit()
				M[i][j] = asymmetric_jaccard(pathways[p1],pathways[p2])
			else:
				if i==j:
					M[i][j] = 0.0
				else:
					M[i][j] = influence_score(pathways[p1],pathways[p2],k)
	return M

def read_files(prefix):
	files = glob.glob(prefix+'*')
	print('%d files' % (len(files)))
	pathways = {}
	for f in files:
		pathway_name = f.replace(prefix,'').replace('_b_relax.txt','')
		print('reading %s' % (pathway_name))
		pathways[pathway_name] = {}
		with open(f) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				row = line.strip().split('\t')
				pathways[pathway_name][int(row[0])] = set(row[2].split(';'))
	return pathways

def influence_score(p1,p2,k):
	initp1 = p1[-1]
	initp2 = p2[-1]
	init_intersection = initp1.intersection(initp2)

	cumulative1 = set()
	cumulative2 = set()
	for dist in range(-1,k+1):
		if dist in p1:
			cumulative1.update(p1[dist])
		if dist in p2:
			cumulative2.update(p2[dist])

	numerator = len(cumulative1.intersection(initp2).difference(init_intersection))
	denominator = len(cumulative1.difference(init_intersection))
	score = numerator/denominator
	if score < 0 or score > 1:
		print(score)
		sys.exit()
	return score

def asymmetric_jaccard(p1,p2):
	initp1 = p1[-1]
	initp2 = p2[-1]
	jaccard = len(initp1.intersection(initp2))/len(initp1)
	return jaccard

if __name__ == '__main__':
	
	if len(sys.argv) != 3:
		print('USAGE: python3 pathway-influence.py <PATHWAY_OUTPREFIX> <OUTPREFIX')
	main(sys.argv[1],sys.argv[2])
