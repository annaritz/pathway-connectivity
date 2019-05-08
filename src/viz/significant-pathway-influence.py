# generate heatmap of b_relaxation survey.
import sys
import os
import matplotlib.pyplot as plt
from math import log
import glob
from matplotlib.colors import LogNorm


# tutorial from http://blog.nextgenetics.net/?e=43
# and http://code.activestate.com/recipes/578175-hierarchical-clustering-heatmap-python/
import numpy, scipy
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist

import copy # for cmap bad_valu overwrite
from matplotlib import cm

from viz_utils import *

TEST = False

def main(inprefix,outprefix,permprefix,perm,swap):
	pathways = read_files(inprefix)
	permuted_pathways = read_permuted_files(permprefix,perm,swap)
	print('%d pathways' % (len(pathways)))
	print('%d permutations' % (len(permuted_pathways)))

	num = len(sorted_pathways) # sorted_pathways is a variable from viz_utils

	if TEST:
		k_range= [3]
	else:
		k_range = [-1,0,1,2,3,4,5,6,7,8,9,10,15,20,30,40]
	all_data = []
	for k in k_range:
		print('k=%d' % (k))
		M = get_data(k,sorted_pathways,pathways,permuted_pathways,num)
		all_data.append(M)
		plot_single(M,k,'%s%d_perms_%d_swaps' % (outprefix,perm,swap),sorted_pathways,num,logvals=True,rev=True)

	#outfile = '%s%d_perms_%d_swaps_full' % (outprefix,perm,swap)
	#make_summary_plot(k_range[1:],all_data[1:],outfile,rev=True)
	
	#make_summary_plot(k_range[1:],all_data[1:],outprefix+'_full_log',logvals=True)

	return


def make_summary_plot(k_range,all_data,outprefix,logvals=False,rev=False):
	nrows = 5
	fig, axes_box = plt.subplots(ncols=int((len(k_range)+1)/nrows), nrows=nrows, figsize=(4,8))
	axes = []
	for a1 in axes_box:
		for a2 in a1:
			axes.append(a2)
	max_val = 0
	max_k = 0
	for i in range(len(all_data)):
		for j in range(len(all_data[i])):
			if max(all_data[i][j]) > max_val:
				max_val = max(all_data[i][j])
				max_k = i
			
			if logvals:
				for k in range(len(all_data[i][j])):
					all_data[i][j][k] = log(all_data[i][j][k],10)

	#print('MAX VAL IS',max_val,'AT K=',max_k)
	for i in range(len(axes)):
		if len(all_data) == i:
			axes[i].axis('off')
			continue
		ax = axes[i]
		k = k_range[i]
		M = all_data[i]
		if rev:
			if logvals:
				ca = ax.matshow(M, aspect='auto', vmin=log(1/100,10), vmax=log(1,10),cmap=plt.get_cmap('Blues_r'))
			else:
				ca = ax.matshow(M, aspect='auto', vmin=0, vmax=1,cmap=plt.get_cmap('Blues_r'))
		else:
			if logvals:
				ca = ax.matshow(M, aspect='auto', vmin=log(1/100,10), vmax=log(1,10),cmap=plt.get_cmap('Blues'))
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

def plot_single(M,k,outprefix,sorted_pathways,num,logvals=False,rev=False):
	# if logvals:
	# 	for i in range(len(M)):
	# 		for j in range(len(M[i])):
	# 			M[i][j] = log(M[i][j]+1,10)

	fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6,6))

	if rev:
		my_cmap = copy.copy(cm.get_cmap('Blues_r')) # copy the default cmap
		#https://stackoverflow.com/questions/9455044/problems-with-zeros-in-matplotlib-colors-lognorm
		my_cmap.set_bad((.5,0,.5))
	else:
		my_cmap = copy.copy(cm.get_cmap('Blues')) # copy the default cmap
	if logvals:
		norm_val = LogNorm()
		vmin = None
		vmax = None
	else:
		norm_val = None
		vmin=0
		vmax=1
			
	#ca = ax.matshow(M, aspect='auto', norm=LogNorm(), cmap=plt.get_cmap('Blues_r'))
	ca = ax.matshow(M, aspect='auto', vmin=vmin, vmax=vmax, norm=norm_val, cmap=my_cmap)
	fig.colorbar(ca,ax=ax)
	if k == -1:
		ax.set_title('Pathway Overlap Significance\n(Should be all 1)')
		pathway_outprefix = outprefix+'_init'
	else:
		ax.set_title('Influence Score $s_{%d}$ Significance' % (k))
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

def get_data(k,sorted_pathways,pathways,perm_pathways,num):
	## get cumulative values for fast influence score calculations:
	cumulative_pathways = {}
	for p in sorted_pathways:
		cumulative_pathways[p] = set()
		for dist in range(-1,k+1):
			if dist in pathways[p]:
				cumulative_pathways[p].update(pathways[p][dist])
	cumulative_perm_pathways = {}
	for perm in perm_pathways: # permutation_id
		cumulative_perm_pathways[perm] = {}
		for p in sorted_pathways:
			cumulative_perm_pathways[perm][p] = set()
			for dist in range(-1,k+1):
				if dist in perm_pathways[perm][p]:
					cumulative_perm_pathways[perm][p].update(perm_pathways[perm][p][dist])
	print('cumulatives done')

	M = []
	for i in range(num):
		M.append([0]*num)
		p1 = sorted_pathways[i]
		for j in range(num):
			p2 = sorted_pathways[j]
			
			if k==-1:
				real_val = asymmetric_jaccard(pathways[p1],pathways[p2])
				num_better = 0
				for p in perm_pathways:
					perm_val = asymmetric_jaccard(perm_pathways[p][p1],perm_pathways[p][p2])
					if perm_val >= real_val:
						num_better +=1
				M[i][j] = num_better/len(perm_pathways)
			else:
				if i==j:
					M[i][j] = 1.0 ## can't compare a pathway to itself; set to 1.
				else:
					if TEST and k==3 and p1 == 'Signaling-by-MST1' and p2=='Signaling-by-MET':
						#real_val = influence_score(pathways[p1],pathways[p2],k)
						real_val = influence_score_fast(pathways[p1][-1],pathways[p2][-1],cumulative_pathways[p1],cumulative_pathways[p2])
						print('REAL VAL:',real_val)
						num_better = 0
						for p in perm_pathways:
							#perm_val = influence_score(perm_pathways[p][p1],perm_pathways[p][p2],k)
							perm_val =  influence_score_fast(perm_pathways[p][p1][-1],perm_pathways[p][p2][-1],cumulative_perm_pathways[p][p1],cumulative_perm_pathways[p][p2])
							print('  PERM VAL:',perm_val)
							if perm_val >= real_val:
								num_better +=1
						print('NUM BETTER:',num_better)
						M[i][j] = num_better/len(perm_pathways)
						print('M[i][j]:',M[i][j])
					#real_val = influence_score(pathways[p1],pathways[p2],k)
					real_val = influence_score_fast(pathways[p1][-1],pathways[p2][-1],cumulative_pathways[p1],cumulative_pathways[p2])
					num_better = 0
					for p in perm_pathways:
						#perm_val = influence_score(perm_pathways[p][p1],perm_pathways[p][p2],k)
						perm_val =  influence_score_fast(perm_pathways[p][p1][-1],perm_pathways[p][p2][-1],cumulative_perm_pathways[p][p1],cumulative_perm_pathways[p][p2])
						if perm_val >= real_val:
							num_better +=1
					M[i][j] = num_better/len(perm_pathways)
	return M 

# def get_data(k,sorted_pathways,pathways,num):
# 	M = []
# 	for i in range(num):
# 		M.append([0]*num)
# 		p1 = sorted_pathways[i]
# 		for j in range(num):
# 			p2 = sorted_pathways[j]
			
# 			if k==-1:
# 				M[i][j] = asymmetric_jaccard(pathways[p1],pathways[p2])
# 			else:
# 				if i==j:
# 					M[i][j] = 0.0
# 				else:
# 					M[i][j] = influence_score(pathways[p1],pathways[p2],k)
# 	return M

def read_files(prefix):
	files = glob.glob(prefix+'*')
	print('%d files' % (len(files)))
	pathways = {}
	for f in files:
		pathway_name = f.replace(prefix,'').replace('_b_relax.txt','').replace('.txt','')
		print('reading %s: %s' % (pathway_name,f))
		pathways[pathway_name] = {}
		with open(f) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				row = line.strip().split('\t')
				pathways[pathway_name][int(row[0])] = set(row[2].split(';'))
	return pathways

def read_permuted_files(permprefix,perm,swap):
	perm_pathways = {}
	for p in range(perm):
		prefix = '%s%d_perms_%d_swaps_' % (permprefix,p,swap)
		#print(prefix)
		files = glob.glob(prefix+'*')
		print('%d files for permutation %d of %d' % (len(files),p,perm))
		pathways = {}
		for f in files:
			pathway_name = f.replace(prefix,'').replace('_b_relax.txt','').replace('.txt','')
			#print('reading %s' % (pathway_name))
			pathways[pathway_name] = {}
			with open(f) as fin:
				for line in fin:
					if line[0] == '#':
						continue
					row = line.strip().split('\t')
					pathways[pathway_name][int(row[0])] = set(row[2].split(';'))
		perm_pathways[p] = pathways
	return perm_pathways

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

def influence_score_fast(initp1,initp2,cumulative1,cumulative2):
	init_intersection = initp1.intersection(initp2)	
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
	
	if len(sys.argv) != 6:
		print('USAGE: python3 pathway-influence.py <PATHWAY_OUTPREFIX> <OUTPREFIX> <PATHWAY_PERM_OUTPREFIX> <PERM> <SWAP>')
	main(sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4]),int(sys.argv[5]))
