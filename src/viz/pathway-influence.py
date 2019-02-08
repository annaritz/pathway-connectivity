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

NAMES = {
'DAG-and-IP3-signaling':'DAG/IP3',
'ERK1-ERK2-pathway':'ERK1/ERK2',
'FasL--CD95L-signaling':'FasL/CD95L',
'Integrin-signaling':'Integrin',
'MAPK6-MAPK4-signaling':'MAPK4/MAPK6',
'mTOR-signalling':'mTOR',
'p75-NTR-receptor-mediated-signalling':'NTR',
'PI3K-AKT-Signaling':'PI3K/AKT',
'Signaling-by-Activin':'Activin',
'Signaling-by-BMP':'BMP',
'Signaling-by-EGFR':'EGFR',
'Signaling-by-ERBB2':'ERBB2',
'Signaling-by-ERBB4':'ERBB4',
'Signaling-by-FGFR':'FGFR',
'Signaling-by-GPCR':'GPCR',
'Signaling-by-Hedgehog':'Hedgehog',
'Signaling-by-Hippo':'Hippo',
'Signaling-by-Insulin-receptor':'Insulin',
'Signaling-by-Leptin':'Leptin',
'Signaling-by-MET':'MET',
'Signaling-by-MST1':'MST1',
'Signaling-by-NOTCH':'Notch',
'Signaling-by-NTRKs':'NTRKs',
'Signaling-by-Nuclear-Receptors':'Nuclear',
'Signaling-by-PDGF':'PDGF',
'Signaling-by-PTK6':'PTK6',
'Signaling-by-Rho-GTPases':'Rho GTPases',
'Signaling-by-SCF-KIT':'SCF-KIT',
'Signaling-by-TGF-beta-Receptor-Complex':'TGFB',
'Signaling-by-Type-1-Insulin-like-Growth-Factor-1-Receptor-(IGF1R)':'IGF1R',
'Signaling-by-VEGF':'VEGF',
'Signaling-by-WNT':'Wnt',
'TNF-signaling':'TNF',
'TRAIL-signaling':'TRAIL'
}

def main(inprefix,outprefix):
	pathways = read_files(inprefix)
	num = len(pathways)
	sorted_pathways = sorted(pathways.keys())
	for k in range(-1,10):
		M = []
		for i in range(num):
			M.append([0]*num)
			p1 = sorted_pathways[i]
			for j in range(num):
				p2 = sorted_pathways[j]
				#print(k,p1,p2)
				if k==-1:
					M[i][j] = asymmetric_jaccard(pathways[p1],pathways[p2])
				else:
					if i==j:
						M[i][j] = 0.0
					else:
						M[i][j] = influence_score(pathways[p1],pathways[p2],k)

		fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(4,4))
		ca = ax.matshow(M, vmin=0,vmax=1, aspect='auto') 
		fig.colorbar(ca,ax=ax)
		if k == -1:
			ax.set_title('Initial Node Overlap')
			pathway_outprefix = outprefix+'_init'
		else:
			ax.set_title('Influence Score $S_%d$' % (k))
			pathway_outprefix = outprefix+'_k_%d' % (k)
		ax.xaxis.set_ticks_position('bottom')

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
		cumulative1.update(p1[dist])
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
		print('USAGE: python3 brelax-survey.py <PATHWAY_OUTPREFIX> <OUTPREFIX')
	main(sys.argv[1],sys.argv[2])
