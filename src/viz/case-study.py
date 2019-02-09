## Copied anna_overlap_analysis_heatmaps.py and modified from there.
## Anna Ritz April 2018

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import sys
import scipy.cluster.hierarchy as sch
import glob

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
sorted_pathways = ['Signaling-by-MET', 'Signaling-by-Type-1-Insulin-like-Growth-Factor-1-Receptor-(IGF1R)', 
'Signaling-by-Insulin-receptor', 
'Signaling-by-EGFR', 'Signaling-by-ERBB2', 'Signaling-by-ERBB4', 'Signaling-by-SCF-KIT', 
'Signaling-by-FGFR', 'ERK1-ERK2-pathway','Signaling-by-GPCR','Signaling-by-PDGF', 
'Signaling-by-VEGF',  'DAG-and-IP3-signaling', 
'Signaling-by-NTRKs',  'PI3K-AKT-Signaling', 'Signaling-by-WNT', 'Integrin-signaling',
'TNF-signaling', 'TRAIL-signaling',  'FasL--CD95L-signaling', 
'Signaling-by-Activin', 'Signaling-by-TGF-beta-Receptor-Complex', 'Signaling-by-NOTCH', 'Signaling-by-PTK6', 
'Signaling-by-Rho-GTPases',  'MAPK6-MAPK4-signaling', 
 'p75-NTR-receptor-mediated-signalling', 'Signaling-by-Hippo',
'mTOR-signalling', 'Signaling-by-Hedgehog',
'Signaling-by-Nuclear-Receptors', 'Signaling-by-Leptin', 'Signaling-by-BMP', 'Signaling-by-MST1']

# from https://htmlcolorcodes.com/color-picker/
COLORS = ['#71FAB5','#FA7171','#95A5D5','#A8381A','#A87F1A','#0008AF','#AF0060','#47AF00','#FAB571','#A27A9B', \
'#71FAB5','#FA7171','#95A5D5','#A8381A','#A87F1A','#0008AF','#AF0060','#47AF00','#FAB571','#A27A9B']
TEXT_COLORS = ['k','k','k','w','w','w','w','k','k','k','k','k','k','w','w','w','w','k','k','k']

def main(inprefix,outprefix):  
	pathways = read_files(inprefix)
	pathway_inits = {p:pathways[p][-1] for p in pathways}
	for p in pathways:
		plist = [pathways[p][k] for k in  range(100) if k in pathways[p]]
		pathways[p] = plist

	max_k = max([len(pathways[p]) for p in sorted_pathways])
	#max_k = 20
	print('max k is ',max_k)
	pathways_to_run = ['Signaling-by-MET','Signaling-by-MST1'] # <-- list of pathways to generate these figs for
	pathway_thres = {'Signaling-by-MET':600,'Signaling-by-MST1':200}

	print('MST1:',len(pathway_inits['Signaling-by-MST1']))
	print('MET:',len(pathway_inits['Signaling-by-MET']))
	#for i in range(0,max_k):
#		influence_score_print(pathway_inits['Signaling-by-MST1'],pathway_inits['Signaling-by-MET'],pathways['Signaling-by-MST1'],pathways['Signaling-by-MET'],i)

	for p in pathways_to_run:
		print('Running ',p)
		overlap = {}
		for name in sorted_pathways:
			if name == p:
				continue
			overlap[name] = []
			cumu_p = set()
			for k in range(max_k+1):
				#print(k,len(pathways[p]),len(pathways[name]))
				if k>=len(pathways[p]) or k>=len(pathways[name]):
					overlap[name].append(overlap[name][k-1])
				else:	
					cumu_p.update(pathways[p][k])
					overlap[name].append(len(cumu_p.intersection(pathway_inits[name])))
				
			if name in pathways_to_run:
				print(name,overlap[name])
		thres = pathway_thres[p]
		pathways_to_highlight = [pa for pa in sorted_pathways if pa in overlap and overlap[pa][-1] >= thres]
		if p == 'Signaling-by-MET' and 'Signaling-by-MST1' not in pathways_to_highlight:
			pathways_to_highlight.append('Signaling-by-MST1')
		if p == 'Signaling-by-MST1' and 'Signaling-by-MET' not in pathways_to_highlight:
			pathways_to_highlight.append('Signaling-by-MET')
		print('highlighting:')
		for pa in pathways_to_highlight:
			print(pa,overlap[pa][-1])
		this_outfile_prefix = outprefix+p+'.png'
		make_figure(pathways,pathway_inits,p,this_outfile_prefix,overlap,pathways_to_highlight,thres)
		
	
	return


def influence_score_print(initp1,initp2,p1,p2,k):
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
	print('k=%d: init_intersect=%d, cumu_intersect=%d, numerator=%f, denom=%f, SCORE=%f' % (k,len(init_intersection),len(cumulative1.intersection(cumulative2)),numerator,denominator,score))
	return score

def make_figure(pathways,pathway_inits,pathway,filename,overlap,pathways_to_highlight,thres):
	if pathway == 'Signaling-by-MST1':
		plt.figure(figsize=(8,6))
		buff = 50
		scale = 0.7
	elif pathway == 'Signaling-by-MET':
		plt.figure(figsize=(8,6))
		buff = 50
		scale=0.7
	else:
		plt.figure(figsize=(7,5))
		buff = 7
		scale=0.8
	ax = plt.subplot(1,1,1)
	# from https://stackoverflow.com/questions/8971834/matplotlib-savefig-with-a-legend-outside-the-plot?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
	#box = ax.get_position()
	#print(box.width,box.height)
	#ax.set_position([.15,1-scale+0.1, box.width, box.height*scale])
	
	max_k = max([len(pathways[p]) for p in sorted_pathways])
	#max_k = 20
	x = range(max_k+1) ## max_k+1 for max_k
	i=0
	text_list = []
	for n in sorted_pathways:
		if n == pathway:
			continue
		num = len(pathway_inits[n])
		perc = int(overlap[n][-1]/num*10000)/100.0
		if n in pathways_to_highlight:
			#print(x)
			#print(overlap[n])
			ax.plot(x,overlap[n],color=COLORS[i],lw=3,label='Target: %s ($n=%d$)' % (NAMES[n],num))
			text_list.append([n,i,perc,overlap[n][-1]])
			i+=1
		else:
			ax.plot(x,overlap[n],color=[0.8,0.8,0.8],lw=1,label='_nolegend_')
	# adjust text_list to be at least 3 spaces apart
	text_list.sort(key=lambda x:x[3])
	median = int(len(text_list)/2)
	for i in range(len(text_list)):
		if i < median:
			while text_list[i][3] > text_list[i+1][3]-buff:
				for j in range(i+1):
					if text_list[j][3] > text_list[i+1][3]-(i-j+1)*buff:
						text_list[j][3] = text_list[j][3] - 1
		if i > median:
			while text_list[i][3] < text_list[i-1][3]+buff:
				for j in range(i,len(text_list)):
					if text_list[j][3] < text_list[i-1][3]+(i-j+1)*buff:
						text_list[j][3] = text_list[j][3] + 1
	for i in range(len(text_list)):
		ax.text(max_k+1,text_list[i][3],'%s (%.1f)' % (NAMES[text_list[i][0]],text_list[i][2])+'%',backgroundcolor=COLORS[text_list[i][1]],color=TEXT_COLORS[text_list[i][1]],fontsize=8)
		ax.plot([max_k-1,max_k+1],[overlap[text_list[i][0]][-1],text_list[i][3]],color=COLORS[text_list[i][1]],label='_nolegend_')

	ax.plot(range(max_k+15),[thres]*len(range(max_k+15)),'--k',label='_nolegend_')
	#ax.legend(ncol=2,bbox_to_anchor=(.9, -.15),fontsize=10)
	ax.set_title('Source Pathway %s' % (NAMES[pathway]),fontsize=14)
	ax.set_xlabel('$k$')
	ax.set_ylabel('# of Target Nodes in $kB$-Connected Set')
	ax.set_xlim(0,max_k+15)
	#plt.tight_layout()
	plt.savefig(filename)
	print('wrote to '+filename)
	filename = filename.replace('png','pdf')
	plt.savefig(filename)
	print('wrote to '+filename)
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
  

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print('USAGE: python3 case-study.py <PATHWAY_OUTPREFIX> <OUTPREFIX')
	main(sys.argv[1],sys.argv[2])