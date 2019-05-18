from __future__ import print_function

## packages
import argparse
import sys
import os
import networkx as nx  ## TODO verify version 
import time
import random
import glob
from multiprocessing import Pool ## https://docs.python.org/3.4/library/multiprocessing.html?highlight=process
import matplotlib.pyplot as plt
from matplotlib import cm

## custom code
import SIF.graph_utils as graph_utils
import hypergraph_code.hgraph_utils as hgraph_utils
import viz.viz_utils as viz_utils
import viz.cumulative_histogram as cumulative_histogram
import viz.connectivity_survey_parameterized as heatmap_viz
import viz.significant_pathway_scores as permutation_viz
import STRING_channels.run_channels as run_channels
import STRING_channels.viz_channels as viz_channels

# TODO cplex is broken for python 3.7 on home laptop
import hypergraph_code.ILP.shortest_hyperpath as shortest_hyperpath 
import hypergraph_code.permutation_test as permutation_test
## halp
from halp.utilities import directed_statistics as stats
from halp.algorithms import directed_paths as hpaths

## TODO: make these options or a config file
SIF_FILE = 'SIF/PathwayCommons10.reactome.hgnc.sif'
SIF_CONV_FILE = 'SIF/conversion-types.txt'
HGRAPH_PREFIX = '../hypergraph/reactome_hypergraph_full/reactome'
HGRAPH_BLACKLIST_PREFIX = '../hypergraph/reactome_hypergraph_full/blacklist_filter'
HGRAPH_SMALLMOL_PREFIX = '../hypergraph/reactome_hypergraph_full/small_molecule_filter'
PATHWAY_HGRAPH_ENTITIES = '../hypergraph/reactome_hypergraphs/'
PATHWAY_HGRAPH_PREFIX = '../hypergraph/reactome_hypergraphs_parsed/'
BLACKLIST_FILE = '../data/blacklist.txt'
PROCESSED_STRING_DIR = '../data/STRING/processed/' ## TODO incorporate code for processing string

## outdirectories
OUT_TXT_DIR = 'out_txt/'
OUT_PERM_DIR = 'out_txt/permutations/'
OUT_PATHWAY_DIR = 'out_txt/pathway_survey/'
OUT_TXT_CHANNEL_DIR = 'out_txt/STRING_channels/'
OUT_VIZ_DIR = 'out_viz/'
OUT_VIZ_CHANNEL_DIR = 'out_viz/STRING_channels/'
DIRS = [OUT_TXT_DIR,OUT_PERM_DIR,OUT_PATHWAY_DIR,OUT_VIZ_DIR,OUT_TXT_CHANNEL_DIR,OUT_VIZ_CHANNEL_DIR]
for d in DIRS:
	if not os.path.isdir(d):
		os.system('mkdir %s' % (d))

## global variables (will be set with parse_options)
FORCE = False
PRINT_ONLY = False

PATHWAYS = viz_utils.sorted_pathways 
PATHWAY_NAMES = viz_utils.NAMES

LARGE_VAL = 10000000

#############################
## MAIN FUNCTION
#############################
def main():

	opts = parse_options()
	print(opts)
	print()

	## get representations no matter what
	sif_graph,compound_graph,bipartite_graph,compact_bipartite_graph,hypergraph,identifier2id,id2identifier = \
		get_representations(opts.small_molecule_filter,opts.blacklist_filter,opts.keep_singletons)

	## TODO: add hypergraph building into this script.
	## TODO CHECK WHAT HAPPENS IF WE KEEP SINGLETON NODES -- SO many singleton nodes oof.

	######### Print statistics for table
	if opts.stats: 
		print()
		print('------ STATS -------')
		print('KEEP_SINGLETONS:',opts.keep_singletons,'SMALL_MOL FILTER:',opts.small_molecule_filter,'BLACKLIST FILTER:',opts.blacklist_filter)
		print('SIF Graph: %d nodes and %d edges' % (nx.number_of_nodes(sif_graph),nx.number_of_edges(sif_graph)))
		print('Bipartite Graph: %d nodes and %d edges' % (nx.number_of_nodes(bipartite_graph),nx.number_of_edges(bipartite_graph)))
		print('Compact Bipartite Graph: %d nodes and %d edges' % (nx.number_of_nodes(compact_bipartite_graph),nx.number_of_edges(compact_bipartite_graph)))
		print('Hypergraph: %d nodes and %d hyperedges' % (stats.number_of_nodes(hypergraph),stats.number_of_hyperedges(hypergraph)))
		print('------ END STATS -------')
		print()

	########## Make survey histograms
	if opts.histograms:
		## make parameterized files if they don't exist (or --force is specified)
		sif_graph_file = make_outfile(opts,OUT_TXT_DIR,'sif-graph-survey')
		survey_graph(sif_graph,sif_graph_file)
		bipartite_graph_file = make_outfile(opts,OUT_TXT_DIR,'bipartite-graph-survey')
		survey_graph(bipartite_graph,bipartite_graph_file)
		## TODO: this is broken. Copy files from work laptop over before continuing.
		hypergraph_file = make_outfile(opts,OUT_TXT_DIR,'hypergraph-survey')
		survey_hgraph(hypergraph,hypergraph_file)
		brelax_file = make_outfile(opts,OUT_TXT_DIR,'hypergraph-brelax-survey')
		survey_hgraph_brelax(hypergraph,brelax_file,id2identifier,identifier2id,opts)

		viz_histograms(sif_graph_file,bipartite_graph_file,bipartite_graph_file,hypergraph_file,brelax_file,opts)

	########## read pathway file from hypergraph in certain cases
	if opts.perm_test or opts.case_studies:
		## Get pathway files
		pathwayfile = make_outfile(opts,OUT_TXT_DIR,'pathways-from-hypergraph')
		if FORCE or not os.path.isfile(pathwayfile):
			make_pathways_from_hypergraph(hypergraph,pathwayfile)
		else:
			force_print_statement(pathwayfile)
		pathway_members = read_pathway_members(pathwayfile)

	########## Run permutation tests
	if opts.perm_test:

		## read SIF pathway file too.
		sif_pathwayfile = make_outfile(opts,OUT_TXT_DIR,'pathways-from-sif')
		if FORCE or not os.path.isfile(sif_pathwayfile):
			make_pathways_from_sif_graph(sif_graph,pathwayfile,sif_pathwayfile)
		else:
			force_print_statement(pathwayfile)
		sif_pathway_members = read_pathway_members(sif_pathwayfile)
		
		# calculate permutation tests with 10,000 swaps each. 
		num_swaps = 10000 
		generate_pathway_permutations(sif_pathway_members,num_swaps,opts,'sif_graph')
		generate_pathway_permutations(pathway_members,num_swaps,opts,'hypergraph')

		####### HYPERGRAPH
		# calculate original pathway survey first
		survey_hgraph_pathways(hypergraph,id2identifier,identifier2id,pathway_members,'hypergraph',opts,verbose=True)
		scores_file = make_outfile(opts,OUT_PATHWAY_DIR,'hypergraph')

		# run permutation surveys
		num_threads = 3
		args_to_run = []
		for permutation in range(opts.perm_test):
			suffix = 'hypergraph_%d_perms_%d_swaps' % (permutation,num_swaps)
			outfile = make_outfile(opts,OUT_PATHWAY_DIR,suffix)
			if FORCE or not os.path.isfile(outfile):
				permfile = make_outfile(opts,OUT_PERM_DIR,'hypergraph_%d_perms_%d_swaps' % (permutation,num_swaps))
				perm_pathway_members = read_pathway_members(permfile)
				args_to_run.append((hypergraph,id2identifier,identifier2id,perm_pathway_members,suffix,opts))
		if len(args_to_run) > 0:
			print('Running %d args with %d threads' % (len(args_to_run),num_threads))
			with Pool(num_threads) as p:
				p.map(survey_hgraph_pathways_threaded,args_to_run)
		else:
			print('not running any args. Use --force to run.')

		viz_permutations(scores_file,'hypergraph',opts.perm_test,opts,k_vals=[0,1,2,3,4,5,10,15,20,25,30])

		####### BIPARTITE GRAPH
		# calculate original bipartite pathway survey first
		survey_graph_pathways(compact_bipartite_graph,pathway_members,'bipartite',opts,verbose=True)
		scores_file = make_outfile(opts,OUT_PATHWAY_DIR,'bipartite')

		args_to_run = []
		for permutation in range(opts.perm_test):
			suffix = 'bipartite_%d_perms_%d_swaps' % (permutation,num_swaps)
			outfile = make_outfile(opts,OUT_PATHWAY_DIR,suffix)
			if FORCE or not os.path.isfile(outfile):
				permfile = make_outfile(opts,OUT_PERM_DIR,'hypergraph_%d_perms_%d_swaps' % (permutation,num_swaps))
				perm_pathway_members = read_pathway_members(permfile)
				args_to_run.append((compact_bipartite_graph,perm_pathway_members,suffix,opts,True))
		if len(args_to_run) > 0:
			print('Running %d args with %d threads' % (len(args_to_run),num_threads))
			with Pool(num_threads) as p:
				p.map(survey_graph_pathways_threaded,args_to_run)
		else:
			print('not running any args. Use --force to run.')

		viz_permutations(scores_file,'bipartite',opts.perm_test,opts,k_vals=[0,1,2,3,4,5,10,20,30,40,50,60])
		########### SIF GRAPH

		# calculate original SIF pathway survey first
		survey_graph_pathways(sif_graph,sif_pathway_members,'sif_graph',opts,verbose=True)
		scores_file = make_outfile(opts,OUT_PATHWAY_DIR,'sif_graph')

		args_to_run = []
		for permutation in range(opts.perm_test):
			suffix = 'sif_graph_%d_perms_%d_swaps' % (permutation,num_swaps)
			outfile = make_outfile(opts,OUT_PATHWAY_DIR,suffix)
			if FORCE or not os.path.isfile(outfile):
				permfile = make_outfile(opts,OUT_PERM_DIR,'sif_graph_%d_perms_%d_swaps' % (permutation,num_swaps))
				perm_pathway_members = read_pathway_members(permfile)
				args_to_run.append((sif_graph,perm_pathway_members,suffix,opts,False))
		if len(args_to_run) > 0:
			print('Running %d args with %d threads' % (len(args_to_run),num_threads))
			with Pool(num_threads) as p:
				p.map(survey_graph_pathways_threaded,args_to_run)
		else:
			print('not running any args. Use --force to run.')

		viz_permutations(scores_file,'sif_graph',opts.perm_test,opts,k_vals=[0,1,2,3,4,5,6,7,8,9,10])

	## visualize case studies for hypergraph
	if opts.case_studies:
		# run two batches
		pathways_to_highlight = ['Signaling-by-Activin', 'Signaling-by-TGF-beta-Receptor-Complex','Signaling-by-BMP']
		viz_influence_histogram(hypergraph,id2identifier,identifier2id,pathway_members,pathways_to_highlight,opts,[50,250,50])

		pathways_to_highlight = ['Signaling-by-MET','Signaling-by-MST1','ERK1-ERK2-pathway','PI3K-AKT-Signaling']
		viz_influence_histogram(hypergraph,id2identifier,identifier2id,pathway_members,pathways_to_highlight,opts,[250,50,250,300])

	## run STRING channels (all 200 signaling pathways)
	if opts.string_channels:
		pathway_nodes,all_pathway_nodes = get_pathways_for_STRING_analysis(opts)

		## add entity set info and get node memberships
		hypergraph = hgraph_utils.add_entity_set_info(hypergraph)
		nodes,node_membership = get_node_memberships(hypergraph)

		# get pathway Identifiers to Uniprot ID
		pc2uniprot,uniprot2pc = hgraph_utils.get_id_map(PATHWAY_HGRAPH_ENTITIES)

		files = glob.glob('%s*.txt' % (PROCESSED_STRING_DIR))
		processed_nodes = {}
		for f in files:
			name,interactions = get_STRING_channel_interactions(f)

			outfile = make_outfile(opts,OUT_TXT_CHANNEL_DIR,name)
			if FORCE or not os.path.isfile(outfile):
				interactions_in_reactome = get_STRING_interactions_in_reactome(interactions,uniprot2pc,name,nodes,opts)
				interactions_in_pathways,interactions_in_same_pathway = run_channels.get_pathway_interactions(interactions_in_reactome,pathway_nodes,all_pathway_nodes)
				print('  %d INTERACTIONS HAVE BOTH NODES IN THE REACTOME PATHWAYS' % (len(interactions_in_pathways)))
				print('  %d INTERACTIONS HAVE BOTH NODES IN SAME REACTOME PATHWAY' % (len(interactions_in_same_pathway)))
				sys.stdout.flush()

				b_visit_dict = get_bvisit_dict(hypergraph,id2identifier,identifier2id,opts)
				brelax_dicts,processed_nodes = run_channels.preprocess_brelax_dicts(hypergraph,interactions_in_pathways,node_membership,b_visit_dict,processed_nodes)
				interactions_brelax = run_channels.get_bconn_interactions(brelax_dicts,interactions_in_pathways,node_membership)
				interactions_bipartite = list(interactions_brelax.keys())
				interactions_bconn = [e for e in interactions_bipartite if interactions_brelax[e] == 0]

				print('  %d INTERACTIONS ARE Bipartite CONNECTED IN REACTOME' % (len(interactions_bipartite)))
				print('  %d INTERACTIONS ARE B-CONNECTED IN REACTOME' % (len(interactions_bconn)))
				sys.stdout.flush()
				
				write_channel_output(interactions_in_reactome,interactions_in_pathways,interactions_in_same_pathway,interactions_brelax,outfile)
				## ensure that variables are the right type for viz
				interactions_in_reactome = {(n1,n2):val for n1,n2,val in interactions_in_reactome} 
				interactions_in_pathways = set(interactions_in_pathways)
				interactions_in_same_pathway = set(interactions_in_same_pathway)
			else:
				force_print_statement(outfile)
				interactions_in_reactome = {}
				interactions_in_pathways = set()
				interactions_in_same_pathway = set()
				interactions_brelax = {}
				with open(outfile) as fin:
					for line in fin:
						if line[0] == '#':
							continue
						row = line.strip().split()
						key = (row[0],row[1])
						#    0      1       2       3           4             5         6
						# '#Node1\tNode2\tScore\tAnyPathway\tSamePathway\tBipartite\tBRelaxDist\n'
						interactions_in_reactome[key] = int(row[2])
						if row[3] == '1':
							interactions_in_pathways.add(key)
						if row[4] == '1':
							interactions_in_same_pathway.add(key)
						if row[5] == '1':
							interactions_brelax[key] = int(row[6])
			connected_set = set(interactions_brelax.keys())
			bconn_set = set([key for key in interactions_brelax.keys() if interactions_brelax[key]==0])
			outfile = make_outfile(opts,OUT_VIZ_CHANNEL_DIR,name,filetype='')
			viz_channels.viz(interactions_in_reactome,[interactions_in_pathways,interactions_in_same_pathway,connected_set,bconn_set],\
			['Pair in\nAny Pathway','Pair in\nSame Pathway','Pair\nConnected','Pair\nB-Connected'],\
			['Pair in Any Pathway','Pair in Same Pathway','Pair Connected','Pair B-Connected'],outfile,name,brelax=interactions_brelax)
	return


#############################
## STRING CHANNEL FUNCTIONS
#############################

def write_channel_output(interactions_in_reactome,interactions_in_pathways,interactions_in_same_pathway,interactions_brelax,outfile):
	out = open(outfile,'w')
	out.write('#Node1\tNode2\tScore\tAnyPathway\tSamePathway\tBipartite\tBRelaxDist\n')
	for n1,n2,val in interactions_in_reactome:
		a = 1 if (n1,n2) in interactions_in_pathways else 0
		b = 1 if (n1,n2) in interactions_in_same_pathway else 0
		c = 1 if (n1,n2) in interactions_brelax else 0
		d = interactions_brelax[(n1,n2)] if (n1,n2) in interactions_brelax else -1
		out.write('%s\t%s\t%s\t%d\t%d\t%d\t%d\n' % (n1,n2,val,a,b,c,d))
	out.close()
	print('  wrote outfile to %s' % (outfile))
	sys.stdout.flush()
	return

def get_node_memberships(H):
	nodes = set() ## get proteins and complex members.
	node_membership = {}
	num_complexes = 0
	num_entitysets = 0
	for n in H.get_node_set():
		attrs = H.get_node_attributes(n)		
		if attrs['is_hypernode']:
			nodes.update(attrs['hypernode_members'])
			for m in attrs['hypernode_members']:
				if m not in node_membership:
					node_membership[m] = set()
				node_membership[m].add(n)
			num_complexes+=1
		if attrs['is_entityset']:
			nodes.update(attrs['entityset_members'])
			for m in attrs['entityset_members']:
				if m not in node_membership:
					node_membership[m] = set()
				node_membership[m].add(n)
			num_entitysets+=1
		nodes.add(n)
		if n not in node_membership:
			node_membership[n] = set([n])
	print('%d complexes and %d entity sets' % (num_complexes,num_entitysets))
	print('%d nodes including hypernode and entity set members' % (len(nodes)))	
	return nodes,node_membership

def get_STRING_interactions_in_reactome(interactions,uniprot2pc,infix,nodes,opts):
	interactions_in_reactome = []
	mismapped = 0
	notinreactome = 0
	missing = {}
	for n1,n2,val in interactions:
		if n1 in uniprot2pc and n2 in uniprot2pc:
			un1 = uniprot2pc[n1]
			un2 = uniprot2pc[n2]
		else:
			if n1 not in uniprot2pc:
				missing[n1] = ('NA','NotInPC')
			if n2 not in uniprot2pc:
				missing[n2] = ('NA','NotInPC')
			mismapped+=1
			continue
		
		if un1 in nodes and un2 in nodes:
			interactions_in_reactome.append([un1,un2,val])
		else:
			if un1 not in nodes:
				missing[n1] = (un1,'NotInHypergraph')
			if un2 not in nodes:
				missing[n2] = (un2,'NotInHypergraph')
			notinreactome+=1

	print('  %d INTERACTIONS HAVE BOTH NODES IN REACTOME\n  %d interactions not in PathwayCommons Reactome mapping\n  %d interactions are not in this hypergraph' % (len(interactions_in_reactome),mismapped,notinreactome))
	outfile = make_outfile(opts,OUT_TXT_CHANNEL_DIR,infix+'-mismapped')
	out = open(outfile,'w')
	out.write('#UniProtID\tPathwayCommonsID\tMismappingReason\n')
	for m in missing:
		out.write('%s\t%s\t%s\n' % (m,missing[m][0],missing[m][1]))
	out.close()
	print('  wrote %d (%.4f) mismapped nodes to %s' % (len(missing),len(missing)/len(interactions),outfile))
	sys.stdout.flush()
	return interactions_in_reactome

def get_STRING_channel_interactions(f):
	print('FILE %s' % (f))
	name = f.replace(PROCESSED_STRING_DIR,'').replace('.txt','')
	print('NAME %s' % (name))
	interactions = []
	missing = {}
	with open(f) as fin:
		for line in fin:
			row = line.strip().split()
			interactions.append([row[2],row[3],int(row[4])])
	print('  %d INTERACTIONS' % (len(interactions)))
	return name,interactions

#############################
## VIZ FUNCTIONS
#############################

def viz_histograms(sif_graph_file,compound_graph_file,bipartite_graph_file,hypergraph_file,brelax_file,opts):

	cumulative_file = make_outfile(opts,OUT_VIZ_DIR,'cumulative-histogram',filetype='')
	cumulative_histogram.cumulative_histogram(sif_graph_file,bipartite_graph_file,bipartite_graph_file,hypergraph_file,cumulative_file)

	heatmap_file_unnorm = make_outfile(opts,OUT_VIZ_DIR,'cumulative_heatmap',filetype='')
	heatmap_file_norm = make_outfile(opts,OUT_VIZ_DIR,'cumulative_heatmap_normalized',filetype='')
	heatmap_viz.three_panel(sif_graph_file,compound_graph_file,bipartite_graph_file,hypergraph_file,heatmap_file_unnorm)
	heatmap_viz.three_panel_percentage(sif_graph_file,compound_graph_file,bipartite_graph_file,hypergraph_file,heatmap_file_norm)

	# make singles
	infiles = [sif_graph_file,compound_graph_file,bipartite_graph_file,hypergraph_file,brelax_file]
	titles = ['Graph Connectivity','Compound Graph Connectivity','Bipartite Graph Connectivity','Hypergraph B-Connectivity','Hypergraph B-Relaxation Distance']
	outprefixes = [make_outfile(opts,OUT_VIZ_DIR,'sif_graph_heatmap',filetype=''),\
		make_outfile(opts,OUT_VIZ_DIR,'compound_graph_heatmap',filetype=''), \
		make_outfile(opts,OUT_VIZ_DIR,'bipartite_graph_heatmap',filetype=''), \
		make_outfile(opts,OUT_VIZ_DIR,'hypergraph_heatmap',filetype=''), \
		make_outfile(opts,OUT_VIZ_DIR,'brelax_heatmap',filetype='')]
	for i in range(len(infiles)):
		heatmap_viz.single_panel(infiles[i],titles[i],outprefixes[i],norm=False)
		heatmap_viz.single_panel(infiles[i],titles[i],outprefixes[i]+'_normalized',norm=True)
	return

def viz_permutations(scores_file,perm_infix,num_perms,opts,k_vals=[0,1,2,3,4,5],viz_jaccard=True):
	# read scores file
	scores = read_influence_scores_file(scores_file)

	## visualize jaccard overlap
	M = [0]*len(PATHWAYS)
	for i in range(len(PATHWAYS)):
		M[i] = [0]*len(PATHWAYS)
		for j in range(len(PATHWAYS)):
			M[i][j] = scores[PATHWAYS[i]][PATHWAYS[j]][-1]
	fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6,6))
	ca = ax.matshow(M, aspect='auto', vmin=0.0, vmax=1,cmap=plt.get_cmap('Blues'))
	fig.colorbar(ca,ax=ax)
	ax.set_title('%s Pathway Overlap' % (perm_infix))
	ax.xaxis.set_ticks_position('bottom')
	ax.set_xticks(range(len(PATHWAYS)))
	ax.set_yticks(range(len(PATHWAYS)))
	ax.set_xticklabels([PATHWAY_NAMES[p] for p in PATHWAYS], rotation=270, fontsize=9)
	ax.set_yticklabels([PATHWAY_NAMES[p] for p in PATHWAYS], fontsize=9)

	plt.tight_layout()
	outprefix = make_outfile(opts,OUT_VIZ_DIR,perm_infix+'_jaccard',filetype='')
	plt.savefig(outprefix+'.png')
	print('saved to '+outprefix+'.png')
	plt.savefig(outprefix+'.pdf')
	os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
	print('saved to '+outprefix+'.pdf')

	## get permutation scores
	X = {}
	for k in k_vals:
		X[k] = [0]*len(PATHWAYS)
		for i in range(len(PATHWAYS)):
			X[k][i] = [0]*len(PATHWAYS)
	for p in range(num_perms):
		perm_outfile = make_outfile(opts,OUT_PATHWAY_DIR,'%s_%d_perms_10000_swaps' % (perm_infix,p))	
		perm_scores = read_influence_scores_file(perm_outfile)
		for i in range(len(PATHWAYS)):
			for j in range(len(PATHWAYS)):
				for k in k_vals:
					if k not in perm_scores[PATHWAYS[i]][PATHWAYS[j]]:
						print('PERM_OUTFILE_ERROR:',p,perm_outfile)
						sys.exit()
					if k not in scores[PATHWAYS[i]][PATHWAYS[j]]:
						print('ORIG_SCORE ERROR:',p.scores_file)
					if perm_scores[PATHWAYS[i]][PATHWAYS[j]][k] >= scores[PATHWAYS[i]][PATHWAYS[j]][k]:
						X[k][i][j] += 1/num_perms

	## visualize significant scores for selected values of k
	for k in k_vals:
		## make grid of scatter plots.
		x = []
		y = []
		areas = []
		colors = []
		for i in range(len(PATHWAYS)):
			for j in range(len(PATHWAYS)):
				## here's the thing: influence score is (pathwayA,pathwayB,k,score)
				## pathwayA is on the ROWS and pathwayB is on the COLUMNS
				## so pathwayA is the Y-AXIS and pathwayB is the X-AXIS
				x.append(j)
				y.append(i)
				areas.append(permutation_viz.get_area(X[k][i][j]))
				colors.append(scores[PATHWAYS[i]][PATHWAYS[j]][k])
		fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6.1,6))
		plt.gca().invert_yaxis()
		plt.gca().invert_xaxis()
		ca = ax.scatter(x,y,s=areas,c=colors,cmap=cm.get_cmap('Blues'),vmin=0,vmax=1.0, edgecolors='k',linewidths=0.1)
		ax.set_xlim(-0.5,len(PATHWAYS)-0.5)
		ax.set_ylim(len(PATHWAYS)-0.5,-0.5)
		ax.set_xticks(range(len(PATHWAYS)))
		ax.set_yticks(range(len(PATHWAYS)))
		ax.set_xticklabels([PATHWAY_NAMES[p] for p in PATHWAYS], rotation=270, fontsize=9)
		ax.set_yticklabels([PATHWAY_NAMES[p] for p in PATHWAYS], fontsize=9)
		fig.colorbar(ca, ax=ax, fraction=0.1, aspect=30)
		plt.tight_layout()
		outprefix = make_outfile(opts,OUT_VIZ_DIR,perm_infix+'_k_%d' % (k),filetype='')
		plt.savefig(outprefix+'.png')
		print('saved to '+outprefix+'.png')
		plt.savefig(outprefix+'.pdf')
		os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
		print('saved to '+outprefix+'.pdf')
		plt.close()
	return
		
def read_influence_scores_file(scores_file):
	scores = {}  # pathwayA pathwayB k: float
	with open(scores_file) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split()
			if row[0] not in scores:
				scores[row[0]] = {}
			if row[1] not in scores[row[0]]:
				scores[row[0]][row[1]] = {}
			for k in range(2,len(row)):
				scores[row[0]][row[1]][k-3] = float(row[k]) # pathwayA pathwayB k: float
	return scores

def viz_influence_histogram(H,id2identifier,identifier2id,pathway_members,pathways_to_highlight,opts,buffs):
	COLORS = {'Signaling-by-MET':'#FA7171',
	'Signaling-by-MST1':'#95A5D5',
	'PI3K-AKT-Signaling':'#A8381A',
	'Signaling-by-BMP':'#0008AF',
	'Signaling-by-ERBB4':'#AF0060',
	'ERK1-ERK2-pathway':'#47AF00',
	'Integrin-signaling':'#A27A9B',
	'Signaling-by-Activin':'#FA7171',
	'Signaling-by-TGF-beta-Receptor-Complex':'#95A5D5'
	}

	for p in range(len(pathways_to_highlight)):
		pathway = pathways_to_highlight[p]
		buff = buffs[p]
		suffix = 'case_study_%s' % (pathway)
		print('  PATHWAY %s' % (pathway))

		brelax_sets = survey_hgraph_single_pathway(H,id2identifier,identifier2id,pathway_members[pathway],suffix,opts)
		overlap,running_tot = compute_brelax_overlaps(brelax_sets,pathway_members)
		max_k = 10
		running_tot = running_tot[:max_k+1]
		for p2 in overlap:
			overlap[p2] = overlap[p2][:max_k+1]

		ax = plt.subplot(1,1,1)
		x = range(max_k+1) ## max_k+1 for max_k
		i=0
		text_list = []
		for n in PATHWAYS:
			if n == pathway:
				continue
			num = len(pathway_members[pathway])
			perc = int(overlap[n][-1]/num*10000)/100.0
			if n in pathways_to_highlight:
				ax.plot(x,overlap[n],color=COLORS[n],lw=3,label='_nolegend_',zorder=2)
				text_list.append([n,i,perc,overlap[n][-1]])
				i+=1
			else:
				ax.plot(x,overlap[n],color=[0.8,0.8,0.8],lw=1,label='_nolegend_',zorder=1)

		# adjust text_list to be at least 3 spaces apart
		text_list = space(text_list,buff)	
		for i in range(len(text_list)):
			label = '%s' % (PATHWAY_NAMES[text_list[i][0]])
			ax.text(max_k+.25,text_list[i][3],label,backgroundcolor=COLORS[text_list[i][0]],color='w',fontsize=12)
			ax.plot([max_k,max_k+.25],[overlap[text_list[i][0]][-1],text_list[i][3]],color=COLORS[text_list[i][0]],label='_nolegend_')

		
		ax.plot(range(len(running_tot)),running_tot,'--k',label='_nolegend_')
		print('RUNNING:TOT',running_tot)
		if max_k <= 5:
			ax.text(max_k+.25,running_tot[-1],'All Nodes',backgroundcolor='k',color='w',fontsize=12)
			ax.plot([max_k,max_k+.25],[running_tot[-1],running_tot[-1]],color='k',label='_nolegend_')
		else:
			ax.text(max_k+1,running_tot[-1],'All Nodes',backgroundcolor='k',color='w',fontsize=12)
			ax.plot([max_k,max_k+1],[running_tot[-1],running_tot[-1]],color='k',label='_nolegend_')

		#ax.legend(ncol=2,bbox_to_anchor=(.9, -.15),fontsize=10)
		ax.set_title('Source Pathway %s' % (PATHWAY_NAMES[pathway]),fontsize=18)
		ax.set_xlabel('$k$',fontsize=14)
		ax.set_xticks(range(max_k+1))
		for tick in ax.xaxis.get_major_ticks():
			tick.label.set_fontsize(12)
		for tick in ax.yaxis.get_major_ticks():
			tick.label.set_fontsize(12) 
		#ax.set_xtick_labels(range(max_k+1))
		ax.set_ylabel('# of Nodes',fontsize=14)
		#plt.tight_layout()
		outprefix = make_outfile(opts,OUT_VIZ_DIR,suffix,filetype='')
		plt.savefig(outprefix+'.png')
		print('  saved to '+outprefix+'.png')
		plt.savefig(outprefix+'.pdf')
		os.system('pdfcrop %s.pdf %s.pdf' % (outprefix,outprefix))
		print('  saved to '+outprefix+'.pdf')
		plt.close()
		
	return

def space(text_list,buff):
	text_list.sort(key=lambda x:x[3])
	median=0
	jump = 1

	for i in range(len(text_list)):
		if i < median:
			while text_list[i][3] > text_list[i+1][3]-buff:
				for j in range(i+1):
					if text_list[j][3] > text_list[i+1][3]-(i-j+1)*buff:
						text_list[j][3] = text_list[j][3] - jump
		if i > median:
			while text_list[i][3] < text_list[i-1][3]+buff:
				for j in range(i,len(text_list)):
					if text_list[j][3] < text_list[i-1][3]+(i-j+1)*buff:
						text_list[j][3] = text_list[j][3] + jump
	return text_list

def compute_brelax_overlaps(brelax_sets,pathway_members):
	max_k = max(brelax_sets.keys())
	overlap = {}
	for pathway in PATHWAYS:
		overlap[pathway] = []
		running_total = [] # this will be recomputed over and over but it's oK!
		cumu_p = set()
		for k in range(max_k+1): # skip -1 (init set)
			cumu_p.update(brelax_sets[k])
			running_total.append(len(cumu_p))
			overlap[pathway].append(len(cumu_p.intersection(set(pathway_members[pathway]))))
	return overlap, running_total


#############################
## SURVEY FUNCTIONS
#############################
def survey_graph(G,outfile):

	if not FORCE and os.path.isfile(outfile):
		force_print_statement(outfile)
		return

	print('surveying graph and writing to %s' % (outfile))

	out = open(outfile,'w')
	out.write('#Name\td=1\td=2\td=3\t...\n')
	i = 0
	processed=0
	print('Printing messages every %d nodes' % (int(nx.number_of_nodes(G)/50)))
	for s in G.nodes():
		i+=1
		if i % int(nx.number_of_nodes(G)/50) == 0:
			print(' %d of %d' % (i,nx.number_of_nodes(G)))
		# if this node is a hyperedge_node (from bipartite graph), ignore.
		attrs = nx.get_node_attributes(G,s)
		if 'hedge_node' in G.node[s] and G.node[s]['hedge_node']: 
			continue
		processed+=1
		dist_sets = graph_utils.bfs_histogram(G,s)
		max_dist = max(dist_sets.keys())
		distances = [str(len(dist_sets[d])) for d in range(1,max_dist+1)] # skip d=0
		out.write('%s\t%s\n' % (s,'\t'.join(distances)))
	out.close()
	print('wrote %d processed nodes to %s' % (processed,outfile))
	return

def survey_graph_pathways_threaded(args):
	G,pathway_members,suffix,opts,doubled = args
	return survey_graph_pathways(G,pathway_members,suffix,opts,double_distances=doubled,verbose=False)

def survey_graph_pathways(G,pathway_members,suffix,opts,double_distances=True,verbose=False):
	if verbose:
		print('%d pathways' % (len(pathway_members)))
	start = time.time()

	outfile = make_outfile(opts,OUT_PATHWAY_DIR,suffix)
	hist_dicts = {}
	max_dist = 0
	if FORCE or not os.path.isfile(outfile):
		for pathway in pathway_members:
			if verbose:
				print('  Pathway "%s"' % (pathway))
		
			dist_dict = graph_connectivity_set(G,pathway_members[pathway])
			hist_dicts[pathway] = graph_utils.dist2hist(dist_dict)
			if double_distances: # for bipartite graph
				doubled = {}
				for d in hist_dicts[pathway]:
					doubled[d*2] = hist_dicts[pathway][d]
				hist_dicts[pathway] = doubled
			max_dist = max(max_dist,max(hist_dicts[pathway].keys()))
			#if verbose:
			#	print('max dist:',max_dist)
			#	print(hist_dicts[pathway].keys())
		out = open(outfile,'w')
		out.write('#PathwayA\tPathwayB\t%s\n' % ('\t'.join(['s_k=%d' % (i) for i in range(-1,max_dist+1)])))
		for pathwayA in pathway_members:
			initA = set(pathway_members[pathwayA])
			for pathwayB in pathway_members:
				initB = set(pathway_members[pathwayB])
				if pathwayA == pathwayB:
					jaccard_val = 1
					influence_scores = [0]*(max_dist+1)
				else:
					jaccard_val = asymmetric_jaccard(initA,initB)
					cumulativeA = set(pathway_members[pathwayA])
					cumulativeB = set(pathway_members[pathwayB])
					influence_scores = [0]*(max_dist+1)
					for i in range(max_dist+1):
						cumulativeA.update(hist_dicts[pathwayA].get(i,set()))
						cumulativeB.update(hist_dicts[pathwayB].get(i,set()))
						influence_scores[i] = influence_score_fast(initA,initB,cumulativeA,cumulativeB)
				#print(influence_scores,max_dist)
				out.write('%s\t%s\t%.4f\t%s\n' % (pathwayA,pathwayB,jaccard_val,'\t'.join(['%.4f' % s for s in influence_scores])))
		out.close()
		end = time.time()
		print('%s time: %f' % (outfile,end-start))
	else:
		force_print_statement(outfile)
	return

def graph_connectivity_set(G,pathway_members):
	## fast BFS
	ss = 'ss'
	for node in pathway_members:
		G.add_edge(ss,node)
	this_dist_dict = graph_utils.bfs(G,ss)
	## decrement all distances by one (since we started at ss)
	dist_dict = {}
	for d in this_dist_dict:
		if d == ss:
			continue
		dist_dict[d]=this_dist_dict[d]-1
		
	for node in pathway_members:
		G.remove_edge(ss,node)
	G.remove_node(ss)
	return dist_dict

def survey_hgraph(H,outfile):
	## TODO cplex doesn't work on home laptop
	if not FORCE and os.path.isfile(outfile):
		force_print_statement(outfile)
		return

	out = open(outfile,'w')
	out.write('#Name\td=1\td=2\td=3\t...\n')
	i = 0
	print('Printing messages every %d nodes' % (int(stats.number_of_nodes(H)/50)))
	for node in H.node_iterator():
		i+=1
		if i % int(stats.number_of_nodes(H)/50) == 0:
			print('node %d of %d' % (i,stats.number_of_nodes(H)))

		bconn, ignore, ignore, ignore = hpaths.b_visit(H,node)
		if len(bconn) == 1: ## only one node; skip
			out.write('%s\t1\n' % (node))
			continue

		G = H.get_induced_subhypergraph(bconn)

		dist_dict = {}
		for t in bconn:		
			# get distance from node to t using ILP
			dist_dict[t] = shortest_hyperpath.runILP(G,node,t)

		## convert distances to cumulative list.
		hist_dict = graph_utils.dist2hist(dist_dict,counts=True)
		max_val = int(max(hist_dict.keys()))
		distlist = [hist_dict.get(i,0) for i in range(max_val+1)]
		cumulist = [sum(distlist[:i+1]) for i in range(max_val+1)]
		out.write('%s\t%s\n' % (node,'\t'.join([str(i) for i in cumulist])))
	out.close()
	print('wrote to %s' % (outfile))
	return

def survey_hedges(H,id2identifier,outfile):
	print('Generating hyperedge information file.')
	## writes a file of bconn,traversed,restrictive for each hyperedge
	out1 = open(outfile,'w')
	out1.write('#HyperedgeIdentifier\tNumBVisit\tBVisitNodes\tTraversedHedges\tRestrictiveHedges\n')
	i = 0
	for hedge_id in H.hyperedge_id_iterator():
		identifier = H.get_hyperedge_attribute(hedge_id,'identifier')
		i+=1
		if i % 500 == 0:
			print('hyperedge %d of %d' % (i,stats.number_of_hyperedges(H)))
		bconn, traversed, restrictive = hpaths.b_visit_restrictive(H,H.get_hyperedge_head(hedge_id))
		out1.write('%s\t%d\t%s\t%s\t%s\n' % \
			(identifier,len(bconn),
				';'.join([n for n in bconn]),
				';'.join([id2identifier[e] for e in traversed]),
				';'.join([id2identifier[e] for e in restrictive])))
	out1.close()
	print('wrote to %s' % (outfile))
	return

def survey_hgraph_brelax(H,outfile,id2identifier,identifier2id,opts):
	## TODO cplex doesn't work on home laptop
	if not FORCE and os.path.isfile(outfile):
		force_print_statement(outfile)
		return
	
	out = open(outfile,'w')
	out.write('#Name\td=1\td=2\td=3\t...\n')
	i = 0
	print('Printing messages every %d nodes' % (int(stats.number_of_nodes(H)/50)))
	for node in H.node_iterator():
		i+=1
		if i % int(stats.number_of_nodes(H)/50) == 0:
			print('node %d of %d' % (i,stats.number_of_nodes(H)))

		b_visit_dict = get_bvisit_dict(H,id2identifier,identifier2id,opts)
		dist_dict,ignore = hpaths.b_relaxation(H,set([node]),b_visit_dict=b_visit_dict)

		## convert distances to cumulative list.
		hist_dict = graph_utils.dist2hist(dist_dict,counts=True)
		max_val = int(max(hist_dict.keys()))
		distlist = [hist_dict.get(i,0) for i in range(max_val+1)]
		cumulist = [sum(distlist[:i+1]) for i in range(max_val+1)]
		out.write('%s\t%s\n' % (node,'\t'.join([str(i) for i in cumulist])))
	out.close()
	print('wrote to %s' % (outfile))
	return

def get_bvisit_dict(H,id2identifier,identifier2id,opts):
	## get b-visit dictionary for every hyperedge. Generate file if not already present.
	hedge_connectivity_file = make_outfile(opts,OUT_TXT_DIR,'hyperedge-connectivity-info')
	if FORCE or not os.path.isfile(hedge_connectivity_file):
		survey_hedges(H,id2identifier,hedge_connectivity_file)
	b_visit_dict = hgraph_utils.make_b_visit_dict(hedge_connectivity_file,identifier2id)
	return b_visit_dict

def survey_hgraph_pathways_threaded(args):
	H,id2identifier,identifier2id,pathway_members,suffix,opts = args
	return survey_hgraph_pathways(H,id2identifier,identifier2id,pathway_members,suffix,opts,verbose=False)

def survey_hgraph_pathways(H,id2identifier,identifier2id,pathway_members,suffix,opts,verbose=False):
	if verbose:
		print('%d pathways' % (len(pathway_members)))
	start = time.time()

	outfile = make_outfile(opts,OUT_PATHWAY_DIR,suffix)
	hist_dicts = {}
	max_dist = 0
	if FORCE or not os.path.isfile(outfile):
		b_visit_dict = get_bvisit_dict(H,id2identifier,identifier2id,opts)
		for pathway in pathway_members:
			if verbose:
				print('  Pathway "%s"' % (pathway))
		
			dist_dict,ignore = hpaths.b_relaxation(H,pathway_members[pathway],b_visit_dict=b_visit_dict)
			hist_dicts[pathway] = graph_utils.dist2hist(dist_dict)
			max_dist = max(max_dist,max(hist_dicts[pathway].keys()))
		out = open(outfile,'w')
		out.write('#PathwayA\tPathwayB\t%s\n' % ('\t'.join(['s_k=%d' % (i) for i in range(-1,max_dist+1)])))
		for pathwayA in pathway_members:
			initA = set(pathway_members[pathwayA])
			for pathwayB in pathway_members:
				initB = set(pathway_members[pathwayB])
				if pathwayA == pathwayB:
					jaccard_val = 1
					influence_scores = [0]*(max_dist+1)
				else:
					jaccard_val = asymmetric_jaccard(initA,initB)
					cumulativeA = set(pathway_members[pathwayA])
					cumulativeB = set(pathway_members[pathwayB])
					influence_scores = [0]*(max_dist+1)
					for i in range(max_dist+1):
						cumulativeA.update(hist_dicts[pathwayA].get(i,set()))
						cumulativeB.update(hist_dicts[pathwayB].get(i,set()))
						influence_scores[i] = influence_score_fast(initA,initB,cumulativeA,cumulativeB)
				out.write('%s\t%s\t%.4f\t%s\n' % (pathwayA,pathwayB,jaccard_val,'\t'.join(['%.4f' % s for s in influence_scores])))
		out.close()
		end = time.time()
		print('%s time: %f' % (outfile,end-start))
	else:
		force_print_statement(outfile)
	return

def survey_hgraph_single_pathway(H,id2identifier,identifier2id,members,suffix,opts):
	outfile = make_outfile(opts,OUT_PATHWAY_DIR,suffix)
	to_return = {}
	if FORCE or not os.path.isfile(outfile):
		b_visit_dict = get_bvisit_dict(H,id2identifier,identifier2id,opts)
		dist_dict,ignore = hpaths.b_relaxation(H,members,b_visit_dict=b_visit_dict)
		hist_dicts = graph_utils.dist2hist(dist_dict)
		max_dist = max(hist_dicts.keys())
		out = open(outfile,'w')
		out.write('#k\tNumInSet\tMembers\n')
		out.write('-1\t%d\t%s\n'% (len(members),';'.join(members)))
		to_return[-1] = set(members)
		for k in range(max_dist+1):
			if k in hist_dicts:
				out.write('%d\t%d\t%s\n' % (k,len(hist_dicts[k]),';'.join(hist_dicts[k])))
				to_return[k] = set(hist_dicts[k])
			else:
				out.write('%d\t0\t\n' % (k))
				to_return[k] = set()
		out.close()
	else:
		force_print_statement(outfile)
		with open(outfile) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				row = line.strip().split()
				to_return[int(row[0])] = set(row[2].split(';'))
	return to_return

#############################
## PATHWAY & PERMUTATION FUNCTIONS
#############################
def make_pathways_from_hypergraph(H,outfile):
	out = open(outfile,'w')
	out.write('#Pathway\tNumMembers\tMembers\n')
	pathway_members = {}
	hgraph_nodes = H.get_node_set()
	print(PATHWAYS)
	for name in PATHWAYS:
		members = set()
		with open('%s/%s-hypernodes.txt' % (PATHWAY_HGRAPH_PREFIX,name)) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				members.add(line.strip().split()[0])
		print(name,'intersection with hypergraph:',len(members),'-->',len(members.intersection(hgraph_nodes)))
		members = members.intersection(hgraph_nodes)
		pathway_members[name] = members
		out.write('%s\t%d\t%s\n' % (name,len(members),';'.join(members)))
	out.close()
	print('wrote to %s' % (outfile))
	return

def make_pathways_from_sif_graph(sif_graph,hypergraph_pathway_file,outfile):
	pcid2hgnc, hgnc2pcid = hgraph_utils.get_id_map(PATHWAY_HGRAPH_ENTITIES,common_name=True)
	orig_pathways = read_pathway_members(hypergraph_pathway_file)

	orig_pathway_names = {p:set() for p in orig_pathways.keys()}
	for pathway in orig_pathways:
		num_found = 0
		num_skipped = 0
		for item in orig_pathways[pathway]:
			if item in pcid2hgnc:
				num_found+=1
				orig_pathway_names[pathway].add(pcid2hgnc[item])
			else:
				num_skipped+=1
		print('%s: %d found and %d skipped (skipped are entity sets or complexes, the members of which are recorded)' % (pathway,num_found,num_skipped))

	sif_nodes = sif_graph.nodes()
	#print('%d total SIF nodes' % (len(sif_nodes)))
	#print(list(sif_nodes)[:20])

	out = open(outfile,'w')
	out.write('#Pathway	NumMembers	Members\n')
	pathway_members = {}
	for name in orig_pathway_names:
		members = orig_pathway_names[name].intersection(sif_nodes)
		print('%s: %d -> %d' % (name,len(orig_pathway_names[name]),len(members)))
		if len(members) == 0:
			print(orig_pathway_names[name])
			print(members)
			sys.exit()
		out.write('%s\t%d\t%s\n' % (name,len(members),';'.join(members)))
	out.close()
	print('wrote to %s' % (outfile))
	return

def generate_pathway_permutations(pathways,num_swaps,opts,infix):
	G,set_names,node_membership = permutation_test.generate_graph(pathways)

	for permutation in range(opts.perm_test):
		fname = make_outfile(opts,OUT_PERM_DIR,'%s_%d_perms_%d_swaps' % (infix,permutation,num_swaps))
		if FORCE or not os.path.isfile(fname):
			print('__PERMUTATION__%d__of__%d__' % (permutation+1,opts.perm_test))
			permutation_test.run_permutation(G,pathways,set_names,fname,num_swaps)
		else:
			force_print_statement(fname)
	print('done with %d permutations' % (permutation))

def influence_score_fast(initp1,initp2,cumulative1,cumulative2):
	init_intersection = initp1.intersection(initp2)	
	numerator = len(cumulative1.intersection(initp2).difference(init_intersection))
	denominator = len(cumulative1.difference(init_intersection))
	score = numerator/denominator
	if score < 0 or score > 1:
		print(score)
		sys.exit()
	return score

def asymmetric_jaccard(initp1,initp2):
	jaccard = len(initp1.intersection(initp2))/len(initp1)
	return jaccard

def get_pathways_for_STRING_analysis(opts):

	# these are the "top lvel" reactome pathways - they are too general. Ignore.
	TO_IGNORE = ['Circadian-Clock', 'Cell-Cycle', 'Disease',  'Programmed-Cell-Death',  'Extracellular-matrix-organization',  'Vesicle-mediated-transport', \
	 'Cellular-responses-to-external-stimuli',  'Organelle-biogenesis-and-maintenance',  'Neuronal-System',  'NICD-traffics-to-nucleus',  'Signaling-Pathways',  \
	 'Metabolism-of-RNA',  'DNA-Repair',  'Metabolism',  'Mitophagy',  'Gene-expression-(Transcription)',  'Developmental-Biology',  'Chromatin-organization', \
	 'Transport-of-small-molecules',  'Immune-System',  'Metabolism-of-proteins',  'Muscle-contraction',  'Digestion-and-absorption',  'Reproduction', \
	  'Hemostasis',  'Cell-Cell-communication']

	pathway_nodes = {}
	outfile = make_outfile(opts,OUT_TXT_DIR,'pathways-for-STRING')
	if FORCE or not os.path.isfile(outfile):
		files = glob.glob('%s/*-hypernodes.txt' % (PATHWAY_HGRAPH_PREFIX))
		print('%d files' % (len(files)))
		for f in files:
			name = f.replace(PATHWAY_HGRAPH_PREFIX,'').replace('-hypernodes.txt','')
			#print(name)
			pathway_nodes[name] = set()
			with open(f) as fin:
				for line in fin:
					if line[0] == '#':
						continue
					row = line.strip().split()
					#print(row)
					pathway_nodes[name].add(row[0])
					if len(row) > 1: ## to be safe, add all hypernodes AND members. Add everything.
						pathway_nodes[name].update(row[1].split(';'))
		print('%d pathways in total' % (len(pathway_nodes)))

		for i in TO_IGNORE:
			del pathway_nodes[i]
		print('%d pathways after removing top-level pathways' % (len(pathway_nodes)))

		#to remove redundants
		to_remove = set()
		pathway_list = list(pathway_nodes.keys())
		for i in pathway_list:
			for j in pathway_list:
				if i != j and len(pathway_nodes[i].intersection(pathway_nodes[j])) == len(pathway_nodes[i]):
					to_remove.add(i)
					break
		print('Removing %d redundant sets' % (len(to_remove)))
		for t in to_remove:
			del pathway_nodes[t]


		print('%d pathways remain' % (len(pathway_nodes)))
		out = open(outfile,'w')
		out.write('#Pathway\tNumMembers\tMembers\n')
		for p in pathway_nodes:
			out.write('%s\t%d\t%s\n' % (p,len(pathway_nodes[p]),';'.join(pathway_nodes[p])))
		out.close()
		print('Wrote to %s' % outfile)

	else:
		force_print_statement(outfile)
		with open(outfile) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				row = line.strip().split()
				pathway_nodes[row[0]] = set(row[2].split(';'))
	
	all_pathway_nodes = set()
	for p in pathway_nodes:
		all_pathway_nodes.update(pathway_nodes[p])

	return pathway_nodes,all_pathway_nodes

#############################
## DATA READERS AND UTILITY FUNCTIONS
#############################
def get_representations(small_molecule_filter,blacklist_filter,keep_singletons):
	## SIF graph is an edges file; by definition there are no singletons.
	sif_graph = graph_utils.read_graph(SIF_FILE,SIF_CONV_FILE)
	nodes = sif_graph.nodes()
	
	pcid2hgnc, hgnc2pcid = hgraph_utils.get_id_map(PATHWAY_HGRAPH_ENTITIES,common_name=True)
	nodes_to_keep = set()

	if small_molecule_filter:
		ADDITIONAL_MOLECULES = set([
	'http://pathwaycommons.org/pc2/Protein_fa92e525bc3eb51cffd7786a1f19e317', # Ub
	'http://pathwaycommons.org/pc2/Complex_6ae49f2fe344df4b6985f2f372910a77', # Nuclear Pore Complex
	'http://pathwaycommons.org/pc2/Protein_80c9e4746b9a9261c9c7b174d2cf8292', # Ub
	])
		for n in nodes:
			if n in hgnc2pcid and hgnc2pcid[n] not in ADDITIONAL_MOLECULES and 'SmallMolecule' not in hgnc2pcid[n]:
				nodes_to_keep.add(n)
			if n in hgnc2pcid and (hgnc2pcid[n] in ADDITIONAL_MOLECULES or 'SmallMolecule' in hgnc2pcid[n]):
				print('skipping small mol node',n)
	elif blacklist_filter:
		blacklisted = set()
		with open(BLACKLIST_FILE) as fin:
			for line in fin: # take first column of blacklist file
				blacklisted.add(line.strip().split()[0])
		for n in nodes:
			if n in hgnc2pcid and hgnc2pcid[n] not in blacklisted:
				nodes_to_keep.add(n)
			if n in hgnc2pcid and hgnc2pcid[n] in blacklisted:
				print('skippig blacklisted node',n)
	else: # make sure nodes are in pcid
		nodes_to_keep = set([n for n in nodes if n in hgnc2pcid.keys()])	
	sif_graph = sif_graph.subgraph(nodes_to_keep)
	print('filtering to keep %d nodes' % (len(nodes_to_keep)))
	print('Graph after filtering by PCID: %d nodes' % (len(sif_graph.nodes())))
	
	# compound graph later
	if small_molecule_filter:
		hypergraph, identifier2id, id2identifier = hgraph_utils.make_hypergraph(HGRAPH_SMALLMOL_PREFIX,keep_singleton_nodes=keep_singletons)
	elif blacklist_filter:
		hypergraph, identifier2id, id2identifier = hgraph_utils.make_hypergraph(HGRAPH_BLACKLIST_PREFIX,keep_singleton_nodes=keep_singletons)
	else:
		hypergraph, identifier2id, id2identifier = hgraph_utils.make_hypergraph(HGRAPH_PREFIX,keep_singleton_nodes=keep_singletons)
	bipartite_graph = hgraph_utils.to_bipartite_graph(hypergraph)
	compact_bipartite_graph = hgraph_utils.to_digraph(hypergraph)
	return sif_graph,None,bipartite_graph,compact_bipartite_graph,hypergraph,identifier2id,id2identifier

def read_pathway_members(pathwayfile):
	pathway_members = {}
	with open(pathwayfile) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split()
			pathway_members[row[0]] = row[2].split(';')
			assert len(pathway_members[row[0]]) == int(row[1])
	return pathway_members

def force_print_statement(fname):
	print("File %s exists. Use --force or remove file to run." % (fname))
	return

def make_outfile(opts,dir_prefix,infix,filetype='.txt'):
	outfile = dir_prefix+infix
	if opts.keep_singletons:
		outfile+='_with-singletons'
	if opts.small_molecule_filter:
		outfile+='_small-molecule-filter'
	if opts.blacklist_filter:
		outfile+='_blacklist-filter'
	outfile += filetype
	return outfile

def compress_and_remove_files(tarfile,regex):
	cmd = 'tar -cvzf %s %s' % (tarfile,regex)
	print(cmd)
	os.system(cmd)
	cmd = 'rm -f %s' % (regex)
	print(cmd)
	os.system(cmd)
	return


#############################
## OPTION PARSER
#############################
def parse_options():
	parser = argparse.ArgumentParser(description='Run experiments for pathway connectivity. The four represenations are "SIF-Graph","Compound Graph","Bipartite Graph", and "Hypergraph"')
	
	## general arguments
	group1 = parser.add_argument_group('General Arguments')
	group1.add_argument('--force', action='store_true', default=False,
		help='force existing files to be overwritten (default=False)')
	group1.add_argument('--printonly',action='store_true',default=False,
		help='print the commands instead of running them (default=False)')

	## datasets
	group2 = parser.add_argument_group('Dataset Arguments')
	group2.add_argument('--keep_singletons',action='store_true',default=False,
		help='Keep singleton nodes. Default False.')
	group2.add_argument('--small_molecule_filter',action='store_true',default=False,
		help='Filter by small molecule nodes')
	group2.add_argument('--blacklist_filter',action='store_true',default=False,
		help='Filter by PathwayCommons blacklist filtered nodes')
	#group2.add_argument('--common_nodes',action='store_true',default=False,
	#	help='Filter representations so they have a similar node set (filter graph).')

	## experiments
	group3 = parser.add_argument_group('Experiment Arguments')
	group3.add_argument('--stats', action='store_true',default=False,
		help='print statistics about each representation.')
	group3.add_argument('--histograms', action='store_true',default=False,
		help='print histograms and heatmaps for each representation.')
	group3.add_argument('--perm_test',metavar='#PERMS',type=int,default=None,
		help='pathway influence permutation test with #PERMS number of permutations.')
	group3.add_argument('--set_seed', action='store_true',default=False,
		help='Set the seed of the random number generator to 123456. Default False.')
	group3.add_argument('--case_studies',  action='store_true',default=False,
		help='pathway influence case studies (hard-coded)')
	group3.add_argument('--string_channels',action='store_true',default=False,
		help='run STRING channel assessment.')
	## TODO hub survey

	opts = parser.parse_args()

	## run checks
	if opts.small_molecule_filter and opts.blacklist_filter:
		sys.exit('ERROR: cannot filter by both small molecules and the blacklisted set.')

	if not opts.perm_test and opts.set_seed:
		sys.exit('ERROR: setting the seed is only useful when --perm_test is specified.')

	## set global variables
	global FORCE,PRINT_ONLY
	if opts.force:
		FORCE = True
	if opts.printonly:
		PRINT_ONLY = True

	if opts.set_seed:
		random.seed(123456)

	return opts

if __name__ == '__main__':
	main()