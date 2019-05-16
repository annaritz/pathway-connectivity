from __future__ import print_function

## packages
import argparse
import sys
import os
import networkx as nx  ## TODO verify version 
import time
import random
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

## outdirectories
OUT_TXT_DIR = 'out_txt/'
if not os.path.isdir(OUT_TXT_DIR):
	os.system('mkdir %s' % (OUT_TXT_DIR))
OUT_PERM_DIR = 'out_txt/permutations/'
if not os.path.isdir(OUT_PERM_DIR):
	os.system('mkdir %s' % (OUT_PERM_DIR))
OUT_PATHWAY_DIR = 'out_txt/pathway_survey/'
if not os.path.isdir(OUT_PATHWAY_DIR):
	os.system('mkdir %s' % (OUT_PATHWAY_DIR))
OUT_VIZ_DIR = 'out_viz/'
if not os.path.isdir(OUT_VIZ_DIR):
	os.system('mkdir %s' % (OUT_VIZ_DIR))

## global variables (will be set with parse_options)
FORCE = False
PRINT_ONLY = False

PATHWAYS = viz_utils.sorted_pathways 
PATHWAY_NAMES = viz_utils.NAMES
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


	########## Run permutation tests
	if opts.perm_test:

		## Get pathway files
		pathwayfile = make_outfile(opts,OUT_TXT_DIR,'pathways-from-hypergraph')
		if FORCE or not os.path.isfile(pathwayfile):
			make_pathways_from_hypergraph(hypergraph,pathwayfile)
		else:
			force_print_statement(pathwayfile)
		pathway_members = read_pathway_members(pathwayfile)

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
	return

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

def survey_hgraph_pathways_old(H,id2identifier,identifier2id,pathway_members,suffix,opts,verbose=False):
	if verbose:
		print('%d pathways' % (len(pathway_members)))
	b_visit_dict = {}
	for pathway in pathway_members:
		if verbose:
			print('  Pathway "%s"' % (pathway))
		outfile = make_outfile(opts,OUT_PATHWAY_DIR,pathway+suffix)
		if FORCE or not os.path.isfile(outfile):
			if len(b_visit_dict) == 0:
				b_visit_dict = get_bvisit_dict(H,id2identifier,identifier2id,opts)
			dist_dict,traversed = hpaths.b_relaxation(H,pathway_members[pathway],b_visit_dict=b_visit_dict)
			hist_dict = graph_utils.dist2hist(dist_dict)
			out = open(outfile,'w')
			out.write('#k\t#Members\tMembers\tIncidentNodes\n')
			out.write('-1\t%d\t%s\t%s\n' % (len(pathway_members[pathway]),';'.join(pathway_members[pathway]),'None')) # write original set
			for dist in range(max(hist_dict)+1):
				if dist in hist_dict:
					out.write('%d\t%d\t%s\t%s\n' % (dist,len(hist_dict[dist]),';'.join(hist_dict[dist]),';'.join(traversed[dist])))
				else:
					out.write('%s\t0\t\t\t\n' % (dist))
			out.close()
			print('wrote to %s' % (outfile))
		else:
			force_print_statement(outfile)
	return 

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
	group3.add_argument('--case_study_A', metavar='NAME', type=str, default=None,
		help='pathway influence case study "upstream" pathway')
	group3.add_argument('--case_study_B', metavar='NAME', type=str, default=None,
		help='pathway influence case study "downstream" pathway')
	group3.add_argument('--string_channels',action='store_true',default=False,
		help='run STRING channel assessment.')

	opts = parser.parse_args()

	## run checks
	if (opts.case_study_A and not opts.case_study_B) or (opts.case_study_B and not opts.case_study_A):
		sys.exit('ERROR: case studies requires both --case_study_A and --case_study_B.')

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