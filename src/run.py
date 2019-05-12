from __future__ import print_function

## packages
import argparse
import sys
import os
import networkx as nx  ## TODO verify version 
import time
import random

## custom code
import SIF.graph_utils as graph_utils
import hypergraph_code.hgraph_utils as hgraph_utils
import viz.viz_utils as viz_utils

# TODO cplex is broken for python 3.7
#import hypergraph_code.ILP.shortest_hyperpath as shortest_hyperpath 

## halp
from halp.utilities import directed_statistics as stats
from halp.algorithms import directed_paths as hpaths

## TODO: make these options or a config file
SIF_FILE = 'SIF/PathwayCommons10.reactome.hgnc.sif'
SIF_CONV_FILE = 'SIF/conversion-types.txt'
HGRAPH_PREFIX = '../hypergraph/reactome_hypergraph_full/reactome'
PATHWAY_HGRAPH_PREFIX = '../hypergraph/reactome_hypergraphs_parsed/'

## outdirectories
OUT_TXT_DIR = 'out_txt/'
if not os.path.isdir(OUT_TXT_DIR):
	os.system('mkdir %s' % (OUT_TXT_DIR))
OUT_PERM_DIR = 'out_txt/permutations/'
if not os.path.isdir(OUT_PERM_DIR):
	os.system('mkdir %s' % (OUT_PERM_DIR))
OUT_VIZ_DIR = 'out_viz/'
if not os.path.isdir(OUT_VIZ_DIR):
	os.system('mkdir %s' % (OUT_TXT_DIR))

## global variables (will be set with parse_options)
FORCE = False
PRINT_ONLY = False

PATHWAYS = viz_utils.sorted_pathways 
#############################
## MAIN FUNCTION
#############################
def main():

	opts = parse_options()
	print(opts)
	print()

	## get representations no matter what
	sif_graph,compound_graph,bipartite_graph,hypergraph = \
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
		print('Hypergraph: %d nodes and %d hyperedges' % (stats.number_of_nodes(hypergraph),stats.number_of_hyperedges(hypergraph)))
		print('------ END STATS -------')
		print()

	########## Make survey histograms
	if opts.histograms:
		## make parameterized files if they don't exist (or --force is specified)
		survey_graph(sif_graph,make_outfile(opts,OUT_TXT_DIR,'sif-graph-survey'))
		survey_graph(bipartite_graph,make_outfile(opts,OUT_TXT_DIR,'bipartite-graph-survey'))
		## TODO: this is broken. Copy files from work laptop over before continuing.
		#survey_hgraph(hypergraph,make_outfile(opts,OUT_TXT_DIR,'hypergraph-survey'))
		## TODO: finish this section.

	########## Run permutation tests
	if opts.perm_test:
		print(len(PATHWAYS))
		outfile = make_outfile(opts,OUT_TXT_DIR,'pathways-from-hypergraph')
		if FORCE or not os.path.isfile(outfile):
			make_pathways_from_hypergraph(hypergraph,outfile)
		else:
			force_print_statement(outfile)

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

def survey_hgraph(H,outfile):
	## TODO cplex doesn't work.
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
		hist_dict = graph_utils.dist2hist(dist_dict)
		max_val = int(max(hist_dict.keys()))
		distlist = [hist_dict.get(i,0) for i in range(max_val+1)]
		cumulist = [sum(distlist[:i+1]) for i in range(max_val+1)]
		out.write('%s\t%s\n' % (node,'\t'.join([str(i) for i in cumulist])))
	out.close()
	print('wrote to %s' % (outfile))
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
		print(len(members),'-->',len(members.intersection(hgraph_nodes)))
		members = members.intersection(hgraph_nodes)
		pathway_members[name] = members
		out.write('%s\t%d\t%s\n' % (name,name,len(members),';'.join(members)))
	out.close()
	print('wrote to %s' % (outfile))
	return

#############################
## DATA READERS AND UTILITY FUNCTIONS
#############################
def get_representations(small_molecule_filter,blacklist_filter,keep_singletons):
	## SIF graph is an edges file; by definition there are no singletons.
	sif_graph = graph_utils.read_graph(SIF_FILE,SIF_CONV_FILE)

	# compound graph later
	hypergraph, identifier2id, id2identifier = hgraph_utils.make_hypergraph(HGRAPH_PREFIX,keep_singleton_nodes=keep_singletons)
	bipartite_graph = hgraph_utils.to_bipartite_graph(hypergraph)
	return sif_graph,None,bipartite_graph,hypergraph

def force_print_statement(fname):
	print("Not running graph survey b/c file %s exists. Use --force or remove file to regenerate." % (fname))
	return

def make_outfile(opts,dir_prefix,infix):
	outfile = dir_prefix+infix
	if opts.keep_singletons:
		outfile+'_with-singletons'
	if opts.small_molecule_filter:
		outfile+'_small-molecule-filter'
	if opts.blacklist_filter:
		outfile+'_blacklist-filter'
	outfile += '.txt'
	return outfile


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