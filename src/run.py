from __future__ import print_function

## packages
import argparse
import sys
import os
import networkx as nx  ## TODO verify version 
import time
import random
from multiprocessing import Pool ## https://docs.python.org/3.4/library/multiprocessing.html?highlight=process

## custom code
import SIF.graph_utils as graph_utils
import hypergraph_code.hgraph_utils as hgraph_utils
import viz.viz_utils as viz_utils

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
PATHWAY_HGRAPH_PREFIX = '../hypergraph/reactome_hypergraphs_parsed/'

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
	sif_graph,compound_graph,bipartite_graph,hypergraph,identifier2id,id2identifier = \
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
		survey_hgraph(hypergraph,make_outfile(opts,OUT_TXT_DIR,'hypergraph-survey'))
		survey_hgraph_brelax(hypergraph,make_outfile(opts,OUT_TXT_DIR,'hypergraph-brelax-survey'),id2identifier,identifier2id,opts)

	########## Run permutation tests
	if opts.perm_test:
		####### HYPERGRAPH
		pathwayfile = make_outfile(opts,OUT_TXT_DIR,'pathways-from-hypergraph')
		if FORCE or not os.path.isfile(pathwayfile):
			make_pathways_from_hypergraph(hypergraph,pathwayfile)
		else:
			force_print_statement(pathwayfile)
		pathway_members = read_pathway_members(pathwayfile)

		# calculate original pathway survey first
		survey_hgraph_pathways(hypergraph,id2identifier,identifier2id,pathway_members,'-hypergraph',opts,verbose=True)

		# calculate permutation tests with 10,000 swaps each. 
		num_swaps = 10000 
		generate_pathway_permutations(pathway_members,num_swaps,opts)

		# run permutation surveys
		tarfile = make_outfile(opts,OUT_PATHWAY_DIR,'permutation_runs_%d_perms_%d_swaps-hypergraph' % (opts.perm_test,num_swaps),filetype='.tgz')
		if FORCE or not os.path.isfile(tarfile):
			num_threads = 4
			args_to_run = []
			for permutation in range(opts.perm_test):
				suffix = '_%d_perms_%d_swaps-hypergraph' % (permutation,num_swaps)
				permfile = make_outfile(opts,OUT_PERM_DIR,'hypergraph_%d_perms_%d_swaps' % (permutation,num_swaps))
				perm_pathway_members = read_pathway_members(permfile)
				args_to_run.append((hypergraph,id2identifier,identifier2id,perm_pathway_members,suffix,opts))
			print('Running %d args with %d threads' % (len(args_to_run),num_threads))
			with Pool(num_threads) as p:
				p.map(survey_hgraph_pathways_threaded,args_to_run)
			print('Done')
			## compress files
			regex = make_outfile(opts,OUT_PATHWAY_DIR,'*_perms_%d_swaps-hypergraph' % (num_swaps))
			compress_and_remove_files(tarfile,regex)
		else:
			force_print_statement(tarfile)
		
		

		####### BIPARTITE GRAPH
		# calculate original bipartite pathway survey first
		survey_graph_pathways(bipartite_graph,pathway_members,'-bipartite',opts,verbose=True)

		tarfile = make_outfile(opts,OUT_PATHWAY_DIR,'permutation_runs_%d_perms_%d_swaps-bipartite' % (opts.perm_test,num_swaps),filetype='.tgz')
		if FORCE or not os.path.isfile(tarfile):
			# run permutation surveys
			num_threads = 4
			args_to_run = []
			
			for permutation in range(opts.perm_test):
				suffix = '_%d_perms_%d_swaps-bipartite' % (permutation,num_swaps)
				example_outfile = make_outfile(opts,OUT_PATHWAY_DIR,example_pathway+suffix)
				if FORCE or not os.path.isfile(example_outfile):
					permfile = make_outfile(opts,OUT_PERM_DIR,'hypergraph_%d_perms_%d' % (permutation,num_swaps))
					perm_pathway_members = read_pathway_members(permfile)
					args_to_run.append((bipartite_graph,perm_pathway_members,suffix,opts))
			print('Running %d args with %d threads' % (len(args_to_run),num_threads))
			with Pool(num_threads) as p:
				p.map(survey_graph_pathways_threaded,args_to_run)
			print('Done')
			regex = make_outfile(opts,OUT_PATHWAY_DIR,'*_perms_%d_swaps-bipartite' % (num_swaps))
			compress_and_remove_files(tarfile,regex)
		else:
			force_print_statement(tarfile)

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
	G,pathway_members,suffix,opts = args
	return survey_graph_pathways(G,pathway_members,suffix,opts,verbose=False)

def survey_graph_pathways(G,pathway_members,suffix,opts,verbose=False):
	if verbose:
		print('%d pathways' % (len(pathway_members)))

	for pathway in pathway_members:
		if verbose:
			print('  Pathway "%s"' % (pathway))
		outfile = make_outfile(opts,OUT_PATHWAY_DIR,pathway+suffix)
		if FORCE or not os.path.isfile(outfile):
			dist_dict = graph_connectivity_set(G,pathway_members)
			hist_dict = graph_utils.dist2hist(dist_dict)
			out = open(outfile,'w')
			out.write('#k\t#Members\tMembers\n')
			out.write('-1\t%d\t%s\n' % (len(pathway_members[pathway]),';'.join(pathway_members[pathway]))) # write original set
			for dist in range(max(hist_dict)+1):
				if dist in hist_dict:
					out.write('%d\t%d\t%s\n' % (dist,len(hist_dict[dist]),';'.join(hist_dict[dist])))
				else:
					out.write('%s\t0\t\n' % (dist))
			out.close()
			print('wrote to %s' % (outfile))
		else:
			force_print_statement(outfile)

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

		b_visit_dict = get_bvisit_dict(H,identifier2id,opts)
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

def get_bvisit_dict(H,identifier2id,opts):
	## get b-visit dictionary for every hyperedge. Generate file if not already present.
	hedge_connectivity_file = make_outfile(opts,OUT_TXT_DIR,'hyperedge-connectivity-info')
	if FORCE or not os.path.isfile(hedge_connectivity_file):
		survey_hedges(H,id2identifier,hedge_connectivity_file)
	else:
		force_print_statement(hedge_connectivity_file)
	b_visit_dict = hgraph_utils.make_b_visit_dict(hedge_connectivity_file,identifier2id)
	return b_visit_dict

def survey_hgraph_pathways_threaded(args):
	H,id2identifier,identifier2id,pathway_members,suffix,opts = args
	return survey_hgraph_pathways(H,id2identifier,identifier2id,pathway_members,suffix,opts,verbose=False)

def survey_hgraph_pathways(H,id2identifier,identifier2id,pathway_members,suffix,opts,verbose=False):
	if verbose:
		print('%d pathways' % (len(pathway_members)))
	b_visit_dict = {}
	for pathway in pathway_members:
		if verbose:
			print('  Pathway "%s"' % (pathway))
		outfile = make_outfile(opts,OUT_PATHWAY_DIR,pathway+suffix)
		if FORCE or not os.path.isfile(outfile):
			if len(b_visit_dict) == 0:
				b_visit_dict = get_bvisit_dict(H,identifier2id,opts)
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

def generate_pathway_permutations(pathways,num_swaps,opts,sif_graph=False):
	G,set_names,node_membership = permutation_test.generate_graph(pathways)

	for permutation in range(opts.perm_test):
		if sif_graph:
			fname = make_outfile(opts,OUT_PERM_DIR,'sif-graph_%d_perms_%d_swaps' % (permutation,num_swaps))
		else:
			fname = make_outfile(opts,OUT_PERM_DIR,'hypergraph_%d_perms_%d_swaps' % (permutation,num_swaps))
		if FORCE or not os.path.isfile(fname):
			print('__PERMUTATION__%d__of__%d__' % (permutation+1,opts.perm_test))
			permutation_test.run_permutation(G,pathways,set_names,fname,num_swaps)
		else:
			force_print_statement(fname)
	print('done with %d permutations' % (permutation))

#############################
## DATA READERS AND UTILITY FUNCTIONS
#############################
def get_representations(small_molecule_filter,blacklist_filter,keep_singletons):
	## SIF graph is an edges file; by definition there are no singletons.
	sif_graph = graph_utils.read_graph(SIF_FILE,SIF_CONV_FILE)
	# compound graph later
	hypergraph, identifier2id, id2identifier = hgraph_utils.make_hypergraph(HGRAPH_PREFIX,keep_singleton_nodes=keep_singletons)
	bipartite_graph = hgraph_utils.to_bipartite_graph(hypergraph)
	return sif_graph,None,bipartite_graph,hypergraph,identifier2id,id2identifier

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
	print("File %s exists, so skipping this step. Use --force or remove file to regenerate." % (fname))
	return

def make_outfile(opts,dir_prefix,infix,filetype='.txt'):
	outfile = dir_prefix+infix
	if opts.keep_singletons:
		outfile+'_with-singletons'
	if opts.small_molecule_filter:
		outfile+'_small-molecule-filter'
	if opts.blacklist_filter:
		outfile+'_blacklist-filter'
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