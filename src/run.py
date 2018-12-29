from __future__ import print_function

import argparse
import sys
import itertools
from halp import directed_hypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform
import matplotlib.pyplot as plt

def main():

	opts = parse_options()
	if opts.test:
		H,names = make_hypergraph('../hypergraph/Signaling-by-WNT')
		example_source = set(['http://pathwaycommons.org/pc2/Protein_a919000a9c31204e20860513b512c9e4'])
		dist = hpaths.b_relaxation(H,example_source)
		hist = make_hist(dist)

		G = transform.to_graph_decomposition(H)
		dist = hpaths.b_relaxation(G,example_source)
		hist = make_hist(dist)

		G2 = to_graph_decomposition_considering_hypernodes(H,names)
		dist = hpaths.b_relaxation(G2,example_source)
		hist = make_hist(dist)

	file_prefixes = []
	if opts.single:
		print("Running Experiments for single pathway with file prefix " + opts.single)
		file_prefixes = [opts.single]
	## todo: if directory is passed get file prefixes.

	if opts.strict:
		print("Running Experiments to compute histograms of strict connectivity defs.")
		for prefix in file_prefixes:
			H,names = make_hypergraph(prefix)
			G = transform.to_graph_decomposition(H)
			G2 = to_graph_decomposition_considering_hypernodes(H,names)
			num_nodes = stats.number_of_nodes(H)
			hgraph_conn = [0]*num_nodes
			graph_conn = [0]*num_nodes
			i = 0
			out = open('tmp-strict.txt','w')
			out.write('HGraph\tGraph\tName\n')
			for node in H.node_iterator():
				bconn, ignore, ignore, ignore = hpaths.b_visit(H,node)
				hgraph_conn[i] = len(bconn)
				bconn, ignore, ignore, ignore = hpaths.b_visit(G,node)
				graph_conn[i] = len(bconn)
				out.write('%d\t%d\t%s\n' % (hgraph_conn[i],graph_conn[i],H.get_node_attribute(node,'name')))
				i+=1
			out.close()
			print('wrote to tmp-strict.txt')
			
			orig_graph_conn = graph_conn
			orig_hgraph_conn = hgraph_conn

			## introduce SuperSource
			ss = 'SuperSource'
			H.add_node(ss,{'is_hypernode':False})
			num_nodes = stats.number_of_nodes(G2)
			hgraph_conn = [0]*num_nodes
			graph_conn = [0]*num_nodes
			i=0
			out = open('tmp-strict-decompose.txt','w')
			out.write('HGraph\tGraph\tNumContainingNodes\tName\n')
			for node in G2.node_iterator():
				# get nodes in H that contain this node
				containing_nodes = set()
				
				for n in H.node_iterator():
					if n == node or (H.get_node_attribute(n,'is_hypernode') and node in H.get_node_attribute(n,'hypernode_members')):	
						containing_nodes.add(n)
				# introduce super source & edges to containing_nodes
				
				hyperedge_ids = set()
				for n in containing_nodes:
					hyperedge_ids.add(H.add_hyperedge(set([ss]),set([n])))

				bconn, ignore, ignore, ignore = hpaths.b_visit(H,ss)
				hgraph_conn[i] = len(bconn)-1 # ignore SS node.		

				## remove nodes
				for hyperedge_id in hyperedge_ids:
					H.remove_hyperedge(hyperedge_id)

				bconn, ignore, ignore, ignore = hpaths.b_visit(G2,node)
				graph_conn[i] = len(bconn)

				out.write('%d\t%d\t%d\t%s\n' % (hgraph_conn[i],graph_conn[i],len(containing_nodes),G2.get_node_attribute(node,'name')))
				i+=1
			out.close()
			print('wrote to tmp-strict-decompose.txt')
			plot_hist(orig_hgraph_conn,orig_graph_conn,hgraph_conn,graph_conn,'Including Complexes','Breaking Apart Complexes','connectivity-histogram')

		






def parse_options():
	parser = argparse.ArgumentParser(description='Run experiments for pathway connectivity.')
	parser.add_argument('--force', action='store_true', default=False,
		help='force existing files to be overwritten (default=False)')
	parser.add_argument('--printonly',action='store_true',default=False,
		help='print the commands instead of running them (default=False)')
	parser.add_argument('--single', metavar='NAME', type=str,
		help='compute connectivity for a single pathway with prefix NAME.')
	parser.add_argument('--test', action='store_true',default=False,
		help='Run a test example for Wnt signaling.')
	parser.add_argument('--strict', action='store_true',default=False,
		help='Calculate histograms of strict distiances for graph and hypergraph versions.')

	opts = parser.parse_args()

	if opts.single == None and opts.test == None:
		sys.exit('ERROR: must include one of --single or --test.')

	if opts.single and not (opts.strict):
		sys.exit('ERROR: --single is specified but no experiments will be run: need one of --strict')

	return opts

def make_hypergraph(file_prefix,delim=';',sep='\t'):
	hypernodes = {}
	with open(file_prefix+'-hypernodes.txt') as fin:
		for line in fin:
			if line[0] == '#': 
				continue
			row = line.strip().split(sep)
			## TODO fix -- this happens when a complex contains other complexes
			if len(row) == 1:
				hypernodes[row[0]] = ['OtherComplexes-FIX']
			else:
				hypernodes[row[0]] = row[1].split(delim)

	names = {'OtherComplexes-FIX':'OtherComplexes-FIX'}
	files_with_names = [file_prefix+suff for suff in ['-complexes.txt','-entitysets.txt','-elements.txt']]
	for f in files_with_names:
		with open(f) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				row = line.strip().split(sep)
				names[row[0]] = row[1]
	
	H = directed_hypergraph.DirectedHypergraph()	
	with open(file_prefix+'-hyperedges.txt') as fin:
		for line in fin:
			if line[0] == '#':
				continue

			row = line.strip().split(sep)
			tail = set()
			head = set()

			## Tail includes tail and regulators.
			## Head includes head.
			if row[0] != 'None':
				tail.update(row[0].split(delim))
			if row[1] != 'None':
				head.update(row[1].split(delim))
			if row[2] != 'None':
				tail.update(row[2].split(delim))
			if row[3] != 'None':
				tail.update(row[3].split(delim))
			hedge_id = row[4]

			if len(tail) > 0 and len(head) > 0:
				H.add_hyperedge(tail,head,identifier=hedge_id)
			else:
				print('WARNING: skipping hyperedge',hedge_id)

	## annotate nodes
	num_hypernodes = 0
	for node in H.get_node_set():
		if node in hypernodes and hypernodes[node] != [node]:
			H.add_node(node,hypernode_members=hypernodes[node],is_hypernode=True)
			num_hypernodes+=1
		else:
			H.add_node(node,is_hypernode=False)
		if node in names:
			H.add_node(node,name=names[node])
		else:
			print('WARNING: Node has no name:',node)
			H.add_node(node,name='None')

	print('Hypergraph has %d hyperedges and %d nodes (%d of these are hypernodes)' % 
		(stats.number_of_hyperedges(H),stats.number_of_nodes(H),num_hypernodes))

	return H, names

def make_hist(dist,verbose=False):
	sorted_vals = sorted(dist.values())
	num_disconnected = sum([v == None for v in sorted_vals])
	sorted_vals = [v for v in sorted_vals if v != None]
	max_val = sorted_vals[-1]
	hist_list = [0]*(max_val+1)
	for v in sorted_vals:
		hist_list[v]+=1
	hist_str = ['dist-%d:%d' % (i,hist_list[i]) for i in range(len(hist_list))]
	if verbose:
		print('%d total nodes: %d connected, %d disconnected' % (len(dist),len(sorted_vals),num_disconnected))
		print(', '.join(hist_str))

	return hist_list, num_disconnected

def to_graph_decomposition_considering_hypernodes(H,names, connect_hypernode_members=True):
	if not isinstance(H, directed_hypergraph.DirectedHypergraph):
		raise TypeError("Transformation only applicable to directed hypergraphs")

	G = directed_hypergraph.DirectedHypergraph()

	nodes = {}
	hypernode_edges = []
	hypernode_members = {}
	for node in H.node_iterator():
		if not H.get_node_attribute(node,'is_hypernode'):
			nodes[node] = H.get_node_attributes(node)
			nodes[node]['hypernode_parents'] = set()
			hypernode_members[node] = set([node])
		else:
			members = H.get_node_attribute(node,'hypernode_members')
			hypernode_members[node] = members
			for member in members:
				if member not in nodes:
					nodes[member] = {'name':names[member],'is_hypernode':False,'hypernode_members':set(),'hypernode_parents':set()}
				nodes[member]['hypernode_parents'].add(node)
			if connect_hypernode_members:
				for t,h in itertools.permutations(members,2):
					hypernode_edges.append( (set([t]),set([h]),{'hedge_type':'complex_membership'}) )

	G.add_nodes((n,nodes[n]) for n in nodes)
	
	print('%d nodes are hypernodes (SHOULD BE 0)' % 
		sum([G.get_node_attribute(n,'is_hypernode') for n in G.node_iterator()]))

   	hyperedge_edges = []
	for hyperedge_id in H.hyperedge_id_iterator():
		for tail_node in H.get_hyperedge_tail(hyperedge_id):
			for head_node in H.get_hyperedge_head(hyperedge_id):
				hyperedge_edges += [ 
					(set([t]),set([h]),{'hedge_type':'hyperedge_membership'})
						for t in hypernode_members[tail_node] 
						for h in hypernode_members[head_node] ]

	G.add_hyperedges(hypernode_edges)
	G.add_hyperedges(hyperedge_edges)
	print('Added %d nodes and %d edges (%d edges from hypernode membership and %d edges from hyperedge connections, many dups)' % (stats.number_of_nodes(G),stats.number_of_hyperedges(G),len(hypernode_edges),len(hyperedge_edges)))
	return G

def plot_hist(hgraph,graph,hgraph_decomp,graph_decomp,title,title_decomp,file_prefix):
	fig, (ax1,ax2) = plt.subplots(ncols=2, nrows=1, figsize=(10,5))
	bin_range = range(0,max(graph+graph_decomp)+10,10)
	ax1.hist(graph,bins=bin_range,color='#4C8DD6',edgecolor='k', linewidth=1.2,label='Graph Connectivity')
	ax1.hist(hgraph,bins=bin_range,color='#F5CBA7', rwidth=0.5, edgecolor='k', linewidth=1.2, label='Hypergraph B-Connectivity')
	ax1.set_ylabel('# of Source Nodes')
	ax1.set_xlabel('# of Reachable Nodes (%d total)' % (len(graph)))
	ax1.set_title(title)
	ax1.legend(loc='best',fontsize=10)

	ax2.hist(graph_decomp,bins=bin_range,color='#4C8DD6',edgecolor='k', linewidth=1.2,label='Graph Connectivity')
	ax2.hist(hgraph_decomp,bins=bin_range,color='#F5CBA7', rwidth=0.5, edgecolor='k', linewidth=1.2, label='Hypergraph B-Connectivity')
	ax2.set_ylabel('# of Source Nodes')
	ax2.set_xlabel('# of Reachable Nodes (%d total)' % (len(graph_decomp)))
	ax2.set_title(title_decomp)
	ax2.legend(loc='best',fontsize=10)

	plt.tight_layout()
	plt.savefig(file_prefix+'.png')
	print('saved to '+file_prefix+'.png')


if __name__ == '__main__':
	main()