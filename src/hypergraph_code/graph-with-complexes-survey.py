## computes the survey on the graph-with-complexes.

import networkx as nx
import sys
import itertools
import hgraph_utils
from halp import directed_hypergraph
from halp.directed_hypergraph import DirectedHypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform

def main(inprefix,outfile,filter_file):
	H, identifier2id, id2identifier = hgraph_utils.make_hypergraph(inprefix)
	G = transform.to_networkx_digraph(H)

	if filter_file:
		#nodes_to_keep = set([n for n in G.nodes() if 'CHEBI' not in n])
		nodes_to_keep = set()
		with open('../../hypergraph/reactome_hypergraphs/proteins.txt') as fin:
			for line in fin:
				nodes_to_keep.add(line.strip())
		#print(nodes_to_keep)
		G = G.subgraph(nodes_to_keep)
		print('Graph has %d nodes and %d edges' % (nx.number_of_nodes(G),nx.number_of_edges(G)))

	nx.write_edgelist(G,outfile+'.graph',delimiter='\t',data=False)
	print('wrote networkx graph to %s.graph' % (outfile))
	print('Graph has %d nodes and %d edges' % (nx.number_of_nodes(G),nx.number_of_edges(G)))

	# print('MAKING BIPARTITE GRAPH')
	# B = nx.DiGraph()
	# B.add_node_set(H.get_nodes())
	# B.add_node_set(H.get_hyperedge_ids())
	# nx.write_edgelist(B,outfile+'.bipartite_graph',delimiter='\t',data=False)
	# print('wrote networkx graph to %s.graph' % (outfile))
	# print('Graph has %d nodes and %d edges' % (nx.number_of_nodes(B),nx.number_of_edges(B)))

	#print('EXITING BEFORE SURVEY')
	#sys.exit()

	survey_graph(G,outfile)
	survey_graph_parameterized(G,outfile+'.parameterized')
	return

def survey_graph(G,outfile):
	out = open(outfile,'w')
	out.write('#Name\tNumConnected\n')
	i = 0
	for n in G.nodes():
		i+=1
		if i % 100 == 0:
			print(' %d of %d' % (i,len(G.nodes())))
		#print(n)
		edges = nx.bfs_edges(G,n)
		nset = set()
		for u,v in edges:
			nset.add(u)
			nset.add(v)
		out.write('%s\t%d\n' % (n,len(nset)))

	out.close()
	print('wrote to %s' % (outfile))

	return

def survey_graph_parameterized(G,outfile):
	out = open(outfile,'w')
	out.write('#Name\td=1\td=2\td=3\t...\n')
	i = 0
	for s in G.nodes():
		i+=1
		if i % 100 == 0:
			print(' %d of %d' % (i,nx.number_of_nodes(G)))
		#print(n)
		succ = nx.bfs_successors(G,s)
		dist_sets = {}
		curr_dist = 0
		to_traverse = set([s])
		while len(to_traverse) > 0:
			#print('CURR DIST',curr_dist)
			#print('TO TRAVERSE',to_traverse)
			dist_sets[curr_dist] = set([t for t in to_traverse])
			if curr_dist > 0:
				dist_sets[curr_dist].update([t for t in dist_sets[curr_dist-1]])
			traverse_next = set()
			for t in to_traverse:
				if t in succ:
					traverse_next.update(succ[t])
				# if t is not in successors, it has no outgoing edges in the BFS tree. This is fine.
			curr_dist += 1
			to_traverse = traverse_next
		distances = [str(len(dist_sets[d])) for d in range(1,curr_dist)] # skip d=0
		out.write('%s\t%s\n' % (s,'\t'.join(distances)))

	out.close()
	print('wrote to %s' % (outfile))

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print('USAGE: python3 GraphSurvey.py <inprefix> <outfile>')
	inprefix = sys.argv[1]
	outfile = sys.argv[2]
	main(inprefix,outfile)