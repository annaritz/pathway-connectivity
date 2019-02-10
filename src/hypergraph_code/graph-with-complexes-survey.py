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

def main(inprefix,outfile):
	H, identifier2id, id2identifier = hgraph_utils.make_hypergraph(inprefix)
	G = transform.to_networkx_digraph(H)
	nx.write_edgelist(G,outfile+'.graph',delimiter='\t',data=False)
	print('wrote networkx graph to %s.graph' % (outfile))
	print('Graph has %d nodes and %d edges' % (nx.number_of_nodes(G),nx.number_of_edges(G)))

	print('MAKING BIPARTITE GRAPH')
	B = nx.DiGraph()
	B.add_node_set(H.get_nodes())
	B.add_node_set(H.get_hyperedge_ids())
	nx.write_edgelist(B,outfile+'.bipartite_graph',delimiter='\t',data=False)
	print('wrote networkx graph to %s.graph' % (outfile))
	print('Graph has %d nodes and %d edges' % (nx.number_of_nodes(B),nx.number_of_edges(B)))

	print('EXITING BEFORE SURVEY')
	sys.exit()
	survey_graph(G,outfile)
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

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print('USAGE: python3 GraphSurvey.py <inprefix> <outfile>')
	inprefix = sys.argv[1]
	outfile = sys.argv[2]
	main(inprefix,outfile)