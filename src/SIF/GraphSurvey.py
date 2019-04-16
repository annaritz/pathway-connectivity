## computes the survey on SIF-formatted graphs.

import networkx as nx
import sys

def main(infile,outfile,conversions,filter_file):
	## use a directed graph.
	G = nx.DiGraph()
	with open(infile) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			# SIF format is <n1 TYPE n2>
			row = line.strip().split('\t')
			n1 = row[0]
			ntype = row[1]
			n2 = row[2]

			if ntype not in conversions:
				print('ERROR:',ntype)
				sys.exit()

			if conversions[ntype] == 'IGNORE':
				continue
			elif conversions[ntype] == 'DIR':
				G.add_edge(n1,n2)
			elif conversions[ntype] == 'UNDIR':
				G.add_edge(n1,n2)
				G.add_edge(n2,n1)
			else:
				print('ERROR:',ntype,conversions[ntype])
				sys.exit()
	print('Graph has %d nodes and %d edges' % (nx.number_of_nodes(G),nx.number_of_edges(G)))

	if filter_file:
		#nodes_to_keep = set([n for n in G.nodes() if 'CHEBI' not in n])
		nodes_to_keep = set()
		with open('../../hypergraph/reactome_hypergraphs/proteins.txt') as fin:
			for line in fin:
				nodes_to_keep.add(line.strip())
		#print(nodes_to_keep)
		G = G.subgraph(nodes_to_keep)
		print('Graph has %d nodes and %d edges' % (nx.number_of_nodes(G),nx.number_of_edges(G)))
	#nx.write_edgelist(G,outfile+'.graph',delimiter='\t',data=False)
	#survey_graph(G,outfile)
	survey_graph_parameterized(G,outfile+'.parameterized')

	return

def survey_graph(G,outfile):
	out = open(outfile,'w')
	out.write('#Name\tNumConnected\n')
	i = 0
	for n in G.nodes():
		i+=1
		if i % 100 == 0:
			print(' %d of %d' % (i,nx.number_of_nodes(G)))
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

	return

if __name__ == '__main__':
	if len(sys.argv) != 4 and len(sys.argv) != 5:
		print('USAGE: python3 GraphSurvey.py <infile> <outfile> <conversion-type-file> <filter-fileOPTIONAL>')
	infile = sys.argv[1]
	outfile = sys.argv[2]
	convfile = sys.argv[3]
	if len(sys.argv)==5:
		filter_file = sys.argv[4]
	else:
		filter_file=None
	conversions = {}
	with open(convfile) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split('\t')
			conversions[row[0]] = row[1]
	print('%d conversions read' % (len(conversions)))

	main(infile,outfile,conversions,filter_file)