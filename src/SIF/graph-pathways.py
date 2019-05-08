## computes the survey on SIF-formatted graphs.

import networkx as nx
import sys

def main(infile,outprefix,conversions,filter_file,pathway_file):
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
	nodes = G.nodes()
	with open(pathway_file) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split('\t')
			pathway_id = row[0]
			pathway_name = row[1]
			#pathway_num = int(row[2])
			pathway_members = set(row[3].split(';'))
			pathway_members = pathway_members.intersection(nodes)
			pathway_num = len(pathway_members)
			print(pathway_name,pathway_id,pathway_num)
			hist_dict = survey_graph_parameterized_fast(G,pathway_members)
						
			outfile = outprefix+pathway_name.replace(' ','-').replace('/','-')+'.txt'					
			out = open(outfile,'w')
			out.write('#k\t#Members\tMembers\n')
			out.write('-1\t%d\t%s\t%s\n' % (len(pathway_members),';'.join(pathway_members),'None')) # write original set
			for dist in range(max(hist_dict)+1):
				if dist in hist_dict:
					out.write('%d\t%d\t%s\n' % (dist,len(hist_dict[dist]),';'.join(hist_dict[dist])))
				else:
					out.write('%s\t0\t\t\n' % (dist))
			out.close()
			print('wrote to %s' % (outfile))

	print('Done')
	return

def survey_graph_parameterized(G,pathway_members):
	dist_dict = {}
	for node in pathway_members:
		this_dist_dict = bfs(G,node)
		if len(dist_dict) == 0:
			dist_dict = this_dist_dict
		else:
			for key in this_dist_dict:
				if key not in dist_dict or this_dist_dict[key] < dist_dict[key]:
					dist_dict[key] = this_dist_dict[key]
	hist_dict = dist2hist(dist_dict)
	return hist_dict

def survey_graph_parameterized_fast(G,pathway_members):
	ss = 'ss'
	for node in pathway_members:
		G.add_edge(ss,node)
	this_dist_dict = bfs(G,node)
	## decrement all distances by one (since we started at ss)
	dist_dict = {}
	for d in this_dist_dict:
		if d == ss:
			continue
		dist_dict[d]=this_dist_dict[d]-1
		
	for node in pathway_members:
		G.remove_edge(ss,node)
	hist_dict = dist2hist(dist_dict)
	return hist_dict



def bfs(G,s):
	succ = nx.bfs_successors(G,s)
	dist_sets = {}
	curr_dist = 0
	to_traverse = set([s])
	while len(to_traverse) > 0:
		#print('CURR DIST',curr_dist)
		#print('TO TRAVERSE',to_traverse)
		for t in to_traverse:
			dist_sets[t] = curr_dist
		traverse_next = set()
		for t in to_traverse:
			if t in succ:
				traverse_next.update(succ[t])
			# if t is not in successors, it has no outgoing edges in the BFS tree. This is fine.
		curr_dist += 1
		to_traverse = traverse_next
	return dist_sets

def dist2hist(dist_dict):
	## get histogram of distances.
	h = {} # distance: # of nodes
	for n,val in dist_dict.items():
		if val == None: # skip 'None' type (these are infinity)
			continue
		if val not in h:
			h[val] = set()
		h[val].add(n)
	return h

if __name__ == '__main__':
	if len(sys.argv) != 5 and len(sys.argv) != 6:
		print('USAGE: python3 GraphSurvey.py <infile> <outprefix> <conversion-type-file> <pathway file> <filter-file-OPTIONAL>')
	infile = sys.argv[1]
	outprefix = sys.argv[2]
	convfile = sys.argv[3]
	pathway_file = sys.argv[4]
	if len(sys.argv)==6:
		filter_file = sys.argv[5]
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

	main(infile,outprefix,conversions,filter_file,pathway_file)