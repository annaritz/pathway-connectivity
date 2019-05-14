import networkx as nx
import sys

def read_graph(infile,convfile):
	## read conversion types
	conversions = {}
	with open(convfile) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split('\t')
			conversions[row[0]] = row[1]

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
	return G

def bfs(G,s):
	# bfs generates a dictionary of {node:dist} from s.
	succ = nx.bfs_successors(G,s)
	node_distances = {}
	curr_dist = 0
	to_traverse = set([s])
	while len(to_traverse) > 0:
		#print('CURR DIST',curr_dist)
		#print('TO TRAVERSE',to_traverse)
		for t in to_traverse:
			node_distances[t] = curr_dist
		traverse_next = set()
		for t in to_traverse:
			if t in succ:
				traverse_next.update(succ[t])
			# if t is not in successors, it has no outgoing edges in the BFS tree. This is fine.
		curr_dist += 1
		to_traverse = traverse_next
	return node_distances

def bfs_histogram(G,s):
	# BFS generates a dictionary of {distance: # within distance} from s.
	succ = nx.bfs_successors(G,s)
	dist_sets = {}
	curr_dist = 0
	to_traverse = set([s])
	while len(to_traverse) > 0:
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
	return dist_sets

def dist2hist(dist_dict,counts=False):
	## get histogram of distances.
	h = {} # distance: # of nodes
	for n,val in dist_dict.items():
		if val == None: # skip 'None' type (these are infinity)
			continue
		if val not in h:
			if counts:
				h[val] = 0
			else:
				h[val] = set()
		if counts:
			h[val]+=1
		else:
			h[val].add(n)
	return h