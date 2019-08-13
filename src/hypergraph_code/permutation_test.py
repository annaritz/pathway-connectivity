# permutation test for pathways.
import sys
import os
from math import log
import glob
import networkx as nx
import random

MAX_TRIES = 1000000
def main(inprefix,outprefix,num_perms,num_swaps):
	pathways = read_files(inprefix)
	print('%d pathways' % (len(pathways)))

	G,set_names,node_membership = generate_graph(pathways)

	for permutation in range(num_perms):
		print('__PERMUTATION__%d__of__%d__' % (permutation+1,num_perms))
		fname = '%s%d_perms_%d_swaps.txt' % (outprefix,permutation,num_swaps)
		run_permutation(G,pathways,set_names,fname,num_swaps)

	return

def run_permutation(G,pathways,set_names,fname,num_swaps,verbose=True):
	# swap edges
	H = swap_edges(G,set_names,num_swaps,verbose)

	# check to make sure they're different. They are!
	#CHECK = nx.difference(G,H)
	#print(nx.info(CHECK))
	#sys.exit()

	# split data 
	perm_pathways = {p:set() for p in pathways.keys()}
	for p in set_names:
		pathway_set = split_set_label(p)
		nodes_to_add = H[p]
		for pathway in pathway_set:
			#print('ADDING %d to PATHWAY %s' % (len(nodes_to_add),pathway))
			perm_pathways[pathway].update(nodes_to_add)

	# write
	out = open(fname,'w')
	out.write('#Pathway\tNumMembers\tMembers\n')
	for p in perm_pathways:
		out.write('%s\t%d\t%s\n' % (p,len(perm_pathways[p]),';'.join(perm_pathways[p])))
	out.close()
	print('wrote to %s' % (fname))
	return

def generate_graph(pathways):
		## how many overlapping sets are there?
	node_membership = {}
	for p in pathways:
		for n in pathways[p]:
			if n not in node_membership:
				node_membership[n] = set()
			node_membership[n].add(p)
	print('%d molecules' % (len(node_membership)))

	sets = {}
	for n in node_membership:
		t = make_set_label(node_membership[n])
		if t not in sets:
			sets[t] = set()
		sets[t].add(n)
	print('%d sets' % (len(sets)))
	set_names = list(sets.keys())
	#for s in sets:
	#		print(len(sets[s]),s)

	## make bipartite graph with this information
	print()
	G = nx.Graph()
	G.add_nodes_from(sets.keys())
	G.add_nodes_from(node_membership.keys())
	print('%d nodes in bipartite graph' % (nx.number_of_nodes(G)))

	for t in sets:
		for n in sets[t]:
			G.add_edge(t,n)
	print('%d edges in bipartite graph' % (nx.number_of_edges(G)))
	return G,set_names,node_membership

def swap_edges(G,set_names,num_swaps,verbose=True):
	
	# start with a new copy of the graph
	H = nx.Graph()
	for n in G.nodes():
		H.add_node(n)
	for e in G.edges():
		H.add_edge(e[0],e[1])

	## swap edges. Can't use swap() methods in networkx because we need to make sure
	## bipartite graph remains bipartite.
	## MOdified from double_edge_swap
	## https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.swap.double_edge_swap.html#networkx.algorithms.swap.double_edge_swap
	## https://networkx.github.io/documentation/stable/_modules/networkx/algorithms/swap.html#double_edge_swap
	successful_swaps = 0
	n = 0
	while successful_swaps < num_swaps:
		n+=1
		if n % int(num_swaps/10) == 0 and verbose:
			print('%d tries & %d successful swaps' % (n,successful_swaps))
		if n > MAX_TRIES and verbose:
			print('MAX TRIES REACHED!! Quitting.')
			sys.exit()

		# pick two random edges without creating edge list
		# choose two source source nodes from pathway nodes.
		u = random.choice(set_names)
		x = random.choice(set_names)
		if u == x: # skip if they're the same node.
			continue
		# choose target uniformly from neighbors
		v = random.choice(list(H[u]))
		y = random.choice(list(H[x]))
		if v == y: # same target; skip
			continue
		# we now have (u,v) and (x,y).  
		if (x in H[u]) or (y in H[v]): # if edge already exists in G, skip.
			continue
		# swapping these edges won't result in parallel edges.
		H.add_edge(u,y)
		H.add_edge(x,v) 
		H.remove_edge(u,v)
		H.remove_edge(x,y)
		successful_swaps+=1
		#print('swapped (%s,%s) and (%s,%s)' % (u,v,x,y))		
	if verbose:
		print('FINAL: %d tries and %d successful swaps' % (n,successful_swaps))
	return H

def make_set_label(pathway_set):
	return '|'.join(sorted(list(pathway_set)))
def split_set_label(label):
	return label.split('|')


def read_files(prefix):
	files = glob.glob(prefix+'*')
	print('%d files' % (len(files)))
	pathways = {}
	for f in files:
		pathway_name = f.replace(prefix,'').replace('_b_relax.txt','')
		print('reading %s' % (pathway_name))
		pathways[pathway_name] = {}
		with open(f) as fin:
			for line in fin:
				if line[0:2] != '-1':
					continue
				row = line.strip().split('\t')
				pathways[pathway_name] = set(row[2].split(';'))
	return pathways

if __name__ == '__main__':
	
	if len(sys.argv) != 5:
		print('USAGE: python3 pathway-influence.py <PATHWAY_OUTPREFIX> <OUTPREFIX> <NUMPERMS> <NUMSWAPS>')

	random.seed(123456)
	main(sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]))
