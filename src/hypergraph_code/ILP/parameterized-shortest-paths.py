import sys
import itertools
import hgraph_utils as hgraph_utils
from halp import directed_hypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform
import time
import shortest_hyperpath

def main(inprefix,outfile):
	H, identifier2id, id2identifier = hgraph_utils.make_hypergraph(inprefix)

	out = open(outfile,'w')
	out.write('#Node\td=1\td=2\td=3\t...\n')
	
	bconn_param(H,out)
		
	out.close()
	print('wrote to %s' % (outfile))

	## old implementation -- go by hyperedges not nodes.
	## all_dist_dicts,times,max_val = b_relaxation_survey(H,b_visit_dict)
	return

def bconn_param(H,out):
	i = 0
	for node in H.node_iterator():
		#print('Node',node)
		i+=1
		if i % 50 == 0:
			print('node %d of %d' % (i,stats.number_of_nodes(H)))

		bconn, ignore, ignore, ignore = hpaths.b_visit(H,node)
		#print('B-conn',bconn)
		if len(bconn) == 1: ## only one node; skip
			out.write('%s\t1\n' % (node))
			continue

		G = H.get_induced_subhypergraph(bconn)
		#print('B-conn hypergraph has %d hyperedges and %d nodes' % (stats.number_of_hyperedges(G),stats.number_of_nodes(G)))

		dist_dict = {}
		for t in bconn:		
			# get distance from node to t using ILP
			dist_dict[t] = shortest_hyperpath.runILP(G,node,t)
			#sys.exit()
		#print(dist_dict)

		
		## convert distances to cumulative list.
		hist_dict = dist2hist(dist_dict)
		max_val = int(max(hist_dict.keys()))
		distlist = [hist_dict.get(i,0) for i in range(max_val+1)]
		cumulist = [sum(distlist[:i+1]) for i in range(max_val+1)]
		#all_dist_dicts[node] = cumulist
		out.write('%s\t%s\n' % (node,'\t'.join([str(i) for i in cumulist])))

		#print(node)
		#print('B-conn hypergraph has %d hyperedges and %d nodes' % (stats.number_of_hyperedges(G),stats.number_of_nodes(G)))
		#print(hist_dict)
		#print(distlist)
		#print(cumulist)
		#sys.exit()

	return 

def dist2hist(dist_dict):
	## get histogram of distances.
	h = {} # distance: # of nodes
	for n,val in dist_dict.items():
		if val == None: # skip 'None' type (these are infinity)
			continue
		if val not in h:
			h[val] = 0
		h[val]+=1
	return h

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print('USAGE: python3 parametereized-shortest-path <INPREFIX> <outfile>')
		sys.exit()
	main(sys.argv[1],sys.argv[2])