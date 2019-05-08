import sys
import itertools
import hgraph_utils
from halp import directed_hypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform
import time
import networkx as nx

def main(inprefix,outprefix,pathway_file):
	H, identifier2id, id2identifier = hgraph_utils.make_hypergraph(inprefix)
	nodes = H.get_node_set()
	G = transform.to_networkx_digraph(H)
	
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
			if 'MST1' in pathway_name:
				v = True
			else:
				v = False
			hist_dict = parameterized_survey_nodes(G,pathway_members,v=v)
			
			outfile = outprefix+pathway_name.replace(' ','-').replace('/','-')+'_parameterized.txt'					
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

def parameterized_survey_nodes(G,pathway_members,v=False):
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
	if len(sys.argv) != 4:
		print('USAGE: python b-relaxation-survey.py <INPREFIX> <OUTPREFIX> <PATHWAYS>')
		sys.exit()
	main(sys.argv[1],sys.argv[2],sys.argv[3])

