import sys
import itertools
import hgraph_utils
from halp import directed_hypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform
import time
## https://docs.python.org/3.4/library/multiprocessing.html?highlight=process
from multiprocessing import Pool


def main(inprefix,outprefix,hedge_connectivity_file,pathway_prefix,num_perms, num_swaps):
	H, identifier2id, id2identifier = hgraph_utils.make_hypergraph(inprefix)
	nodes = H.get_node_set()
	b_visit_dict = hgraph_utils.make_b_visit_dict(hedge_connectivity_file,identifier2id)

	params = [(pathway_prefix,outprefix,H,nodes,b_visit_dict,p,num_swaps) for p in range(num_perms)]
	print('Running %d Params' % (len(params)))
	print('Example:',len(params[0]),type(params[0]))
	with Pool(4) as p:
		p.map(run_instance,params)
	print('Done')
	return

def run_instance(param_tuple):
	pathway_prefix,outprefix,H,nodes,b_visit_dict,perm,swap = param_tuple
	pathway_file = '%s%d_perms_%d_swaps.txt' % (pathway_prefix,perm,swap)

	with open(pathway_file) as fin: # each line lists a pathway.
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
			hist_dict,traversed = b_relaxation_survey_nodes(H,b_visit_dict,pathway_members,v=False)
			
			outfile = '%s%d_perms_%d_swaps_%s_b_relax.txt' % (outprefix,perm,swap,pathway_name.replace(' ','-').replace('/','-'))
			out = open(outfile,'w')
			out.write('#k\t#Members\tMembers\tIncidentNodes\n')
			out.write('-1\t%d\t%s\t%s\n' % (len(pathway_members),';'.join(pathway_members),'None')) # write original set
			for dist in range(max(hist_dict)+1):
				if dist in hist_dict:
					out.write('%d\t%d\t%s\t%s\n' % (dist,len(hist_dict[dist]),';'.join(hist_dict[dist]),';'.join(traversed[dist])))
				else:
					out.write('%s\t0\t\t\t\n' % (dist))
			out.close()
			print('wrote to %s' % (outfile))
	return

def b_relaxation_survey_nodes(H,b_visit_dict,pathway_members,v=False):
	i = 0
	max_val = 0
	all_dist_dicts = {}
	times = {}
	dist_dict,traversed = hpaths.b_relaxation(H,pathway_members,b_visit_dict=b_visit_dict,verbose=v)
	hist_dict = dist2hist(dist_dict)
	return hist_dict,traversed

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
	if len(sys.argv) != 7:
		print('USAGE: python b-relaxation-permutation.py <INPREFIX> <OUTPREFIX> <HEDGE_CONNECTIVITY_FILE> <PATHWAY_PREFIX>  <NUMPERMS> <NUMSWAPS>')
		sys.exit()
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],int(sys.argv[5]),int(sys.argv[6]))

