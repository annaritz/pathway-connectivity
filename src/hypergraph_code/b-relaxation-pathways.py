import sys
import itertools
import hgraph_utils
from halp import directed_hypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform
import time

def main(inprefix,outprefix,hedge_connectivity_file,pathway_file):
	H, identifier2id, id2identifier = hgraph_utils.make_hypergraph(inprefix)
	nodes = H.get_node_set()
	b_visit_dict = hgraph_utils.make_b_visit_dict(hedge_connectivity_file,identifier2id)
	
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

			hist_dict = b_relaxation_survey_nodes(H,b_visit_dict,pathway_members)

			outfile = outprefix+pathway_name.replace(' ','-').replace('/','-')+'_b_relax.txt'					
			out = open(outfile,'w')
			out.write('#k\t#Members\tMembers\n')
			out.write('-1\t%d\t%s\n' % (len(pathway_members),';'.join(pathway_members))) # write original set
			for dist in range(max(hist_dict)+1):
				if dist in hist_dict:
					out.write('%d\t%d\t%s\n' % (dist,len(hist_dict[dist]),';'.join(hist_dict[dist])))
				else:
					out.write('%s\t0\t\t\n' % (dist))
			out.close()
			print('wrote to %s' % (outfile))

	print('Done')
	return

def b_relaxation_survey_nodes(H,b_visit_dict,pathway_members):
	i = 0
	max_val = 0
	all_dist_dicts = {}
	times = {}
	dist_dict = hpaths.b_relaxation(H,pathway_members,b_visit_dict=b_visit_dict)
	hist_dict = dist2hist(dist_dict)
	return hist_dict

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
	if len(sys.argv) != 5:
		print('USAGE: python b-relaxation-survey.py <INPREFIX> <OUTPREFIX> <HEDGE_CONNECTIVITY_FILE> <PATHWAYS>')
		sys.exit()
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

