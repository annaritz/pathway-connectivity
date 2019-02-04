import sys
import itertools
import hgraph_utils
from halp import directed_hypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform
import time

def main(inprefix,outfile,hedge_connectivity_file):
	H, identifier2id, id2identifier = hgraph_utils.make_hypergraph(inprefix)
	b_visit_dict = hgraph_utils.make_b_visit_dict(hedge_connectivity_file,identifier2id)

	#b_relaxation_survey_nodes(H,b_visit_dict,report_single='http://pathwaycommons.org/pc2/Rna_5298e36420325d9483e9edee07ee17af')
	all_dist_dicts,times,max_val = b_relaxation_survey_nodes(H,b_visit_dict)

	out = open(outfile,'w')
	out.write('#Node\tTime\t%s\n' % ('\t'.join([str(i) for i in range(max_val+1)])))
	for node in all_dist_dicts:
		#print(node)
		#print(all_dist_dicts[node])
		distlist = [all_dist_dicts[identifier].get(i,0) for i in range(max_val+1)]
		#print(distlist)
		cumulist = [sum(distlist[:i]) for i in range(max_val+1)]
		#print(cumulist)
		out.write('%s\t%f\t%s\n' % (node,times[node],'\t'.join([str(i) for i in cumulist])))
	out.close()
	print('wrote to %s' % (outfile))

	## old implementation -- go by hyperedges not nodes.
	## all_dist_dicts,times,max_val = b_relaxation_survey(H,b_visit_dict)
	return

def b_relaxation_survey(H,b_visit_dict,report_single=None):
	i = 0
	max_val = 0
	all_dist_dicts = {}
	times = {}
	for hedge_id in H.hyperedge_id_iterator():
		identifier = H.get_hyperedge_attribute(hedge_id,'identifier')
		if report_single and hedge_id != report_single:
			continue
		if report_single:
			print(hedge_id,identifier)
			print('TAILS:',H.get_hyperedge_tail(hedge_id))
			print('HEADS:',H.get_hyperedge_head(hedge_id))
			print('B VISIT DICT:',b_visit_dict[hedge_id])
		i+=1
		if i % 50 == 0:
			print('hyperedge %d of %d' % (i,stats.number_of_hyperedges(H)))
			break
		start_time = time.time()
		dist_dict = hpaths.b_relaxation(H,H.get_hyperedge_tail(hedge_id),b_visit_dict=b_visit_dict)
		end_time = time.time()
		print(hedge_id,end_time-start_time)
		times[identifier] = end_time-start_time
		hist_dict = dist2hist(dist_dict)
		if hist_dict != {}:
			max_val = max(max(hist_dict.keys()),max_val)
		all_dist_dicts[identifier] = hist_dict
		if report_single:
			for i in range(1,3):
				print('CONNECTED K=%d' % (i))
				for d in dist_dict:
					if dist_dict[d] == i:
						print(d,dist_dict[d])
			sys.exit()

	return all_dist_dicts,times, max_val

def b_relaxation_survey_nodes(H,b_visit_dict,report_single=None):
	i = 0
	max_val = 0
	all_dist_dicts = {}
	times = {}
	for node in H.node_iterator():
		if report_single and node != report_single:
			continue
		if report_single:
			print(node)
		i+=1
		if i % 50 == 0:
			print('node %d of %d' % (i,stats.number_of_nodes(H)))
			break
		start_time = time.time()
		dist_dict = hpaths.b_relaxation(H,set([node]),b_visit_dict=b_visit_dict)
		end_time = time.time()
		print(node,end_time-start_time)
		times[node] = end_time-start_time
		hist_dict = dist2hist(dist_dict)
		if hist_dict != {}:
			max_val = max(max(hist_dict.keys()),max_val)
		all_dist_dicts[node] = hist_dict
		if report_single:
			for i in range(1,3):
				print('CONNECTED K=%d' % (i))
				for d in dist_dict:
					if dist_dict[d] == i:
						print(d,dist_dict[d])
			sys.exit()

	return all_dist_dicts,times, max_val

def dist2hist(dist_dict):
	## get histogram of distances.
	h = {} # distance: # of nodes
	for n,val in dist_dict.items():
		if not val: # skip 'none' type
			continue
		if val not in h:
			h[val] = 0
		h[val]+=1
	return h

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print('USAGE: python b-relaxation-survey.py <INPREFIX> <outfile> <HEDGE_CONNECTIVITY_FILE')
		sys.exit()
	main(sys.argv[1],sys.argv[2],sys.argv[3])