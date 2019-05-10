import sys
import hgraph_utils
from halp import directed_hypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform
import glob
# library
import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn3

LARGE_VAL = 10000000

def main(inprefix,hedge_connectivity_file,pathway_prefix):

	# get pathway Identifiers to Uniprot ID
	pc2uniprot,uniprot2pc = hgraph_utils.get_id_map('../../hypergraph/reactome_hypergraphs/')

	interactions = read_functional_interactions()
	print('\n%d INTERACTIONS\n' % (len(interactions)))


	test_set = set([tuple([n1,n2]) for n1,n2,val in interactions])
	print(list(test_set)[:10])
	print('%d in test set' % (len(test_set)))

	H, identifier2id, id2identifier = hgraph_utils.make_hypergraph(inprefix)
	b_visit_dict = hgraph_utils.make_b_visit_dict(hedge_connectivity_file,identifier2id)

	nodes = set() ## get proteins and complex members.
	for n in H.get_node_set():
		attrs = H.get_node_attributes(n)
		if attrs['is_hypernode']:
			nodes.update(attrs['hypernode_members'])
		else:
			nodes.add(n)

	print('%d nodes including hypernode members' % (len(nodes)))	

	interactions_in_reactome = []
	mismapped = 0
	notinreactome = 0
	for n1,n2,val in interactions:
		if n1 in uniprot2pc and n2 in uniprot2pc:
			un1 = uniprot2pc[n1]
			un2 = uniprot2pc[n2]
		else:
			mismapped+=1
			continue
		
		if un1 in nodes and un2 in nodes:
			interactions_in_reactome.append([un1,un2,val])
		else:
			#print('*',n1,un1,un2 in nodes)
			#print('*',n2,un2,un2 in nodes)
			#sys.exit()
			notinreactome+=1
	print('\n%d INTERACTIONS HAVE BOTH NODES IN REACTOME\n%d interactions not in Reactome mapping; %d not in this hypergraph\n' % (len(interactions_in_reactome),mismapped,notinreactome))

	## TODO I suspect the interactions that have proteins IN Reactome mapping but NOT IN this hypergraph are due to 
	## hypernode/complexes; which are not considered.  Expand to look in hypernodes (complexes and entity sets).

	## NEW STEP - get interactions in pathway.
	pathways = read_pathways_from_brelax(pathway_prefix)
	print('%d pathways' % (len(pathways)))
	all_pathway_nodes = set()
	pathway_nodes = {}
	for p in pathways:
		pathway_nodes[p] = set()
		for n in pathways[p]:
			attrs = H.get_node_attributes(n)
			if attrs['is_hypernode']:
				pathway_nodes[p].update(attrs['hypernode_members'])
			else:
				pathway_nodes[p].add(n)
		all_pathway_nodes.update(pathway_nodes[p])
	print('%d pathway nodes (including hypernode members)' % (len(all_pathway_nodes)))
	print(list(all_pathway_nodes)[:10])

	interactions_in_pathways = []
	for n1,n2,val in interactions_in_reactome:
		if n1 in all_pathway_nodes and n2 in all_pathway_nodes:
			interactions_in_pathways.append([n1,n2,val])
	print('\n%d INTERACTIONS HAVE BOTH NODES IN REACTOME PATHWAYS\n' % (len(interactions_in_pathways)))

	
	pos_1_set = set()
	for n1,n2,val in interactions_in_pathways:
		for p in pathway_nodes:
			if n1 in pathway_nodes[p] and n2 in pathway_nodes[p]:
				pos_1_set.add(tuple([n1,n2]))
				
				break
	print('\n%d INTERACTIONS HAVE BOTH NODES IN SAME REACTOME PATHWAY\n' % (len(pos_1_set)))
	
	#precision(interactions_in_pathways,[pos_1_set],['test'])

	final_nodes = {}
	for n1,n2,val in interactions_in_pathways:
		if n1 not in final_nodes:
			final_nodes[n1] = set()
		if n2 not in final_nodes:
			final_nodes[n2] = set()
	final_node_set = set(final_nodes.keys())
	print('\n%d NODES IN THE PATHWAY INTERACTIONS' % (len(final_nodes)))

	### get mapper of proteins to hypernodes (e.g. complexes.)
	for n in H.get_node_set():
		if n in final_node_set:
			final_nodes[n].add(n)
		else:
			attrs = H.get_node_attributes(n)
			if attrs['is_hypernode']:
				for member in attrs['hypernode_members']:
					if member in final_node_set:
						final_nodes[member].add(n)
	#for n in final_nodes:
	#	print(n,len(final_nodes[n]))

	interactions_in_pathways.sort(key=lambda x: x[2],reverse=True)
	i = 0
	pos_2_set = set()
	pos_3_set = set()
	for n1,n2,val in interactions_in_pathways:
		i+=1
		
		n1_nodes = final_nodes[n1]
		n2_nodes = final_nodes[n2]
		dist_dict_n1,ignore = hpaths.b_relaxation(H,n1_nodes,b_visit_dict=b_visit_dict)
		dist_dict_n2,ignore = hpaths.b_relaxation(H,n2_nodes,b_visit_dict=b_visit_dict)
		score = get_min_dist(n1_nodes,dist_dict_n1,n2_nodes,dist_dict_n2)
		#print(score,val,n1,n2)
		insamepathway=False
		for p in pathway_nodes:
			if n1 in pathway_nodes[p] and n2 in pathway_nodes[p]:
				insamepathway=True
		if score < LARGE_VAL:
			#if not insamepathway:
			#	print('HOLY COW! NOT IN SAME PATHWAY!')
			pos_2_set.add(tuple([n1,n2]))
			if score == 0:
				pos_3_set.add(tuple([n1,n2]))
		#else:
		#	if insamepathway:
		#		print('HOLY COW! IN SAME PATHWAY!')
	print('\n%d INTERACTIONS ARE Bipartite CONNECTED IN REACTOME\n' % (len(pos_2_set)))
	print('\n%d INTERACTIONS ARE B-CONNECTED IN REACTOME\n' % (len(pos_3_set)))
	
	venn(pos_1_set,pos_2_set,pos_3_set)
	precision(interactions_in_pathways,[pos_1_set,pos_2_set,pos_3_set],['Pair in Same Pathway','Pair Connected (Bipartite)','Pair B-Connected'])

def venn(pos_1_set,pos_2_set,pos_3_set):
	plt.figure(figsize=(8,6))
	venn3([pos_1_set,pos_2_set,pos_3_set],('Pair in Same Pathway','Pair Connected (Bipartite)','Pair B-Connected'))
	filename = 'venn.png'
	plt.savefig(filename)
	print('wrote to '+filename)
	filename = filename.replace('png','pdf')
	plt.savefig(filename)
	print('wrote to '+filename)

def precision(interactions_in_pathways,pos_sets,pos_names):
	plt.figure(figsize=(8,6))
	ax = plt.subplot(1,1,1)
	interactions_in_pathways.sort(key=lambda x: x[2],reverse=True)
	x = range(len(interactions_in_pathways))
	ys = []
	for p in pos_sets:
		ys.append([])
	for n1,n2,val in interactions_in_pathways:
		node = tuple([n1,n2])
		for i in range(len(pos_sets)):
			if node in pos_sets[i]:
				if len(ys[i]) == 0:
					ys[i].append(1)
				else:
					ys[i].append(ys[i][-1]+1)
			else:
				if len(ys[0])==0:
					ys[i].append(0)
				else:
					ys[i].append(ys[i][-1])
	
	# normalize ys
	for j in range(len(ys)):
		for i in range(len(x)):
			ys[j][i] = ys[j][i]/float(i+1)

	for i in range(len(pos_sets)):
		print(x)
		print(ys[i])
		print(pos_names[i])
		ax.plot(x,ys[i],lw=2,label=pos_names[i])

	ax.set_xlabel('Reactome Interactions (Both Members are in some Pathway)')
	ax.set_ylabel('Precision')
	ax.set_title('Precision of the SAME ranked interactions on DIFFERENT positive sets')
	plt.legend()
	plt.tight_layout()
	filename = 'precision.png'
	plt.savefig(filename)
	print('wrote to '+filename)
	filename = filename.replace('png','pdf')
	plt.savefig(filename)
	print('wrote to '+filename)

def get_min_dist(n1_nodes,dist_dict_n1,n2_nodes,dist_dict_n2):
	val = LARGE_VAL
	for n in n2_nodes:
		if dist_dict_n1[n] != None:
			val = min(dist_dict_n1[n],val)
	for n in n1_nodes:
		if dist_dict_n2[n] != None:
			val = min(dist_dict_n2[n],val)
	return val

def read_functional_interactions():
	svd_file = '../../data/STRING/96066.protein.links.cooccurence.v11.0.txt.mapped'
	interactions = []
	with open(svd_file) as fin:
		for line in fin:
			if line[0] == '#': 
				continue
			row = line.strip().split()
			interactions.append([row[0],row[1],int(row[2])])
	return interactions

def read_pathways_from_brelax(prefix):
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
	if len(sys.argv) != 2:
		print('Usage: python3 run-benchmark.py <HYPERGRAPH_PREFIX> <hedge_connectivity_file> <pathway_brelax_prefix>')
	main(sys.argv[1],sys.argv[2],sys.argv[3])
