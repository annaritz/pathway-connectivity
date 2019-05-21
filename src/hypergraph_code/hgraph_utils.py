from halp import directed_hypergraph
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform
import glob
import networkx as nx
import sys

def make_hypergraph(file_prefix,delim=';',sep='\t',keep_singleton_nodes=False):
	hypernodes = {}
	with open(file_prefix+'-hypernodes.txt') as fin:
		for line in fin:
			if line[0] == '#': 
				continue
			row = line.strip().split(sep)
			## TODO fix -- this happens when a complex contains other complexes
			if len(row) == 1:
				hypernodes[row[0]] = ['OtherComplexes-FIX']
			else:
				hypernodes[row[0]] = row[1].split(delim)
	print('%d hypernodes from hypernodes file' % (len(hypernodes)))
	identifier2id = {}
	id2identifier = {}
	H = directed_hypergraph.DirectedHypergraph()
	if keep_singleton_nodes:
		for n in hypernodes:
			H.add_node(n)

	skipped1 = 0
	skipped2 = 0
	with open(file_prefix+'-hyperedges.txt') as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split(sep)
			tail = set()
			head = set()

			## Tail includes tail and regulators.
			## Head includes head.
			if row[0] != 'None':
				tail.update(row[0].split(delim))
			if row[1] != 'None':
				head.update(row[1].split(delim))
			if row[2] != 'None':
				tail.update(row[2].split(delim))
			if row[3] != 'None':
				tail.update(row[3].split(delim))
			hedge_id = row[4]

			## THIS IS A HACK FOR NOW ( should be incorporated in the make-hypergraph.py code)
			## IGnore any reactions that have a Reactome Identifier (e.g. has "HSA") instead of 
			## a PAthway Commons identifier.
			if any(['HSA' in s for s in tail]+['HSA' in s for s in head]):
				skipped1+=1
			elif len(tail)==0 or len(head)==0:
				skipped2+=1
			else:
				hid = H.add_hyperedge(tail,head,identifier=hedge_id)
				identifier2id[hedge_id] = hid 
				id2identifier[hid] = hedge_id
				
	print('%d reactions skipped because of Reactome identifier' % (skipped1))
	print('%d reactions skipped because of an empty tail or head' % (skipped2))
	## annotate nodes
	num_hypernodes = 0
	for node in H.get_node_set():
		if node in hypernodes and hypernodes[node] != [node]:
			H.add_node(node,hypernode_members=hypernodes[node],is_hypernode=True)
			num_hypernodes+=1
		else:
			H.add_node(node,is_hypernode=False,hypernode_members=[])

		H.add_node(node)

	print('Hypergraph has %d hyperedges and %d nodes (%d of these are hypernodes)' % 
		(stats.number_of_hyperedges(H),stats.number_of_nodes(H),num_hypernodes))

	return H, identifier2id, id2identifier

def check_singleton_nodes(H,file_prefix,entityset_prefix,outfile=None,sep='\t'):
	# get entity set 
	
	candidate_nodes = set()
	candidate_members = set()
	files = glob.glob('%s*-entitysets.txt' % (entityset_prefix))
	#print('checking %s*-entitysets.txt' % (entityset_prefix))
	for f in files:
		with open(f) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				row = line.strip().split()
				candidate_nodes.add(row[0])
				candidate_members.update(set(row[2].split(';')))
	files = glob.glob('%s*-complexes.txt' % (entityset_prefix))
	#print('checking %s*-entitysets.txt' % (entityset_prefix))
	for f in files:
		with open(f) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				row = line.strip().split()
				candidate_nodes.add(row[0])
				candidate_members.update(set(row[4].split(';')))

	singleton_nodes = set()
	nodes = H.get_node_set()
	i=0
	num_members=0
	num_rsa_id = 0
	others=0
	with open(file_prefix+'-hypernodes.txt') as fin:
		for line in fin:
			i+=1
			if line[0] == '#': 
				continue
			row = line.strip().split(sep)
			n = row[0]
			if n not in nodes or (n in nodes and len(H.get_backward_star(n))==0 and len(H.get_backward_star(n))==0):
				## singleton candidate. in ES?  
				if n in candidate_nodes:
					singleton_nodes.add(n)
				elif n in candidate_members:
					num_members +=1
				elif 'R-HSA' in n:
					num_rsa_id+=1
				else:
					
					others+=1

	print('%d singleton nodes (%.2f)' % (len(singleton_nodes),len(singleton_nodes)/i))
	print('%d nodes are members of entity sets or complexes but NOT entity sets themselves.' % (num_members))
	print('%d nodes have R-HSA ids and will be ignored (these are reactions, not entities' % (num_rsa_id))
	print('%d are "other"' % (others))
	if outfile:
		out = open(outfile,'w')
		for n in singleton_nodes:
			out.write('%s\n' % n)
		out.close()
		print('wrote to %s' % (outfile))
	return


def add_entity_set_info(H,prefix):
	## adds entity set information (like complexes)
	## hard -coded in for now. 
	files = glob.glob('%s*-entitysets.txt' % (prefix))
	es = {}
	for f in files:
		with open(f) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				row = line.strip().split()
				es[row[0]] = row[2].split(';')
	for node in H.get_node_set():
		if node in es:
			H.add_node(node,is_entityset=True,entityset_members=es[row[0]])
		else:
			H.add_node(node,is_entityset=False,entityset_members=[])
	return H

def to_digraph(H):
	## simple function so we convert to graph by calling an hgraph_utils function.
	G = transform.to_networkx_digraph(H)
	return G

def to_bipartite_graph(H):
	G = nx.DiGraph()

	for node in H.node_iterator():
		G.add_node(node, hedge_node = False)

	for hyperedge_id in H.hyperedge_id_iterator():
		# add hyperedge id node
		G.add_node(hyperedge_id, hedge_node = True)
		for tail_node in H.get_hyperedge_tail(hyperedge_id):
			G.add_edge(tail_node,hyperedge_id)
		for head_node in H.get_hyperedge_head(hyperedge_id):
			G.add_edge(hyperedge_id,head_node)
	return G

def make_b_visit_dict(hedge_connectivity_file,identifier2id):
	# make a dictionary where the key is a hyperedge ID and the value is a 
	# 3-tuple of (bconnected nodes, traversed hedges, restrictive hedges). 
	b_visit_dict = {}
	with open(hedge_connectivity_file) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split('\t')
			name = row[0]
			bconn = row[2].split(';')
			if len(row) > 3 and row[3] != '':
				traversed = row[3].split(';')
			else:
				traversed = []
			if len(row) > 4 and row[4] != '':
				restrictive = row[4].split(';')
			else:
				restrictive = []

			## convert REACTOME identifier (e.g. R-HSA-XXXX) to hyperedge ID.  Do this here
			## to make sure that hyperedges aren't re-ordered during multiple make_hypergraph calls.
			name = identifier2id[name]
			traversed = [identifier2id[n] for n in traversed]
			restrictive = [identifier2id[n] for n in restrictive]

			b_visit_dict[name] = (bconn,traversed,restrictive)

	#print('%d hyperedges processed' % (len(b_visit_dict)))
	return b_visit_dict

def filter_by_blacklisted_entities(H,blacklist_file,outfile=None):
	blacklisted = set()
	with open(blacklist_file) as fin:
		for line in fin: # take first column of blacklist file
			blacklisted.add(line.strip().split()[0])
	print('%d entities blacklisted' % (len(blacklisted)))

	## remove blacklisted entities, retaining hyperedges if they have nonempty
	## heads and tails. Can't use the halp remove_node() function because that
	## removes any hyperedge that contains the node.
	## also have to check hypernodes, removing them if all members are blacklisted.
	nodes_to_remove = set()
	for node in H.get_node_set():
		if node in blacklisted: # node is in blacklisted set
			nodes_to_remove.add(node)
		if H.get_node_attribute(node,'is_hypernode'):
			#print('checking hypernode %s' % (node))
			members = H.get_node_attribute(node,'hypernode_members')
			num_blacklisted = sum([e in blacklisted for e in members])
			if num_blacklisted == len(members):
				nodes_to_remove.add(node)
				print('adding hypernode to removed nodes.')
	print('set to remove %d nodes' % (len(nodes_to_remove)))
	if outfile:
		out = open(outfile,'w')
		for n in nodes_to_remove:
			out.write('%s\n'% (n))
		out.close()
		print('wrote to %s' % (outfile))

	return H


def get_id_map(prefix,common_name=False):
	## reads ALL of reactome pathways and generates a PC2 identifier to uniprot.  
	pcid2uid = {}
	uid2pcid = {}
	files = glob.glob('%s*-elements.txt' % (prefix))
	for f in files:
		with open(f) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				row = line.strip().split()
				if row[0] in pcid2uid:
					continue
				items = row[3].split(';')
				for i in items:
					if common_name:
						if 'hgnc-symbol' in i:
							#print(row[3].split(':'))
							hgnc = i.split(':')[1]
							pcid2uid[row[0]] = hgnc
							uid2pcid[hgnc] = row[0]
					else:
						if 'uniprot-knowledgebase' in i:
							#print(row[3].split(':'))
							uid = i.split(':')[1]
							pcid2uid[row[0]] = uid
							uid2pcid[uid] = row[0]
	if common_name:
		print('Read mappings for %d pathway commons identifiers to HGNC symbol (common name)' % (len(pcid2uid)))
	else:
		print('Read mappings for %d pathway commons identifiers to UniProt ID' % (len(pcid2uid)))
	return pcid2uid, uid2pcid