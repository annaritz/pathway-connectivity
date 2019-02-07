from halp import directed_hypergraph
from halp.utilities import directed_statistics as stats

def make_hypergraph(file_prefix,delim=';',sep='\t'):
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
	identifier2id = {}
	id2identifier = {}
	H = directed_hypergraph.DirectedHypergraph()	
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

		## only add node if it currently exists
		#if not H.has_node(node):
	#		continue

		if node in hypernodes and hypernodes[node] != [node]:
			H.add_node(node,hypernode_members=hypernodes[node],is_hypernode=True)
			num_hypernodes+=1
		else:
			H.add_node(node,is_hypernode=False,hypernode_members=[])

		H.add_node(node)

	print('Hypergraph has %d hyperedges and %d nodes (%d of these are hypernodes)' % 
		(stats.number_of_hyperedges(H),stats.number_of_nodes(H),num_hypernodes))

	return H, identifier2id, id2identifier

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

	print('%d hyperedges processed' % (len(b_visit_dict)))
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