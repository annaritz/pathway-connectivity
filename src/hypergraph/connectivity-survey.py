import sys
import itertools
from halp import directed_hypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform

def main(prefix,outprefix):
	H = make_hypergraph(prefix)
	#survey_nodes(H,outprefix+'reactome.txt')
	survey_hedges(H,outprefix+'reactome_hedges.txt')
	return

def survey_nodes(H,outfile):
	out = open(outfile,'w')
	out.write('#Name\tNumConnected\n')
	i = 0
	for node in H.node_iterator():
		i+=1
		if i % 100 == 0:
			print('node %d of %d' % (i,stats.number_of_nodes(H)))

		bconn, ignore, ignore, ignore = hpaths.b_visit(H,node)
		out.write('%s\t%d\n' % (node,len(bconn)))
	out.close()
	print('wrote to %s' % (outfile))
	return

def survey_hedges(H,outfile):
	out1 = open(outfile,'w')
	out1.write('#Name\tNumBVisit\tBVisitNodes\n')
	i = 0
	H.add_node('SUPERSOURCE')
	for hedge_id in H.hyperedge_id_iterator():
		i+=1
		if i % 100 == 0:
			print('hyperedge %d of %d' % (i,stats.number_of_hyperedges(H)))
		hid = H.add_hyperedge('SUPERSOURCE',H.get_hyperedge_tail(hedge_id))
		bconn, ignore, ignore, ignore = hpaths.b_visit(H,'SUPERSOURCE')
		H.remove_hyperedge(hid)

		out1.write('%s\t%d\t%s\n' % (hedge_id,len(bconn),';'.join([n for n in bconn])))
	out1.close()
	print('wrote to %s' % (outfile))
	return

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
	
	H = directed_hypergraph.DirectedHypergraph()	
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

			if len(tail) > 0 and len(head) > 0:
				H.add_hyperedge(tail,head,identifier=hedge_id)
			else:
				print('WARNING: skipping hyperedge',hedge_id)

	## annotate nodes
	num_hypernodes = 0
	for node in H.get_node_set():
		if node in hypernodes and hypernodes[node] != [node]:
			H.add_node(node,hypernode_members=hypernodes[node],is_hypernode=True)
			num_hypernodes+=1
		else:
			H.add_node(node,is_hypernode=False)

		H.add_node(node)

	print('Hypergraph has %d hyperedges and %d nodes (%d of these are hypernodes)' % 
		(stats.number_of_hyperedges(H),stats.number_of_nodes(H),num_hypernodes))

	return H


if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2])