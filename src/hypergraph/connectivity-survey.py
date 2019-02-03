import sys
import itertools
from halp import directed_hypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform
import hgraph_utils

def main(prefix,outprefix):
	H = hgraph_utils.make_hypergraph(prefix)
	print('skipping both. uncomment one of the surveys: nodes or hedges.')
	#survey_nodes(H,outprefix+'reactome.txt')
	#survey_hedges(H,outprefix+'reactome_hedges.txt')

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
	out1.write('#Name\tNumBVisit\tBVisitNodes\tTraversedHedges\tRestrictiveHedges\n')
	i = 0
	for hedge_id in H.hyperedge_id_iterator():
		i+=1
		if i % 100 == 0:
			print('hyperedge %d of %d' % (i,stats.number_of_hyperedges(H)))
		bconn, traversed, restrictive = hpaths.b_visit_restrictive(H,H.get_hyperedge_tail(hedge_id))
		out1.write('%s\t%d\t%s\t%s\t%s\n' % \
			(hedge_id,len(bconn),
				';'.join([n for n in bconn]),
				';'.join([e for e in traversed]),
				';'.join([e for e in restrictive])))
	out1.close()
	print('wrote to %s' % (outfile))
	return




if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2])