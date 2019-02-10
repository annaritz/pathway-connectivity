## computes the survey on SIF-formatted graphs.

import networkx as nx
import sys

def main(infile,outfile,conversions):
	## use a directed graph.
	G = nx.DiGraph()
	with open(infile) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			# SIF format is <n1 TYPE n2>
			row = line.strip().split('\t')
			n1 = row[0]
			ntype = row[1]
			n2 = row[2]

			if ntype not in conversions:
				print('ERROR:',ntype)
				sys.exit()

			if conversions[ntype] == 'IGNORE':
				continue
			elif conversions[ntype] == 'DIR':
				G.add_edge(n1,n2)
			elif conversions[ntype] == 'UNDIR':
				G.add_edge(n1,n2)
				G.add_edge(n2,n1)
			else:
				print('ERROR:',ntype,conversions[ntype])
				sys.exit()
	print('Graph has %d nodes and %d edges' % (nx.number_of_nodes(G),nx.number_of_edges(G)))
	nx.write_edgelist(G,outfile+'.graph',delimiter='\t',data=False)
	survey_graph(G,outfile)

	return

def survey_graph(G,outfile):
	out = open(outfile,'w')
	out.write('#Name\tNumConnected\n')
	i = 0
	for n in G.nodes():
		i+=1
		if i % 100 == 0:
			print(' %d of %d' % (i,nx.number_of_nodes(G)))
		#print(n)
		edges = nx.bfs_edges(G,n)
		nset = set()
		for u,v in edges:
			nset.add(u)
			nset.add(v)
		out.write('%s\t%d\n' % (n,len(nset)))

	out.close()
	print('wrote to %s' % (outfile))

	return

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print('USAGE: python3 GraphSurvey.py <infile> <outfile> <conversion-type-file>')
	infile = sys.argv[1]
	outfile = sys.argv[2]
	convfile = sys.argv[3]
	conversions = {}
	with open(convfile) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split('\t')
			conversions[row[0]] = row[1]
	print('%d conversions read' % (len(conversions)))
	main(infile,outfile,conversions)