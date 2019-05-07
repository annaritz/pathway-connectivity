import hgraph_utils
import sys

## get pathways with UniProt IDs

HYPERGRAPH_DIR = '../../hypergraph/reactome_hypergraphs/'
SIF_FILE = '../SIF/PathwayCommons10.reactome.hgnc.sif'


pcid2hgnc, hgnc2pcid = hgraph_utils.get_id_map_common_name(HYPERGRAPH_DIR)

orig_pathway_names = {}
with open('reactome-pathways-from-hypergraphs.txt') as fin:
	for line in fin:
		if line[0] == '#':
			continue
		row = line.strip().split('\t')
		orig_pathway_names[row[1]] = set()
		num_found = 0
		num_skipped = 0
		for item in row[3].split(';'):
			if item in pcid2hgnc:
				num_found+=1
				orig_pathway_names[row[1]].add(pcid2hgnc[item])
			else:
				num_skipped+=1
				#print('SKIPPED:',item)
		print('%s: %d found and %d skipped (skipped are entity sets or complexes, the members of which are recorded)' % (row[1],num_found,num_skipped))

sif_nodes = set()
with open(SIF_FILE) as fin:
	for line in fin:
		row = line.strip().split()
		sif_nodes.add(row[0])
		sif_nodes.add(row[2])
print('%d total SIF nodes' % (len(sif_nodes)))
print(list(sif_nodes)[:20])

out = open('reactome-pathways-from-SIF.txt','w')
out.write('#Pathway	Name	NumMembers	Members\n')
pathway_members = {}
for name in orig_pathway_names:
	identifier = name.replace(' ','-').replace('/','-')
	members = orig_pathway_names[name].intersection(sif_nodes)
	print('%s: %d -> %d' % (name,len(orig_pathway_names[name]),len(members)))
	if len(members) == 0:
		print(orig_pathway_names[name])
		print(members)
		sys.exit()
	out.write('%s\t%s\t%d\t%s\n' % (name,name,len(members),';'.join(members)))
out.close()
print('wrote to reactome-pathways-from-SIF.txt')
