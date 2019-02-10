HYPERGRAPH_DIR = '../../hypergraph/reactome_hypergraphs_parsed/'

# sorted_pathways = ['Signaling-by-EGFR', 'Signaling-by-ERBB2', 'Signaling-by-ERBB4', 'Signaling-by-SCF-KIT', 
# 'Signaling-by-FGFR', 'ERK1-ERK2-pathway','Signaling-by-GPCR','Signaling-by-PDGF', 
# 'Signaling-by-VEGF',  'DAG-and-IP3-signaling', 
# 'Signaling-by-NTRKs',  'PI3K-AKT-Signaling', 'Signaling-by-WNT', 
# 'Integrin-signaling',
# 'Signaling-by-MET', 'Signaling-by-Type-1-Insulin-like-Growth-Factor-1-Receptor-(IGF1R)', 
# 'Signaling-by-Insulin-receptor', 
# 'TNF-signaling', #'TRAIL-signaling',  'FasL--CD95L-signaling', 
# 'Signaling-by-Activin', 'Signaling-by-TGF-beta-Receptor-Complex', 'Signaling-by-NOTCH', 'Signaling-by-PTK6', 
# 'Signaling-by-Rho-GTPases',  'MAPK6-MAPK4-signaling', 
#  'p75-NTR-receptor-mediated-signalling', 'Signaling-by-MST1',
# 'mTOR-signalling', 'Signaling-by-Hedgehog',
# 'Signaling-by-Nuclear-Receptors', 'Signaling-by-Leptin', 'Signaling-by-BMP', 'Signaling-by-Hippo']

orig_pathway_names = {}
with open('reactome-pathways.txt') as fin:
	for line in fin:
		if line[0] == '#':
			continue
		row = line.strip().split('\t')
		orig_pathway_names[row[1]] = row[0]

out = open('reactome-pathways-from-hypergraphs.txt','w')
out.write('#Pathway	Name	NumMembers	Members\n')
pathway_members = {}
for name in orig_pathway_names:
	identifier = name.replace(' ','-').replace('/','-')
	members = set()
	with open('%s/%s-hypernodes.txt' % (HYPERGRAPH_DIR,identifier)) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			members.add(line.strip().split()[0])
	pathway_members[identifier] = members
	out.write('%s\t%s\t%d\t%s\n' % (orig_pathway_names[name],name,len(members),';'.join(members)))
out.close()

for p1 in pathway_members:
	for p2 in pathway_members:
		if p1 == p2:
			continue
		if 'Leptin' in p1 and 'IGF' in p2:
			print('**')
		print(p1,p2,len(pathway_members[p1]),len(pathway_members[p2]),len(pathway_members[p1].intersection(pathway_members[p2])),len(pathway_members[p1].intersection(pathway_members[p2]))/len(pathway_members[p2]))
print('DONE. Wrote to reactome-pathways-from-hypergraphs.txt')
