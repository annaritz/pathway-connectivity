import sys

# hard-code files for now. 
prefix = 'hypergraph/Disassembly-of-the-destruction-complex-and-recruitment-of-AXIN-to-the-membrane'
hedge_file = prefix+'-hyperedges.txt'
hnode_file = prefix+'-hypernodes.txt'
sif_file = 'destruction-complex-disassembly.extended.sif'

outprefix = 'consistent-destruction-complex-assembly'

element_file = prefix+'-elements.txt'

def main():
	# read hyperedges
	print(hedge_file)
	hedges = open(hedge_file).readlines()
	hedges = hedges[1:] # skip header
	hedges = [line.strip().split() for line in hedges]
	for i in range(len(hedges)):
		for j in range(len(hedges[i])):
			hedges[i][j] = hedges[i][j].split(';')
	print(' %d hyperedges' % (len(hedges)))

	# read hypernodes
	print(hnode_file)
	hnodes = open(hnode_file).readlines()
	hnodes = hnodes[1:] # skip header
	hnodes = [line.strip().split() for line in hnodes]
	for i in range(len(hnodes)):
		hnodes[i][1] = hnodes[i][1].split(';')
	print(' %d hypernodes' % (len(hnodes)))

	# read SIF file
	print(sif_file)
	simple_interactions = open(sif_file).readlines()
	simple_interactions = [line.strip().split() for line in simple_interactions]
	print(' %d simple interactions' % (len(simple_interactions)))
	sif_types = {}
	for row in simple_interactions:
		if row[1] not in sif_types:
			sif_types[row[1]] = 0
		sif_types[row[1]]+=1
	for t in sif_types:
		print(' --> %d %s' % (sif_types[t],t))

	# read mapping files
	print(element_file)
	uniprot2id = {}
	id2uniprot = {}
	uniprot2name = {}
	id2name = {}
	with open(element_file) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split()
			this_id = row[0]
			this_name = row[1]
			id2name[this_id] = this_name
			this_uniprot = None
			for e in row[3].split(';'):
				if 'uniprot-knowledgebase' in e:
					this_uniprot = e.split('uniprot-knowledgebase:')[1]
					break
			if this_uniprot:
				if this_uniprot not in uniprot2id:
					uniprot2id[this_uniprot] = set()
				uniprot2id[this_uniprot].add(this_id)
				if this_id not in id2uniprot:
					id2uniprot[this_id] = set()
				id2uniprot[this_id].add(this_uniprot)
				if this_uniprot not in uniprot2name:
					uniprot2name[this_uniprot] = set()
				uniprot2name[this_uniprot].add(this_name)
	print(' --> %d uniprot, %d ids' % (len(uniprot2id),len(id2uniprot)))

	

	return

if __name__ == '__main__':
	main()