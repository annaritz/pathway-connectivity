## takes hyperedges file and a blacklist file and returns a hyperedges filtered file.

import sys
ADDITIONAL_MOLECULES = set([
	'http://pathwaycommons.org/pc2/Protein_fa92e525bc3eb51cffd7786a1f19e317', # Ub
	'http://pathwaycommons.org/pc2/Complex_6ae49f2fe344df4b6985f2f372910a77', # Nuclear Pore Complex
	'http://pathwaycommons.org/pc2/Protein_80c9e4746b9a9261c9c7b174d2cf8292', # Ub
	])

def main(infile,blacklist_file):
	blacklisted = set()
	with open(blacklist_file) as fin:
		for line in fin: # take first column of blacklist file
			blacklisted.add(line.strip().split()[0])
	print('%d entities blacklisted' % (len(blacklisted)))

	## NOTE THAT this only checks NODES (not hypernodes). However, another
	## check confirmed that there are no hypernodes that consist solely of
	## blacklisted molecules.
	out = open(infile+'.blacklist_filter','w')
	num_altered = 0
	blacklisted_nodes = set()
	with open(infile) as fin:
		for line in fin:
			if line[0] == '#':
				out.write(line[0])
			else:
				affected = False
				row = line.strip().split('\t') 
				#print(row)
				for i in range(len(row)):
					if row[i] == 'None':
						continue
					items = row[i].split(';')
					newitems = [item for item in items if item not in blacklisted]
					blacklisted_nodes.update(set([item for item in items if item in blacklisted]))
					if len(items) != len(newitems):
						affected = True
					if len(newitems) == 0:
						row[i] = 'None'
					else:
						row[i] = ';'.join(newitems)
				out.write('\t'.join(row)+'\n')
				if affected:
					num_altered+=1
	out.close()
	print('%d reactions pruned to filter blacklisted molecules' % (num_altered))
	print('%d blacklisted nodes in Reactome signailng pathways' % (len(blacklisted_nodes)))
	print('wrote to '+infile+'.blacklist_filter')
	return


def small_mols(infile):

	## NOTE THAT this only checks NODES (not hypernodes). However, another
	## check confirmed that there are no hypernodes that consist solely of
	## blacklisted molecules.
	out = open(infile+'.small_molecule_filter','w')
	num_altered = 0
	blacklisted_nodes = set([m for m in ADDITIONAL_MOLECULES])
	with open(infile) as fin:
		for line in fin:
			if line[0] == '#':
				out.write(line[0])
			else:
				affected = False
				row = line.strip().split('\t') 
				#print(row)
				for i in range(len(row)):
					if row[i] == 'None':
						continue
					items = row[i].split(';')
					newitems = [item for item in items if 'SmallMolecule' not in item and item not in ADDITIONAL_MOLECULES]
					blacklisted_nodes.update(set([item for item in items if'SmallMolecule' in item]))
					if len(items) != len(newitems):
						affected = True
					if len(newitems) == 0:
						row[i] = 'None'
					else:
						row[i] = ';'.join(newitems)
				out.write('\t'.join(row)+'\n')
				if affected:
					num_altered+=1
	out.close()
	print('%d reactions pruned to filter small molecules' % (num_altered))
	print('%d small nodes in Reactome signailng pathways' % (len(blacklisted_nodes)))
	print('wrote to '+infile+'.small_molecule_filter')
	return

if __name__ == '__main__':
	if len(sys.argv) != 3 and len(sys.argv) != 2:
		print('USAGE: python filter-signaling-hypergraph.py <HEDGES FILE> <BLACKLIST FILE - OPTIONAL>')

	if len(sys.argv) == 3:
		main(sys.argv[1],sys.argv[2])
	else:
		small_mols(sys.argv[1])