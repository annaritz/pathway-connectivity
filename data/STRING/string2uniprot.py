import sys

def main(infile):
	string2uniprot,string2name = read_mapping_file()

	num_skipped = 0
	skipped_set = set()
	i=0
	outfile = infile+'.mapped'
	out = open(outfile,'w')
	with open(infile) as fin:
		for line in fin:
			if 'protein1' in line: # header
				out.write('#%s\tname1\tname2\n' % (line.strip()))
				continue
			row = line.strip().split()
			if row[0] in string2uniprot and row[1] in string2uniprot:
				out.write('%s\t%s\t%s\t%s\t%s\n' % (string2uniprot[row[0]],string2uniprot[row[1]],'\t'.join(row[2:]),string2name[row[0]],string2name[row[1]]))
				i+=1
			else:
				num_skipped+=1
				if row[0] not in string2uniprot:
					skipped_set.add(row[0])
				if row[1] not in string2uniprot:
					skipped_set.add(row[1])
				#sys.exit()


	print('%d written, %d interactions skipped, and %d nodes skipped' % (i,num_skipped,len(skipped_set)))
	for s in skipped_set:
		print('  ',s)
	print('wrote to %s' % (outfile))
	return

def read_mapping_file():

	mapfiles = ['human.uniprot_2_string.2018.tsv','human.uniprot_2_string.MISSING.tsv']
	string2uniprot = {}
	string2name = {}
	for mapfile in mapfiles:
		with open(mapfile) as fin:
			for line in fin:
				row = line.strip().split()
				string = row[2]
				item = row[1].split('|')
				uniprot = item[0]
				name = item[1]
				string2uniprot[string] = uniprot
				string2name[string] = name
		print('Mapped %d String identifiers' % (len(string2uniprot)))
		#print('Mapped %d String identifiers to names' % (len(string2name)))
	return string2uniprot,string2name

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print('USAGE: python3 string2uniprot.py <INFILE>')

	main(sys.argv[1])