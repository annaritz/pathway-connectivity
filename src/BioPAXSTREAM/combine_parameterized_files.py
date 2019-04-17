# number of distances processed (eventually will be 30)
max_val = 30

## get ID to CommonNAme mapping
namefile = 'output/reactome_limit1_filtered.txt.names'
names = {}
with open(namefile) as fin:
	for line in fin:
		row = line.strip().split()
		names[row[0]] = row[1]

## initialize data dictionary to be keys of common names and 
## lists sized by the processed distances.
data = {}
for name in names.values():
	data[name] = [0]*max_val

## for every distance, populate each source with the number of items within that distance
## of the source. Cumulative
for i in range(1,max_val+1):
	infile = 'output/reactome_limit%d_filtered.txt' % (i)
	with open(infile) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split()
			name = names[row[0]]
			data[name][i-1] = row[1]

## write values to file
outfile = 'output/reactome_parameterized_filtered.txt' 
out = open(outfile,'w')
out.write('#Name\td=1\td=2\td=3\t...\n')
for name in data:
	out.write('%s\t%s\n' % (name,'\t'.join(data[name])))
out.close()
print('wrote to '+outfile)
