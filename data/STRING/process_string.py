import sys

#actions_file = '9606.protein.actions.v11.0.txt'
links_file = '9606.protein.links.detailed.v11.0.txt'
file_prefix = 'processed/'
map_file = 'uniprot_to_string.tab'

def main():

	string2uni = read_map_file()

	headers = []
	files = []
	counts = []
	with open(links_file) as fin:
		for line in fin:
			if len(headers) == 0:
				headers = line.strip().split()[2:]
				counts = [0]*len(headers)
				for h in headers:
					files.append(open('%s/%s.txt' % (file_prefix,h),'w'))
				continue
			row = line.strip().split()
			protein1 = row[0]
			protein2 = row[1]
			for i in range(len(files)):
				if int(row[i+2]) > 0:
					counts[i]+=1
					files[i].write('%s\t%s\t%s\t%s\t%s\n' % (protein1,protein2,string2uni.get(protein1,'NA'),string2uni.get(protein2,'NA'),row[i+2]))
	for f in files:
		f.close()

	for i in range(len(headers)):
		print(headers[i],counts[i])


def read_map_file():
	string2uni = {}
	start = True
	with open(map_file) as fin:
		for line in fin:
			if start:
				start = False
				continue # skip header
			row = line.strip().split()
			# take first string ID
			uniprot_id = row[0]
			string_id = row[2].split(';')[0]
			string2uni[string_id] = uniprot_id
	return string2uni



if __name__ == '__main__':
	main()