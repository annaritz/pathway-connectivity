import hgraph_utils
from halp import directed_hypergraph
from halp.directed_hypergraph import DirectedHypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform
import sys
from viz_utils import *

HGRAPH_PREFIX = '../../hypergraph/reactome_hypergraph_full/reactome'
hypergraph, identifier2id, id2identifier = hgraph_utils.make_hypergraph(HGRAPH_PREFIX)
all_nodes = hypergraph.get_node_set()
SMALLMOL_PREFIX = '../../hypergraph/reactome_hypergraph_full/small_molecule_filter'
filtered_hypergraph, identifier2id, id2identifier = hgraph_utils.make_hypergraph(SMALLMOL_PREFIX)
filtered_nodes = filtered_hypergraph.get_node_set()

PATHWAY_FILE = '../../data/pathways/reactome-pathways-from-hypergraphs.txt'
pathway_info = {}
with(open(PATHWAY_FILE)) as fin:
	for line in fin:
		if line[0] == '#':
			continue
		row = line.strip().split('\t')
		rid = row[0].split('/')[-1]
		name = row[1]
		identifier = row[1].replace(' ','-').replace('/','-')
		#if identifier not in sorted_pathways:
	#		sys.exit(identifier)
		members = set(row[3].split(';'))
		pathway_info[identifier] = (rid,name,members)

print('\\begin{table}[h]')
print('\\centering')
print('\\begin{tabular}{|ll|ccc|} \\hline')
print('& & \\# in & \\# in & \\# in Filtered \\\\ ')
print(' Signaling Pathway & Reactome ID & Pathway & Hypergraph & Hypergraph \\\\ \\hline')
for pathway in sorted_pathways:
	rid,name,members = pathway_info[pathway]
	print(' %s & %s & %d & %d & %d \\\\' % (NAMES[pathway],rid,len(members),
		len(members.intersection(all_nodes)),len(members.intersection(filtered_nodes))))
print('\\hline')
print('\\end{tabular}')
print('\\caption{\\textbf{%d Reactome signaling pathways considered for the pathway influence analysis.}  Members that are not part of any hyperedge are ignored from the hypergraph.  The filtered hypergraph has removed all small molecules, two forms of Ubiquitinase, and the Nuclear Pore Complex from the hyperedges.}' % (len(sorted_pathways)))
print('\\label{tab:pathways}')
print('\\end{table}')