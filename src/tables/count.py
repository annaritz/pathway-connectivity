# count number of elements in each type of graph.  
import networkx as nx
import sys
import hgraph_utils
from halp import directed_hypergraph
from halp.directed_hypergraph import DirectedHypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform

GRAPH_FILE = '../SIF/outfiles/reactome.txt.graph'
HGRAPH_PREFIX = '../../hypergraph/reactome_hypergraph_full/reactome'
SMALLMOL_PREFIX = '../../hypergraph/reactome_hypergraph_full/small_molecule_filter'

hypergraph, identifier2id, id2identifier = hgraph_utils.make_hypergraph(HGRAPH_PREFIX)
hgraph_nodes = stats.number_of_nodes(hypergraph)
hgraph_edges = stats.number_of_hyperedges(hypergraph)

filtered_hypergraph, identifier2id, id2identifier = hgraph_utils.make_hypergraph(SMALLMOL_PREFIX)
filtered_nodes = stats.number_of_nodes(filtered_hypergraph)
filtered_edges = stats.number_of_hyperedges(filtered_hypergraph)

compound_graph = transform.to_networkx_digraph(hypergraph)
compound_nodes = nx.number_of_nodes(compound_graph)
compound_edges = nx.number_of_edges(compound_graph)

graph = nx.read_edgelist(GRAPH_FILE)
graph_nodes = nx.number_of_nodes(graph)
graph_edges = nx.number_of_edges(graph)

bipartite_nodes = hgraph_nodes+hgraph_edges
bipartite_edges = 0
for hedge_id in hypergraph.hyperedge_id_iterator():
	bipartite_edges += len(hypergraph.get_hyperedge_head(hedge_id))+len(hypergraph.get_hyperedge_tail(hedge_id))

print('\\begin{table}[h]')
print('\\centering')
print('\\begin{tabular}{|r||c|c|c|cc|} \\hline')
print('& Directed  & Bipartite  & Compound  & \\multicolumn{2}{c|}{Hypergraph} \\\\ ')
print('&  Graph &  Graph &  Graph & Full & Filtered \\\\ \\hline')
print('\\# Nodes & %d &  %d &  %d & %d & %d \\\\' % (graph_nodes,bipartite_nodes,compound_nodes,hgraph_nodes,filtered_nodes))
print('\\# Edges/Hyperedges & %d &  %d &  %d & %d & %d \\\\ \\hline' % (graph_edges,bipartite_edges,compound_edges,hgraph_edges,filtered_edges))
print('\\end{tabular}')
print('\\caption{\\textbf{Representations of the Reactome pathway database.}  The filtered hypergraph has removed all small molecules, two forms of Ubiquitinase, and the Nuclear Pore Complex from the hyperedges.}')
print('\\label{tab:representations}')
print('\\end{table}')