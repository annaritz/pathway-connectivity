
python connectivity-survey.py ../../hypergraph/reactome_hypergraph_full/reactome output/
#40 reactions skipped because of Reactome identifier
#0 reactions skipped because of an empty tail or head
#Hypergraph has 11125 hyperedges and 19650 nodes (19650 of these are hypernodes)

python3 b-relaxation-survey.py ../../hypergraph/reactome_hypergraph_full/reactome output/b_relax_reactome.txt output/reactome_hedges.txt 

python3 graph-with-complexes-survey.py ../../hypergraph/reactome_hypergraph_full/reactome output/graph-reactome.txt

## BLACKLIST STUFF
python filter-signaling-hypergraph.py ../../hypergraph/reactome_hypergraph_full/reactome-hyperedges.txt ../../data/blacklist.txt
#3286 reactions pruned to filter blacklisted molecules
#155 blacklisted nodes in Reactome signailng pathways

cp ../../hypergraph/reactome_hypergraph_full/reactome-hyperedges.txt.blacklist_filter ../../hypergraph/reactome_hypergraph_full/blacklist_filter-hyperedges.txt
cp ../../hypergraph/reactome_hypergraph_full/reactome-hypernodes.txt ../../hypergraph/reactome_hypergraph_full/blacklist_filter-hypernodes.txt 

python connectivity-survey.py ../../hypergraph/reactome_hypergraph_full/blacklist_filter output/blacklist_filter-
#40 reactions skipped because of Reactome identifier
#39 reactions skipped because of an empty tail or head
#Hypergraph has 11085 hyperedges and 19469 nodes (8904 of these are hypernodes)

python3 b-relaxation-survey.py ../../hypergraph/reactome_hypergraph_full/blacklist_filter output/blacklist_filter_b_relax.txt output/blacklist_filter-reactome_hedges.txt

## SMALL MOLECULE STUFF
python filter-signaling-hypergraph.py ../../hypergraph/reactome_hypergraph_full/reactome-hyperedges.txt 
##4975 reactions pruned to filter small molecules
##2778 small nodes in Reactome signailng pathways

cp ../../hypergraph/reactome_hypergraph_full/reactome-hyperedges.txt.small_molecule_filter ../../hypergraph/reactome_hypergraph_full/small_molecule_filter-hyperedges.txt
cp ../../hypergraph/reactome_hypergraph_full/reactome-hypernodes.txt ../../hypergraph/reactome_hypergraph_full/small_molecule_filter-hypernodes.txt 

python connectivity-survey.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter output/small_molecule_filter-
#40 reactions skipped because of Reactome identifier
#2326 reactions skipped because of an empty tail or head
#Hypergraph has 8790 hyperedges and 15446 nodes (8328 of these are hypernodes)

python3 b-relaxation-survey.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter output/small_molecule_filter_b_relax.txt output/small_molecule_filter-reactome_hedges.txt

## PATHWAY SURVEY
python3 b-relaxation-pathways.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter output/pathways/small_molecule_filter_ output/small_molecule_filter-reactome_hedges.txt  ../../data/pathways/reactome-pathways.txt


python3 b-relaxation-pathways.py ../../hypergraph/reactome_hypergraph_full/blacklist_filter output/pathways/blacklist_filter_ output/blacklist_filter-reactome_hedges.txt  ../../data/pathways/reactome-pathways.txt