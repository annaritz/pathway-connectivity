
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
-- BEFORE ADDING 3 ADDITIONAL MOLECULES
##4975 reactions pruned to filter small molecules
##2778 small nodes in Reactome signailng pathways
-- AFTER ADDING 3 ADDITIONAL MOLECULES (Ub, Ub, NPC)
##5180 reactions pruned to filter small molecules
##2781 small nodes in Reactome signailng pathways

cp ../../hypergraph/reactome_hypergraph_full/reactome-hyperedges.txt.small_molecule_filter ../../hypergraph/reactome_hypergraph_full/small_molecule_filter-hyperedges.txt
cp ../../hypergraph/reactome_hypergraph_full/reactome-hypernodes.txt ../../hypergraph/reactome_hypergraph_full/small_molecule_filter-hypernodes.txt 

python connectivity-survey.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter output/small_molecule_filter-
#40 reactions skipped because of Reactome identifier
#2342 reactions skipped because of an empty tail or head
#Hypergraph has 8773 hyperedges and 15440 nodes (8327 of these are hypernodes)

python3 b-relaxation-survey.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter output/small_molecule_filter_b_relax.txt output/small_molecule_filter-reactome_hedges.txt

## PATHWAY SURVEY

## traversal-based pathways
#python3 b-relaxation-pathways.py ../../hypergraph/reactome_hypergraph_full/reactome output/pathways/full_reactome_ output/reactome_hedges.txt  ../../data/pathways/reactome-pathways.txt

#python3 b-relaxation-pathways.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter output/pathways/small_molecule_filter_ output/small_molecule_filter-reactome_hedges.txt  ../../data/pathways/reactome-pathways.txt

#python3 b-relaxation-pathways.py ../../hypergraph/reactome_hypergraph_full/blacklist_filter output/pathways/blacklist_filter_ output/blacklist_filter-reactome_hedges.txt  ../../data/pathways/reactome-pathways.txt

## hypergraph-based pathways (PARSING)
python3 b-relaxation-pathways.py ../../hypergraph/reactome_hypergraph_full/reactome output/pathways/full_reactome_ output/reactome_hedges.txt   ../../data/pathways/reactome-pathways-from-hypergraphs.txt

python3 b-relaxation-pathways.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter output/pathways/small_molecule_filter_ output/small_molecule_filter-reactome_hedges.txt  ../../data/pathways/reactome-pathways-from-hypergraphs.txt

### HUB STUFF (from connectivity-sruvey.py)
sort -grk3 output/reactome_hubs.txt | head -n 50 | grep -v SmallMolecule | cut -f 1 > entities_in_top_50_hubs.txt
grep -f entities_in_top_50_hubs.txt ../BioPAXSTREAM/output/reactome_limit20.txt.names
#Ub http://pathwaycommons.org/pc2/Protein_fa92e525bc3eb51cffd7786a1f19e317	104	136
#NPC http://pathwaycommons.org/pc2/Complex_6ae49f2fe344df4b6985f2f372910a77	35	81
#Ub http://pathwaycommons.org/pc2/Protein_80c9e4746b9a9261c9c7b174d2cf8292	31	73



## filter Graphw with Complexes by small molecules
python3 graph-with-complexes-survey.py ../../hypergraph/reactome_hypergraph_full/reactome output/graph-reactome.txt

########## PERMUTATION TEST ############

## 1000 swaps
python3 permutation-test.py output/pathways/small_molecule_filter_ permutations/small_molecule_filter_ 100 1000
python3 b-relaxation-permutation.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter output/pathways/permutations/small_molecule_filter_ output/small_molecule_filter-reactome_hedges.txt  permutations/small_molecule_filter_  100 1000

## 10000 swaps
python3 permutation-test.py output/pathways/small_molecule_filter_ permutations/small_molecule_filter_ 100 10000
python3 b-relaxation-permutation.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter output/pathways/permutations/small_molecule_filter_ output/small_molecule_filter-reactome_hedges.txt  permutations/small_molecule_filter_  100 10000

## 10000 swaps. 1000 iterations
python3 permutation-test.py output/pathways/small_molecule_filter_ permutations/small_molecule_filter_ 1000 10000
python3 b-relaxation-permutation.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter output/pathways/permutations/small_molecule_filter_ output/small_molecule_filter-reactome_hedges.txt  permutations/small_molecule_filter_  1000 10000