## run with 34 original pathways
python3 run-channels.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter ../hypergraph_code/output/small_molecule_filter-reactome_hedges.txt ../../hypergraph/reactome_hypergraphs_parsed/ small_molecule_filter > small_molecule_filter.out
python3 run-channels.py ../../hypergraph/reactome_hypergraph_full/reactome ../hypergraph_code/output/reactome_hedges.txt ../../hypergraph/reactome_hypergraphs_parsed/ full > full.out

python3 viz-channels.py small_molecule_filter
python3 viz-channels.py full

## run with ~200 second-level pathways
python3 run-channels.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter ../hypergraph_code/output/small_molecule_filter-reactome_hedges.txt ../../hypergraph/reactome_hypergraphs_parsed/ small_molecule_filter_allpathways run_all > small_molecule_filter_allpathways.out 
python3 viz-channels.py small_molecule_filter_allpathways

###########
# old
#python3 run-channels.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter ../hypergraph_code/output/small_molecule_filter-reactome_hedges.txt ../hypergraph_code/output/pathways/small_molecule_filter_
#python3 run-channels.py ../../hypergraph/reactome_hypergraph_full/reactome ../hypergraph_code/output/reactome_hedges.txt ../hypergraph_code/output/pathways/reactome_

