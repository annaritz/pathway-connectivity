## old
##python3 connectivity-survey.py ../SIF/outfiles/reactome.txt ../BioPAXSTREAM/output/reactome_limit20.txt ../hypergraph/output/reactome.txt ../hypergraph/output/reactome_hedges.txt connectivity-histogram

python3 connectivity-survey.py ../SIF/outfiles/reactome.txt ../hypergraph_code/output/graph-reactome.txt ../BioPAXSTREAM/output/reactome_limit30_filtered.txt ../hypergraph_code/output/reactome.txt connectivity-histogram


python3 brelax-survey.py ../hypergraph_code/output/b_relax_reactome.txt  ../hypergraph_code/output/blacklist_filter_b_relax.txt b_relax-histogram

python3 brelax-survey.py ../hypergraph_code/output/b_relax_reactome.txt  ../hypergraph_code/output/blacklist_filter_b_relax.txt ../hypergraph_code/output/small_molecule_filter_b_relax.txt b_relax-histogram_threepanel

python3 connectivity-transformations.py ../hypergraph/output/reactome.txt ../hypergraph_code/output/graph-reactome.txt transformed-histogram

python3 pathway-influence.py ../hypergraph_code/output/pathways/small_molecule_filter_ influence-output/small_molecule_filter