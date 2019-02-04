python3 connectivity-survey.py ../SIF/outfiles/reactome.txt ../BioPAXSTREAM/output/reactome_limit20.txt ../hypergraph/output/reactome.txt ../hypergraph/output/reactome_hedges.txt connectivity-histogram

python3 brelax-survey.py ../hypergraph/output/b_relax_reactome.txt b_relax-histogram

 python3 connectivity-transformations.py ../hypergraph/output/reactome.txt ../hypergraph/output/graph-reactome.txt transformed-histogram