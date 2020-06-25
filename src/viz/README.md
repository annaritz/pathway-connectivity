##CONNECTIVITY
## old
##python3 connectivity-survey.py ../SIF/outfiles/reactome.txt ../BioPAXSTREAM/output/reactome_limit20.txt ../hypergraph/output/reactome.txt ../hypergraph/output/reactome_hedges.txt connectivity-histogram

python3 connectivity-survey.py ../SIF/outfiles/reactome.txt ../hypergraph_code/output/graph-reactome.txt ../BioPAXSTREAM/output/reactome_limit30_filtered.txt ../hypergraph_code/output/reactome.txt connectivity-histogram

python3 connectivity-transformations.py ../hypergraph/output/reactome.txt ../hypergraph_code/output/graph-reactome.txt transformed-histogram

## HUB SURVEY
python3 hub-survey.py ../hypergraph_code/output/reactome_hubs.txt hub-histogram

## BRELAX
python3 brelax-survey.py ../hypergraph_code/output/b_relax_reactome.txt  ../hypergraph_code/output/ubiquitous_filter_b_relax.txt b_relax-histogram

python3 brelax-survey.py ../hypergraph_code/output/b_relax_reactome.txt  ../hypergraph_code/output/ubiquitous_filter_b_relax.txt ../hypergraph_code/output/small_molecule_filter_b_relax.txt b_relax-histogram_threepanel

## PATHWAY INFLUENCE
python3 pathway-influence.py ../hypergraph_code/output/pathways/full_reactome_ influence-output/full_reactome
python3 pathway-influence.py ../hypergraph_code/output/pathways/small_molecule_filter_ influence-output/small_molecule_filter
##python3 pathway-influence.py ../hypergraph_code/output/pathways/ubiquitous_filter_ influence-output/ubiquitous_filter

## PERMUTATION TEST!!
python3 significant-pathway-influence.py ../hypergraph_code/output/pathways/small_molecule_filter_ influence-output/significant_small_molecule_filter_ ../hypergraph_code/output/pathways/permutations/small_molecule_filter_ 100 1000

######
## COMBINED SCORES & SIGS - single plot
python3 significant-pathway-scores.py single influence-output/small_molecule_filter_k_03.txt influence-output/significant_small_molecule_filter_100_perms_10000_swaps_k_03.txt influence-output/circles_small_molecule_filter_100_perms_10000_swaps_k_03

python3 significant-pathway-scores.py summary influence-output/small_molecule_filter_k_ influence-output/significant_small_molecule_filter_100_perms_10000_swaps_k_ influence-output/circles_small_molecule_filter_100_perms_10000_swaps_summary

#### for graph-with-complexes
python3 pathway-influence.py ../hypergraph_code/output/pathways/small_molecule_filter_graph_with_complexes_ influence-output/small_molecule_filter_graph_with_complexes
python3 significant-pathway-influence.py ../hypergraph_code/output/pathways/small_molecule_filter_ influence-output/significant_small_molecule_filter_graph_with_complexes_ ../hypergraph_code/output/pathways/permutations/small_molecule_filter_graph_with_complexes_ 100 10000
python3 significant-pathway-scores.py single influence-output/small_molecule_filter_graph_with_complexes_k_03.txt influence-output/significant_small_molecule_filter_graph_with_complexes_100_perms_10000_swaps_k_03.txt influence-output/circles_small_molecule_filter_graph_with_complexes_100_perms_10000_swaps_k_03
python3 significant-pathway-scores.py summary influence-output/small_molecule_filter_graph_with_complexes_k_ influence-output/significant_small_molecule_filter_graph_with_complexes_100_perms_10000_swaps_k_ influence-output/circles_small_molecule_filter_graph_with_complexes_100_perms_10000_swaps_summary


#### for graphs
python3 pathway-influence.py ../SIF/outfiles/pathways/reactome_filtered_ influence-output/reactome_filtered_graph
python3 significant-pathway-influence.py ../SIF/outfiles/pathways/reactome_filtered_ influence-output/significant_reactome_filtered_graph ../SIF/outfiles/pathways/permutations/reactome_filtered_ 100 10000
python3 significant-pathway-scores.py summary influence-output/reactome_filtered_graph_k_ influence-output/significant_reactome_filtered_graph100_perms_10000_swaps_k_ influence-output/circles_reactome_filtered_graph_100_perms_10000_swaps_summary
python3 significant-pathway-scores.py single influence-output/reactome_filtered_graph_k_01.txt influence-output/significant_reactome_filtered_graph100_perms_10000_swaps_k_01.txt influence-output/circles_reactome_filtered_graph_100_perms_10000_swaps_k_01
python3 significant-pathway-scores.py single influence-output/reactome_filtered_graph_k_03.txt influence-output/significant_reactome_filtered_graph100_perms_10000_swaps_k_03.txt influence-output/circles_reactome_filtered_graph_100_perms_10000_swaps_k_03

######

## CASE STUDY
python3 case-study.py ../hypergraph_code/output/pathways/small_molecule_filter_ case-study-output/small_molecule_filter_
