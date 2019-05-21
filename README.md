## Connectivity Measures for Signaling Pathway Topologies

Data and source code for (nearly) all figures and experiments from 

Connectivity Measures for Signaling Pathway Topologies
Franzese, Groce, Murali, and Ritz
[bioRxiv 2019](https://doi.org/10.1101/593913](https://doi.org/10.1101/593913)

_May 2019:_ This repo is still a little rough, Anna plans to refactor code and update instructions for data acquisition and processing in the next few weeks. Email her if you have questions.

### Dependencies

- Python (tested with v3.5.1)
- [Hypergraph Algorithms Library in Python (HALP)](http://murali-group.github.io/halp/)
- [NetworkX](https://networkx.github.io/) (v1.11, should work with v2.1)
- [IBM ILOG CPLEX optimizer Python module](https://www.ibm.com/analytics/cplex-optimizer) (for shortest B-hyperpath - academic license available)

### Datasets

Reactome Pathway available via [PathwayCommons](http://www.pathwaycommons.org/).  PC v10 files available [here](http://www.pathwaycommons.org/archives/PC2/v10/).  Download the SIF and BioPAX formats of Reactome. More details coming...

### 
```
python3 run.py -h
usage: run.py [-h] [--force] [--printonly] [--keep_singletons]
              [--small_molecule_filter] [--blacklist_filter] [--stats]
              [--histograms] [--perm_test #PERMS] [--set_seed]
              [--case_studies] [--string_channels]

Run experiments for pathway connectivity. The four represenations are "SIF-
Graph","Compound Graph","Bipartite Graph", and "Hypergraph"

optional arguments:
  -h, --help            show this help message and exit

General Arguments:
  --force               force existing files to be overwritten (default=False)
  --printonly           print the commands instead of running them
                        (default=False)

Dataset Arguments:
  --keep_singletons     Keep singleton nodes. Default False.
  --small_molecule_filter
                        Filter by small molecule nodes
  --blacklist_filter    Filter by PathwayCommons blacklist filtered nodes

Experiment Arguments:
  --stats               print statistics about each representation.
  --histograms          print histograms and heatmaps for each representation.
  --perm_test #PERMS    pathway influence permutation test with #PERMS number
                        of permutations.
  --set_seed            Set the seed of the random number generator to 123456.
                        Default False.
  --case_studies        pathway influence case studies (hard-coded)
  --string_channels     run STRING channel assessment.
  ```