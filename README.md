## Hypergraph-based connectivity measures for signaling pathway topologies

Data and source code for (nearly) all figures and experiments from a paper on [bioRxiv](https://doi.org/10.1101/593913):

**Hypergraph-based connectivity measures for signaling pathway topologies**  
Nick Franzese, Adam Groce, T.M. Murali, and Anna Ritz  

This work was presented at GLBio 2019 (accepted paper) ISMB 2019 (accepted poster), and is under review at PLOS Computational Biology.

_Aug 2019:_ Full dependencies, code for generating _almost_ all pathway representations (hypergraphs, bipartite graphs, compound graphs, and directed graphs) are now here.  Third-party datasets are also documented.  Email Anna with any questions.

TODO: dockerize this whole thing.

### Dependencies

Most of the experiment options require these tools:
- Python (tested with v3.5.1)
- [Hypergraph Algorithms Library in Python (HALP)](http://murali-group.github.io/halp/)
  - Check out `annabranch` (this passes all `pytest` tests; there are currently some issues with Travis CI)
- [NetworkX](https://networkx.github.io/) (v1.11, should work with v2.1)
- [xlsxwriter](https://xlsxwriter.readthedocs.io/) to write permutation test information to an Excel file.

If you want to implement shortest B-hyperpath on directed graphs, it is implemented as a `cplex` integer linear program.
- [IBM ILOG CPLEX optimizer Python module](https://www.ibm.com/analytics/cplex-optimizer) (academic license available)

### Datasets

The data is generated from the Reactome Pathway, available via [PathwayCommons](http://www.pathwaycommons.org/).  PC v10 files available [here](http://www.pathwaycommons.org/archives/PC2/v10/).  See details in the `data/` directory for more information.

To run all analyses, you must generate the following datasets:
- Simple Interaction Format (SIF, instructions in `data/`)
- Hypergraph, Bipartite Graph, and Compound Graph (OWL, instructions in `data/`)

The four data representations (hypergraph, bipartite graph, directed graph, and compound graph) are available in this repo.  To re-generate all representations, you can run the scripts outlined in the `src/` directory.  

### Usage Information

The main function is `run.py` in the `src/` directory:

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
Refer to the `README.md` files for more commands to generate data, run experiments, etc.
