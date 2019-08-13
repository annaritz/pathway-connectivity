## Hypergraph code and utility functions

These functions are used exclusively in the `run.py` script.  These all rely on `halp`.

1. `hgraph_utils.py`

Utilities for reading, writing, and maninpulating hypergraphs.

2. `permutation_test.py`

Functions for running permutation tests on hypergraphs.

3. `ILP.shortest_hyperpath.py`

Functions for computing the shortest s-t B-hyperpath (from [Ritz et al., TCBB 2017](https://www.ncbi.nlm.nih.gov/pubmed/28991726)).  This requires `cplex`.