## Posting graphs to GraphSpace

This code posts hypergraphs to [GraphSpace](http://graphspace.org/), a collaborative graph-sharing platform.  This is not integrated with `run.py`. If you want to post graphs to GraphSpace,

1. Make an account on [GraphSpace](http://graphspace.org/).

2. Install the [Python GraphSpace module](http://manual.graphspace.org/projects/graphspace-python/en/latest/).

3. Add your username and password to line 2 of `post_to_graphspace.py` (note, the default uploads the graph to `compbio@reed.edu` with password `compbio`).

4. Run the following code:
```
python post_to_graphspace.py ../../hypergraph/reactome_hypergraph_full/small_molecule_filter ../hypergraph_code/output/pathways/ outfile
```