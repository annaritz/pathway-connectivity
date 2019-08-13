## Pathway Lists for Reactome

The list of pathways depends on the data representation.

### Pathways from Hypergraph
This script outputs `reactome-pathways-from-hypergraphs.txt` (in the repo).  This is used by the hypergraph, compound graph, and biparite grpah representations.

```
python3 get_pathways_from_hypergraphs.py
```

### Pathways from SIF-formatted Graph
This script outputs `reactome-pathways-from-SIF.txt` (in the repo). This is used by the graph representation (which has different molecules than the other representations, hence the difference in number of sources in the  connectivity surveys).

```
python3 get-uniprot-pathways.py 
```

### BioPAX Parser output (deprecated)
`reactome-pathways.txt` is from BioPAX Parser.  
