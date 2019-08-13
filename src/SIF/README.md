
## Get SIF Conversion Types
See details about downloading the SIF-formatted graph in `data/` directory.  Here, we get the conversion types:

```
cut -f 2 ../../data/SIF/PathwayCommons10.reactome.hgnc.sif | sort | uniq -c
```

This gives us the number and type of relationships (May 2019 update):

```
263342 catalysis-precedes
15473 chemical-affects
8079 consumption-controlled-by
3733 controls-expression-of
3184 controls-phosphorylation-of
8052 controls-production-of
113558 controls-state-change-of
5232 controls-transport-of
5105 controls-transport-of-chemical
140302 in-complex-with
1922 reacts-with
5888 used-to-produce
```

This is used to make a `conversion-types.txt` file that manually determines which types are directed, undirected, or ignored.  

## Run Graph Survey on the SIF-Formatted graph
```
mkdir outfiles
python3 GraphSurvey.py ../../data/SIF/PathwayCommons10.reactome.hgnc.sif outfiles/reactome.txt conversion-types.txt 
```

We also filter by small molecules:

```
python3 GraphSurvey.py ../../data/SIF/PathwayCommons10.reactome.hgnc.sif outfiles/reactome_filtered.txt conversion-types.txt true
```

Check how many small molecules were removed (in `hypergraph/reactome_hypergraphs/`):

```
cat *-elements.txt | cut -f 4 | sort -u | grep hgnc | sed 's/.*symbol://g' | sed 's/;.*//g' | sort -u > proteins.txt
```
Use `proteins.txt` to determine the number of small proteins removed:

```
Graph has 12086 nodes and 444204 edges
Graph has 10225 nodes and 393980 edges
1861 small molecules removed.
```

## Graph Utils is used in `run.py`
