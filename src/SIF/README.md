To run on Wnt pathway:


```
python3 GraphSurvey.py ../../..//pathway-connectivity/data/SIF/reactome-signaling-by-wnt.extended.sif outfiles/reactome-signaling-by-wnt.txt conversion-types.txt 
```

Get PC v10 files available [here](http://www.pathwaycommons.org/archives/PC2/v10/) -- download /PathwayCommons10.reactome.hgnc.sif.gz (version 2018-05-07)

```
gunzip PathwayCommons10.reactome.hgnc.sif.gz 
cut -f 2 PathwayCommons10.reactome.hgnc.sif | sort | uniq -c
```

This gives us the number and type of relationships:

```
49435 catalysis-precedes
12219 chemical-affects
4914 consumption-controlled-by
3730 controls-expression-of
3184 controls-phosphorylation-of
5266 controls-production-of
106958 controls-state-change-of
5218 controls-transport-of
1882 controls-transport-of-chemical
140097 in-complex-with
 568 reacts-with
3271 used-to-produce
```

```
python3 GraphSurvey.py PathwayCommons10.reactome.hgnc.sif outfiles/reactome.txt conversion-types.txt 
```


## FILTER BY SMALL MOLECULES
python3 GraphSurvey.py PathwayCommons10.reactome.hgnc.sif outfiles/reactome_filtered.txt conversion-types.txt true

## to get proteins.txt from elements.txt file... in hypergraphs_parsed/
cat *-elements.txt | cut -f 4 | sort -u | grep hgnc | sed 's/.*symbol://g' | sed 's/;.*//g' | sort -u > proteins.txt

Graph has 12086 nodes and 444204 edges
Graph has 10225 nodes and 393980 edges

## 1861 small molecules removed.