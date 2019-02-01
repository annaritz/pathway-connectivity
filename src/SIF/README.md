To run on Wnt pathway:


```
python3 GraphSurvey.py ../../reactome-signaling-by-wnt.extended.sif outfiles/reactome-signaling-by-wnt.txt conversion-types.txt 
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