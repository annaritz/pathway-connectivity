
## Datasets

Reactome Pathway available via [PathwayCommons](http://www.pathwaycommons.org/).  PC v10 files available [here](http://www.pathwaycommons.org/archives/PC2/v10/).

### OWL

The BioPAX OWL file of Reactome, available from PathwayCommons. In the `OWL/` directory:

```
wget http://www.pathwaycommons.org/archives/PC2/v10/PathwayCommons10.reactome.BIOPAX.owl.gz
unzip PathwayCommons10.reactome.BIOPAX.owl.gz
```


### SIF

Simple Interaction Format (SIF) file of Reactome, available from PathwayCommons. In the `SIF/` directory:

```
wget http://www.pathwaycommons.org/archives/PC2/v10/PathwayCommons10.reactome.hgnc.sif.gz
unzip PathwayCommons10.reactome.hgnc.sif.gz
```

### Pathways

This directory contains list of Reactome pathways according to different reprentations. See the `README.md` file there for more information.

### STRING

This directory contains interactions from the [STRING database](https://string-db.org/); see the `README.md` file there for more information.

### Blacklisted Molecules

The `blacklist.txt` file contains ubiquitous molecules from PathwayCommons:

```
wget http://www.pathwaycommons.org/archives/PC2/v10/blacklist.txt
```
