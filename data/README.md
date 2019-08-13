
## Datasets

Reactome Pathway available via [PathwayCommons](http://www.pathwaycommons.org/) ([Cerami et al., NAR 2011](http://nar.oxfordjournals.org/content/early/2010/11/10/nar.gkq1039.abstract)).  PC v10 files available [here](http://www.pathwaycommons.org/archives/PC2/v10/) - all files downloaded with the last modified date of `2019-02-19`.

### OWL

The BioPAX OWL file of Reactome, available from PathwayCommons. Make an `OWL/` directory and download the BioPAX File:

```
wget http://www.pathwaycommons.org/archives/PC2/v10/PathwayCommons10.reactome.BIOPAX.owl.gz
gunzip PathwayCommons10.reactome.BIOPAX.owl.gz
```

### SIF

Simple Interaction Format (SIF) file of Reactome, available from PathwayCommons. The `SIF/` directory already contains the following file:

```
wget http://www.pathwaycommons.org/archives/PC2/v10/PathwayCommons10.reactome.hgnc.sif.gz
gunzip PathwayCommons10.reactome.hgnc.sif.gz
```

### Pathways

This directory contains list of Reactome pathways according to different reprentations. See the `README.md` file there for more information.

### STRING

This directory contains interactions from the [STRING database](https://string-db.org/) ([von Mering et al., NAR 2005](https://www.ncbi.nlm.nih.gov/pubmed/15608232)); ([Szklarczyk et al., NAR 2019](https://www.ncbi.nlm.nih.gov/pubmed/30476243)); see the `README.md` file there for more information.

### Blacklisted Molecules

The `blacklist.txt` file contains ubiquitous molecules from PathwayCommons:

```
wget http://www.pathwaycommons.org/archives/PC2/v10/blacklist.txt
```
