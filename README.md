# pathway-connectivity

## [Pathway Commons](https://www.pathwaycommons.org/)

PC v10 files available [here](http://www.pathwaycommons.org/archives/PC2/v10/).  Start by looking at 'Signaling by Wnt' Reactome files.  this will use the [RESTful webservice API](http://www.pathwaycommons.org/pc2/#get). In a browser, 

```
http://www.pathwaycommons.org/pc2/get?uri=http://identifiers.org/reactome/R-HSA-195721&subpw=true
```

where `subpw=true` will report all sub-pathways of the pathway.  Move this file to `reactome-signaling-by-wnt.owl` for now. Copied `paxtools.jar` from v10 files, to convert to all other possible files:

```
java -jar paxtools.jar toSIF reactome-signaling-by-wnt.owl reactome-signaling-by-wnt.extended -extended  -andSIF seqDb=uniprot 
java -jar paxtools.jar toSBGN reactome-signaling-by-wnt.owl reactome-signaling-by-wnt.sbgn
java -jar paxtools.jar toGSEA reactome-signaling-by-wnt.owl reactome-signaling-by-wnt.gmt uniprot -subPathways
```

Example of the inference rules for SIF are on the [PathwayCommons site](http://www.pathwaycommons.org/pc2/formats#sif_relations).

## Parsing Hypergraph

In `/Users/aritz/git/biopax-parsers/java/bin`:

```
java PathwayCommonsParser /Users/aritz/Documents/github/pathway-connectivity/reactome-signaling-by-wnt.owl /Users/aritz/Documents/github/pathway-connectivity/hypergraph/ http://identifiers.org/reactome/R-HSA-195721 verbose
```

Back in this directory:

```
python src/make-signaling-hypergraph.py hypergraph/ hypergraph/
```

Makes `Signaling-by-WNT-hyperedges.txt` and `Signaling-by-WNT-hypernodes.txt` files in the `hypergraph/` directory.

## Get example from the HALP poster

```
http://www.pathwaycommons.org/pc2/get?uri=http://identifiers.org/reactome/R-HSA-4641262&subpw=true
mv ~/Downloads/get destruction-complex-disassembly.owl
java -jar paxtools.jar toSIF destruction-complex-disassembly.owl destruction-complex-disassembly.extended -extended  -andSIF seqDb=uniprot 
cd /Users/aritz/git/biopax-parsers/java/bin
java PathwayCommonsParser /Users/aritz/Documents/github/pathway-connectivity/destruction-complex-disassembly.owl /Users/aritz/Documents/github/pathway-connectivity/hypergraph/ http://identifiers.org/reactome/R-HSA-4641262 verbose
cd /Users/aritz/Documents/github/pathway-connectivity/
python src/make-signaling-hypergraph.py hypergraph/ hypergraph/
```

Make the OWL and SIF formats consistent (remove CHEBI elements; remove things with non-UniProt)

```
python3 src/make-formats-consistent.py hypergraph/Disassembly-of-the-destruction-complex-and-recruitment-of-AXIN-to-the-membrane-hyperedges.txt hypergraph/Disassembly-of-the-destruction-complex-and-recruitment-of-AXIN-to-the-membrane-hypernodes.txt destruction-complex-disassembly.extended.sif consistent-destruction-complex-assembly
```

## Get single reaction

http://www.pathwaycommons.org/pc2/get?uri=http://identifiers.org/reactome/R-HSA-201717&subpw=true
mv ~/Downloads/get R-HSA-201717.owl
java -jar paxtools.jar toSIF R-HSA-201717.owl R-HSA-201717.extended -extended  -andSIF seqDb=uniprot 
cd /Users/aritz/git/biopax-parsers/java/bin
java PathwayCommonsParser /Users/aritz/Documents/github/pathway-connectivity/R-HSA-201717.owl /Users/aritz/Documents/github/pathway-connectivity/hypergraph/ http://identifiers.org/reactome/R-HSA-201717 verbose
cd /Users/aritz/Documents/github/pathway-connectivity/
python src/make-signaling-hypergraph.py hypergraph/ hypergraph/

## Receptors

Looking at Reactome `*-entitysets.txt` and `*-elements.txt` files; taking one ID for each protein/group.  Read through the description of the pathway (e.g. [Reactome's description of Wnt signaling](https://reactome.org/content/detail/R-HSA-195721)). Considering proteins at the cell membrane, if there are options.  Not considering ubiquitinated proteins/groups.

## Multiple Data Sources (in the future)

The `pathways.txt` file contains two parts: the pathways and then the pathways broken by data source.  First get the point when the pathways are broken down by data source:
```
$ grep -n 'URI' pathways.txt 
1:PATHWAY_URI	DISPLAY_NAME	DIRECT_SUB_PATHWAY_URIS	ALL_SUB_PATHWAY_URIS
53763:PATHWAY_URI	DATASOURCE	DISPLAY_NAME	ALL_NAMES	NUM_DIRECT_COMPONENT_OR_STEP_PROCESSES
```
Then count all the pathways listed:
```
$ tail -n +53763 pathways.txt | cut -f 2 | sort | uniq -c
   1 DATASOURCE
 242 HumanCyc
 774 INOH
 122 KEGG
  27 NetPath
 295 PANTHER
2216 Reactome
49006 SMPDB
 333 WikiPathways
 745 pid
```
I am interested in all pathway sources except for the [Small Molecule Pathway Database](http://smpdb.ca/) (SMPDB) and [Integrating Network Objects with Hierarchies](http://inoh.hgc.jp/), which has a broken link.