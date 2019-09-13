## Hypergraph parser

This code parsers PathwayCommons v10 files (specifically Reactome .OWL files). It is built upon the Java BioPAX Parser called PaxTools. See the paper [Using Biological Pathway Data with Paxtools
](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003194) (Demir et al., PLOS Computational Biology, 2013). For more information, refere to the [PaxTools API](http://biopax.github.io/Paxtools/5.1.0/apidocs/), and the [Users Guide](https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=info:doi/10.1371/journal.pcbi.1003194.s001).

(Note that steps 1 and 2 are also from the BioPAXSTREAM/ directory; however they require different PaxTools JAR files.)

1. [Download the JAR file](https://www.biopax.org/Paxtools/) (paxtools-4.3.1-no-jena.jar is available on [SourceForge](https://sourceforge.net/projects/biopax/files//paxtools/))

2. Add jar to classpath
```
export CLASSPATH=$CLASSPATH:<PATH>/paxtools-4.3.1-no-jena.jar
```

3. compile the source files in this directory.
```
javac *.java
```

This compiles three java files:
- `BioPaxParser.java`: the main program
- `BioPaxVisitor.java`: implements the BioPaxVisitor interface provided by paxtools. This traverses files and stores relevant information along the way.
- `DB.java`: an `enum` that provides different databases that are implemented.  Currently, only the `PathwayCommonsParser` is provided. A long-term goal is to release these parsers (or determine that PathwayCommons is all you need to pull relevant databases). 

## Usage Information for Parsing PathwayCommons Files


Running `java PathwayCommonsParser` without arguments give the following usage information:

```
java PathwayCommonsParser <Reactome.owl> <outprefix> [<pathway-id>] [<verbose>]
	<Reactome.owl>: BioPAX Level 3 file (.owl format).
	<outprefix>: Output prefix of all signaling pathways.
	<pathway-id> (optional): if specified, only runs that pathway id.
	<verbose> (optional): the word 'verbose' to print additional information to screen. Optional.

Note: to run from Java source, PathwayCommonsParser requires paxtools-4.2.1-no-jena.jar in your CLASSPATH. If running from the JAR file (using java -jar), it should not be required
```

### Example Run

Before running this parser, you need to (a) download the PathwayCommons OWL file (see the [README in the `data/` directory](https://github.com/annaritz/pathway-connectivity/tree/master/data)) and (b) made an empty directory to dump parsed interactions (I used `hypergraph/` at the of this repo).  

To print _all_ pathways in _verbose_ mode, call the method like the following:

```
java PathwayCommonsParser ../../../data/OWL/PathwayCommons10.reactome.BIOPAX.owl ../../../hypergraph/ verbose
```

This will output a collection of files _for each pathway_:
- `*-elements.txt`: all bioentities that have a single interpretation from this pathway.
- `*-entitysets.txt`: bioentities that may represent a family of interchangeable entities or members of a protein complex. 
- `*-subpathways.txt`: sub-pathways that are considered "nodes" in the pathway.
- `*-reactions.txt`: biochemical reactions for the pathway.
- `*-controls.txt`: controlling bioentities for the pathway (that reference the reaction IDs that they control).

The `verbose` option collects the number of objects and object types that are traversed, as well as objects that were ignored during the file traversal.

## Usage Information for Converting Parsed Files into Hypergraphs

The python script `make-signaling-hypergraph.py` converts the flat files above into a hyperedge file and a hypernode file for each pathway.  Running `python make-signaling-hypergraph.py` prints the following usage information:

```
USAGE: make-signaling-hypergraph.py <DATADIRECTORY> <OUTDIR>
	<DATADIRECTORY>: directory containing the files parsed by BioPAX Paxtools parser called ReactomeParser.jar
	<OUTDIR>: output directory for signaling hyergraph files.
```

### Example Run

I made a new directory for these files: (`../../../hypergraph/reactome_hypergraphs_parsed/`), then called the program (python2.7):

```
python make-signaling-hypergraph.py ../../../hypergraph/ ../../../hypergraph/reactome_hypergraphs_parsed/
```

This script also reports the number of reactions written, the number of reactions skipped, adn the number of hypernodes.

### Developer Notes

(This is copied/merged from Anna's `biopax-parsers` private git repo.  See `biopax-parsers/java/src/`.)


