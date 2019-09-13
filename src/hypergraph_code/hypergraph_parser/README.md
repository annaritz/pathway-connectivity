## Hypergraph parser

This code parsers PathwayCommons v10 files (specifically Reactome .OWL files). It is built upon the Java BioPAX Parser called PaxTools. See the paper [Using Biological Pathway Data with Paxtools
](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003194) (Demir et al., PLOS Computational Biology, 2013). For more information, refere to the [PaxTools API](http://biopax.github.io/Paxtools/5.1.0/apidocs/), and the [Users Guide](https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=info:doi/10.1371/journal.pcbi.1003194.s001).

(Note that steps 1 and 2 are also from the BioPAXSTREAM/ directory; however they require different PaxTools JAR files.)

1. [Download the JAR file](https://www.biopax.org/Paxtools/) (paxtools-4.3.1-no-jena.jar is available on [SourceForge](https://sourceforge.net/projects/biopax/files//paxtools/))

2. Add jar to classpath
```
export CLASSPATH=$CLASSPATH:<PATH>/paxtools-4.3.1-no-jena.jar
```

3. compile Java files
```
javac *.java
```

## Run `ParameterizedStreamSurvey.java`

### Developer Notes

(This is copied/merged from Anna's `biopax-parsers` private git repo.  See `biopax-parsers/java/src/`.)


