
This is the BioPAXSTREAM code from [Using Biological Pathway Data with Paxtools
](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003194) (Demir et al., PLOS Computational Biology, 2013). See the [PaxTools API](http://biopax.github.io/Paxtools/5.1.0/apidocs/), and the [Users Guide](https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=info:doi/10.1371/journal.pcbi.1003194.s001) for more information.


1. [Download the JAR file](https://www.biopax.org/Paxtools/) (Paxtools-5.1.0.jar is available on [SourceForge](https://sourceforge.net/projects/biopax/files//paxtools/))

2. Add jar to classpath
```
export CLASSPATH=$CLASSPATH:<PATH>/paxtools-5.1.0.jar
```

3. compile Java files
```
javac *.java
```

## Run `ParameterizedStreamSurvey.java`

_Note_: this code works in Eclipse after adding some environment variables.  Contact Anna if you have trouble running it.

Change hard-coded filenames to yours.  The filterfile ensures that only entities that make it in the hypergraph are parsed here (extraneous reactions, etc., are skipped to ensure consistency across representations).

## Legacy code:
- `StreamSurvey2.java`: non-paramterized version of compound graph connectivity
- `GetPathwayMembers.java`: enumerates pathway members of compound graph.