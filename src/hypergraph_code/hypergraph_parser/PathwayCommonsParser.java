import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.biopax.paxtools.controller.EditorMap;
import org.biopax.paxtools.controller.PropertyEditor;
import org.biopax.paxtools.controller.SimpleEditorMap;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.level3.Process;

/**
 * 
 * This class parses PathwayCommons v10 files
 * Required JAR File: paxtools-4.2.1-no-jena.jar
 * 
 * @author annaritz
 * @version 2019-08-01
 * 
 * Useful Examples:
 * paxtools.pdf
 * PaxtoolsMain.java (see api)
 * example code: https://groups.google.com/forum/#!topic/biopax-discuss/SN4ETE6gSCc 
 * 
 * Note: This is copied/merged from Anna's `biopax-parsers` private git repo.  See `biopax-parsers/java/src/`.
 * It now is part of the pathway-connectivity github repo, with some cosmetic changes (argument parsing, etc.).
 *
 */
public class PathwayCommonsParser {

	/* Static Class Variables */
	public static String filename;
	public static String fileprefix;
	public static boolean verbose = false;
	public static String USAGE = "java PathwayCommonsParser <Reactome.owl> <outprefix> [<pathway-id>] [<verbose>]\n\t<Reactome.owl>: BioPAX Level 3 file (.owl format).\n\t<outprefix>: Output prefix of all signaling pathways.\n\t<pathway-id> (optional): if specified, only runs that pathway id.\n\t<verbose> (optional): the word 'verbose' to print additional information to screen. Optional.\n\nNote: to run from Java source, ReactomeParser requires paxtools-4.2.1-no-jena.jar in your CLASSPATH. If running from the JAR file (using java -jar), it should not be required\n";

	/**
	 * Main method.
	 *
	 * 
	 * @param args - currently empty (filenames are hard-coded)
	 * 
	 * @throws NoSuchFieldException - If an expected field in the BioPAX file is missing.
	 * @throws IOException - Issue with output files.
	 * 
	 */
	public static void main(String[] args) throws NoSuchFieldException, IOException {
		
		if(args.length < 2) {
			System.out.println(USAGE);
			System.err.println("\nThere are "+args.length+" arguments. Exiting. ");
			System.exit(-1);
		} 
		
		filename = args[0];
		fileprefix = args[1];
		
		System.out.println("BIOPAX FILE:" + filename);
		System.out.println("OUTFILE PREFIX:" + fileprefix);
		
		String pathwayID = null;
		if (args.length >= 3 && ! (args[2].equalsIgnoreCase("verbose"))) 
			pathwayID = args[2];
		if (args.length == 4 && ! (args[3].equalsIgnoreCase("verbose")))
			pathwayID = args[3];
		
		if (args.length >= 3 && args[2].equalsIgnoreCase("verbose")) 
			verbose = true;
		if (args.length == 4 && args[3].equalsIgnoreCase("verbose"))
			verbose = true;
		if (args.length == 4 && !verbose) {
			System.out.println(USAGE);
			System.err.println("\nERROR: There are four arguments but one isn't verbose. Exiting.");
			System.exit(-1);
		}
		if (pathwayID != null) 
			System.out.println("Pathway ID " + pathwayID);
		System.out.println("VERBOSE? " + verbose);
		
		PathwayCommonsParser parser = new PathwayCommonsParser();
		parser.readPC(pathwayID);
		System.out.println("Done Parsing PathwayCommons.");
	}
	
	/**
	 * Constructor.
	 */
	public PathwayCommonsParser() { }
	
	/**
	 * Reads the Reactome BioPAX file and parses all pathways
	 * 
	 * NOTE: PAXTools emits many warnings (possibly because I set mergDuplicates
	 * to true).  However, there are no Errors so I assume this is OK.
	 * 
	 * From this thread: https://groups.google.com/forum/#!topic/biopax-discuss/-DwJGfqDHDE
	 * I get many "Redundant attempt to set the inverse link..." warnings (feature - featureOf - protein). 
	 * So what happens here is that a ModificationFeature object is added to the same Protein more than once: using either 'feature' or 'notFeature' biopax property. Such state does not make any sense, and looks like a curation error.
	 *
	 *What you can see there in the Reactome BioPAX L3 file (Caenorhabditis elegans.owl; but I am note sure which version you use, because in my local file Protein3178 does not have this issue, anyway...)
	 *  
	<bp:Protein rdf:ID="Protein394">
      <bp:displayName rdf:datatype="http://www.w3.org/2001/XMLSchema#string">P91918</bp:displayName>
      ...
      <bp:entityReference rdf:resource="#ProteinReference264" />
      <bp:feature rdf:resource="#ModificationFeature48" />
      <bp:feature rdf:resource="#ModificationFeature48" />
      <bp:feature rdf:resource="#FragmentFeature335" />
      ...
	 </bp:Protein>
	 * 
	 * Checked that ModificationFeature100 appears multiple times for Protein761. Above is the reason why this emits a warning.
	 *
	 *
	 *See this code in case there are UnificationXRef Errors: http://www.biopax.org/m2site/paxtools-4.2.0/apidocs/src-html/org/biopax/paxtools/examples/ReactomeEntitySetUnificationXrefFix.html
	 *
	 * 
	 * @throws NoSuchFieldException - If an expected field in the BioPAX file is missing
	 * @throws IOException - Issue with output files
	 * 
	 */
	public void readPC(String pathwayID) throws NoSuchFieldException, IOException {
		
		// Initialize IOHandler. Currently using SimpleIOHandler;
		// JenaIOHandler is more complex (will need to download a different
		// jar).
		SimpleIOHandler io = new SimpleIOHandler(BioPAXLevel.L3);
		io.mergeDuplicates(true);
		
		// Convert the BioPAX file to a PAXTools Model.
		Model model = io.convertFromOWL(new FileInputStream(filename));
		
		// Get implemented visitor specific to Reactome.
		BioPaxVisitor visitor = new BioPaxVisitor(SimpleEditorMap.L3,model,fileprefix,DB.PC,verbose);
		
		Set<Pathway> pathways;		
		if (pathwayID == null) {
			// get list of all pathways
			pathways =  model.getObjects(Pathway.class);
			System.out.println(pathways.size() + " pathways to process.");
			for (Process pathway : pathways) {
				System.out.println("  Pathway " + pathway.getDisplayName());
			}
		} else {
			pathways = new HashSet<Pathway>();
			BioPAXElement pathwayelement = model.getByID(pathwayID);
			if (pathwayelement == null) {
				System.out.println(USAGE);
				System.err.println("\nERROR: Pathway ID '"+pathwayID+"' is not found in the OWL file. Exiting.");
				System.exit(-1);
			} else
				System.out.println("Adding pathway " + ((Pathway)pathwayelement).getDisplayName());
			pathways.add((Pathway)pathwayelement);
		}
		
		// Parse each pathway
		for (Process pathway : pathways) {
			String pathwayfilename = pathway.getDisplayName().replace(' ','-').replace('/','-');
			//if (pathwayfilename.contains("Wnt")) {
			visitor.reset();
			visitor.traverseAndCount((Pathway)pathway,pathwayfilename);
			//}
		}
	}	
}
