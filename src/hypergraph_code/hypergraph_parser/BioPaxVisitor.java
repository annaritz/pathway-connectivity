import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.biopax.paxtools.controller.EditorMap;
import org.biopax.paxtools.controller.PropertyEditor;
import org.biopax.paxtools.controller.SimpleEditorMap;
import org.biopax.paxtools.controller.Traverser;
import org.biopax.paxtools.controller.Visitor;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;

/**
 * This visitor class is used to traverse BioPAX files.
 * See the Visitor class in PAXTools.
 * 
 * @author annaritz
 * @version 2013-11-19
 *
 */
public class BioPaxVisitor implements Visitor {

	/* Class Variables */
	
	// delimiter to separate elemnets in a single column
	private static final String DELIM = ";";

	// enum describing the database
	private DB database;
	
	// boolean whether to be verbose or not
	private boolean verbose;
	
	// step of the recursion (may be output if verbose = true)
	public static int RECURSIONSTEP = 0;

	// RECTOME_DB is the current version of the Reactome database.
	//private static final String REACTOME_RELEASE = "Reactome Database ID Release 51";

	// ALT_IDS is the set of other IDs to include in the elements.txt file.  
	// Every element SHOULD have at least one ALT_ID
	// In Reactome, some elements are labeled PhysicalEntity: these only have ReactomeIDs.
	private static final String[] ALT_IDS = {"Reactome","UniProt","UniProt Isoform","ChEBI","ENSEMBL","GO","Pubmed","hgnc symbol","uniprot knowledgebase"};

	// origModel is the original BioPAX object that is traversed.
	private Model origModel;

	// Objects for the Visitor interface. See biopax.pdf for an example.
	private Traverser traverser;
	private EditorMap editorMap;

	// Keep track of the number of properties and objects we record.
	public HashSet<String> visitedProperties;
	private HashSet<String> visitedObjects;

	// Keep track of the number of objects we record and skip.
	private HashMap<Class<?>,Integer> recordedCounts;
	private HashMap<Class<?>,Integer> skippedCounts;

	// Variables for writing output
	private String outputprefix;
	private BufferedWriter elementwriter,reactionwriter,controlwriter,subpathwaywriter, complexwriter, entitysetwriter;

	// These are all the properties I've seen so far:
	// the properties that are NOT commented out are ignored
	// in the visit() method (we don't need to traverse).
	Set<String> propertiesToIgnore = new HashSet<String> (Arrays.asList(
			"entityReference",
			//"stepProcess",
			"dataSource",
			"modificationType",
			//"pathwayComponent",
			//"controlled",
			"sequenceIntervalEnd",
			//"memberPhysicalEntity", // NOT COMMENTED OUT FOR NCI-PID!! Modified in constructor.
			"organism",
			"componentStoichiometry",
			"relationshipType",
			//"physicalEntity",
			//"component",
			"featureLocation",
			"xref",
			"participantStoichiometry",
			//"left",
			"feature",
			//"controller",
			//"pathwayOrder",
			"cellularLocation",
			"nextStep", // COMMENTED OUT FOR NCI-PID!! Modified in constructor.
			//"right",
			"sequenceIntervalBegin"));
	

	/* ***************************** */
	/* Constructor and Reset Methods */
	/* ***************************** */

	/**
	 * Constructor.
	 * @param em - EditorMap
	 * @param om - Original Model
	 * @param prefix - output file prefix
	 */
	public BioPaxVisitor(EditorMap em, Model om, String prefix, DB db,boolean v) {
		editorMap = em;
		origModel = om;
		traverser = new Traverser(editorMap,this);
		outputprefix = prefix;
		database = db;
		verbose = v;

		// Initialize variables for tracking objects and properties.
		recordedCounts = new HashMap<Class<?>,Integer>();
		visitedProperties = new HashSet<String>();
		visitedObjects = new HashSet<String>();
		skippedCounts = new HashMap<Class<?>,Integer>();
		
		// if NCI-PID, modify the propertiesToIgnore set.
		if (db == DB.NCIPID) {
			// NCI-PID uses memberPhysicalEntity
			// NCI-PID does NOT iterate through subpathways with the nextStep tag.
			// Make these changes to the propertiesToIgnore set.
			propertiesToIgnore.add("memberPhysicalEntity");
			propertiesToIgnore.remove("nextStep");
		}
	}

	/**
	 * Resets the variables for tracking objects and properties.
	 */
	public void reset() {
		recordedCounts.clear();
		visitedProperties.clear();
		visitedObjects.clear();
		skippedCounts.clear();
	}

	/* *************************** */
	/* Traversal and Visit Methods */
	/* *************************** */

	/**
	 * @Override
	 *
	 * This method is called on every element that is traversed.  
	 * 
	 * @param domain - Parent object
	 * @param element - Currently-traversed object (visited object)
	 * @param model - Model to traverse
	 * @param editor - PropertyEditor
	 * 
	 */
	public void visit(BioPAXElement domain, Object element, Model model,
			PropertyEditor<?, ?> editor) {
		// We are only interested in the BioPAXElements since
		// primitive fields are always copied by value (see paxtools.pdf)
		if (element == null || !(element instanceof BioPAXElement)) 
			return;

		// Convert element to type BioPAXElement.
		BioPAXElement bpelement = (BioPAXElement) element;

		// If we have seen this element, skip.
		if (visitedObjects.contains(bpelement.getRDFId())) 
			return;
		

		// If this property is one we wish to ignore, return.
		if (propertiesToIgnore.contains(editor.getProperty())) 
			return;
		
		if (verbose) {
			System.out.println("STARTING WITH "+ bpelement.getRDFId());
			RECURSIONSTEP++;
		}
		
		// Add this property to the list of visited properties
		visitedProperties.add(editor.getProperty());	

		// Add this element to the list of visited objects.
		visitedObjects.add(bpelement.getRDFId());

		try {
			// write element to file
			writeBioPAXElement(bpelement);

			// increment the counts for recording this class.
			if (!recordedCounts.containsKey(bpelement.getClass())) 
				recordedCounts.put(bpelement.getClass(),0);
			recordedCounts.put(bpelement.getClass(),recordedCounts.get(bpelement.getClass())+1);

		} catch (NoSuchFieldException | IOException e) {
			System.out.println(e.getMessage());
			System.err.println(e.getMessage());
		}

		// Recurse. Note this is ONLY if we haven't seen the element yet.
		traverser.traverse(bpelement, model);

		if (verbose) {
			System.out.println("DONE WITH " + bpelement.getRDFId());
			RECURSIONSTEP--;
		}
	}

	/**
	 * Traverses a pathway BioPAX object, calling the visit() method 
	 * for each traversed element.
	 * 
	 * @param pathway - BioPAX Pathway to traverse
	 * @throws IOException - if output files cannot be written.
	 * 
	 */
	public void traverseAndCount(Pathway pathway,String pathwayname) throws IOException {
		// New Dec 2016: replace any DELIMs in pathwayname with semicolon
		if (DELIM.equals(";")) // if delim is semicolon, replace all with colons.
			pathwayname = pathwayname.replace(DELIM,":");
		else // in case delimiter is something other than semicolom
			pathwayname = pathwayname.replace(DELIM,";");
		if (pathwayname.contains(DELIM)) {
			System.out.println("Error: "+pathwayname);
			System.exit(-1);
		}

		// Reset counts (allows us to call traverse() on different pathways)
		reset();

		System.out.println("\n--------------------------------------------------------------------------------");
		if (pathway.getDisplayName() != null)
			System.out.println(" Traversing pathway '"+pathway.getDisplayName()+"'");
		else if (pathway.getStandardName() != null)
			System.out.println(" Traversing pathway '"+pathway.getStandardName()+"'");
		else if (pathway.getName().size() > 0)
			System.out.println(" Traversing pathway '"+StringUtils.join(pathway.getName(),DELIM)+"'");
		else
			System.out.println("Traversing Unnamed Pathway");

		// Initialize filewriters.
		initializeFileWriters(pathwayname);

		// Call the traverse() function on the Traverser object with, 
		// traversing the original model specified in the Constructor.
		traverser.traverse(pathway, origModel);

		// Close filewriters.
		closeFileWriters();

		// Print information to console.
		printProperties();
		printRecordedElements();
		printSkippedElements();
		System.out.println("Pathway written to files with prefix: " + pathwayname);
		System.out.println("---------------------------------------------------------------------------------\n");
	}	
	
	/**
	 * Traverses a pathway BioPAX object, calling the visit() method 
	 * for each traversed element.
	 * 
	 * @param pathway - BioPAX Pathway to traverse
	 * @throws IOException - if output files cannot be written.
	 * @throws NoSuchFieldException 
	 * 
	 */
	public void traverseAndCount(Model model, String name) throws IOException, NoSuchFieldException {

		// Reset counts (allows us to call traverse() on different pathways)
		reset();

		System.out.println("\n---------------------------------------------------------------------------------");
		System.out.println("Traversing pathway '" + name+"'");

		// Initialize filewriters.
		String pathwayfilename = name.replace(' ','-').replace('/','-');
		initializeFileWriters(pathwayfilename);
	
		// Call the traverse() function on the Traverser object with, 
		// traversing the original model specified in the Constructor.
		Set<BioPAXElement> elements = model.getObjects();
		for (BioPAXElement e : elements) {

			// First, write this element! Traverser traverses children only.
			if (!visitedObjects.contains(e.getRDFId())) {
				writeBioPAXElement(e);
				visitedObjects.add(e.getRDFId());
				// NOTE: properties are NOT checked and added here.
			}
			traverser.traverse(e,origModel);
		}

		// Close filewriters.
		closeFileWriters();

		// Print information to console.
		printProperties();
		printRecordedElements();
		printSkippedElements();
		System.out.println("Pathway written to following files with prefix " + pathwayfilename);
		System.out.println("---------------------------------------------------------------------------------\n");
	}
	
	/**
	 * initialize all file writers
	 * @param pathwayfilename
	 * @throws IOException
	 */
	public void initializeFileWriters(String pathwayfilename) throws IOException {
		String reactionfilename = outputprefix+File.separator+pathwayfilename+"-reactions.txt";	
		String elementfilename = outputprefix+File.separator+pathwayfilename+"-elements.txt";
		String complexfilename = outputprefix+File.separator+pathwayfilename+"-complexes.txt";
		String controlfilename = outputprefix+File.separator+pathwayfilename+"-controls.txt";
		String subpathwayfilename = outputprefix+File.separator+pathwayfilename+"-subpathways.txt";
		reactionwriter = new BufferedWriter(new FileWriter(reactionfilename));
		reactionwriter.write("#id\tE\tF\tR\tinteractionType\tE&Fstoichiometry\tspontaneous\tevidence\txref\n");
		elementwriter = new BufferedWriter(new FileWriter(elementfilename));
		elementwriter.write("#id\tdisplayname\taltnames\taltIDs\telementtype\tfeatures\tcellularLocation\n");
		complexwriter = new BufferedWriter(new FileWriter(complexfilename));
		complexwriter.write("#id\tname\tfeatures\tcelularLocation\telements\tstoichiometry\n");
		controlwriter = new BufferedWriter(new FileWriter(controlfilename));
		controlwriter.write("#id\tcontrollingElements\tcontrolledElements\tcontrolType\tcatalysisdir\tevidence\txref\n");
		subpathwaywriter = new BufferedWriter(new FileWriter(subpathwayfilename));
		subpathwaywriter.write("#id\tname\txref\n");
		
		// Only write entity set file for Reactome or NCI-PID.
		if (database == DB.REACTOME || database == DB.NCIPID || database==DB.PC) {
			String entitysetfilename = outputprefix+File.separator+pathwayfilename+"-entitysets.txt";
			entitysetwriter = new BufferedWriter(new FileWriter(entitysetfilename));
			entitysetwriter.write("#id\tname\tentitySetIDs\tfeatures\tcellularLocation\n");
		}
	}
	
	/**
	 * close all file writers
	 * @throws IOException
	 */
	public void closeFileWriters() throws IOException {
		reactionwriter.close();
		elementwriter.close();
		controlwriter.close();
		complexwriter.close();
		subpathwaywriter.close();
		if (database == DB.REACTOME || database == DB.NCIPID || database == DB.PC) 
			entitysetwriter.close();
	}

	/* *************************** */
	/* Methods to Print Statistics */
	/* *************************** */

	/**
	 * Prints the list of properties recorded.
	 * This is useful to compare against the propertiesToIgnore variable.
	 */
	public void printProperties() {
		System.out.println("\n"+visitedProperties.size()+" properties encountered.");
		for (String property : visitedProperties) 
			System.out.println(" "+property);
		System.out.println();
	}

	/**
	 * Prints the number of each class that were recorded.
	 */
	public void printRecordedElements() {
		System.out.println("\n"+recordedCounts.size()+" object types traversed.");
		System.out.println("Class\t# Traversed");
		for (Map.Entry<Class<?>, Integer> entry : recordedCounts.entrySet()) 
			System.out.println(entry.getKey()+"\t"+entry.getValue());
		System.out.println();
	}

	/**
	 * Prints the number of each class that were skipped.
	 */
	public void printSkippedElements() {
		System.out.println("\n"+skippedCounts.size()+" object types were not written to file.");
		for (Map.Entry<Class<?>, Integer> entry : skippedCounts.entrySet())
			System.out.println(entry.getKey()+"\t"+entry.getValue());
		System.out.println();
	}

	/* *********************** */
	/* Methods to Write Output */
	/* *********************** */

	/**
	 * Determines the type of BioPAXElement and calls the corresponding
	 * write() function.
	 * 
	 * @param bpelement - BioPAX element to write to file
	 * @throws NoSuchFieldException - an expected field is missing
	 * @throws IOException - an issue with writing to output files
	 * 
	 */
	private void writeBioPAXElement(BioPAXElement bpelement) throws NoSuchFieldException, IOException {

		// Check Entity Sets right away.  Complexes, Proteins, etc. may be Entity Sets.
		boolean isES = false;
		
		// In Reactome,there is a comment that says "Converted from EntitySet in Reactome".
		// EntitySets may be proteins, small molecules, complexes, etc.
		// http://wiki.reactome.org/index.php/EntitySet
		// http://wiki.reactome.org/index.php/Glossary_Data_Model#EntitySet_.5C.5Bsuperclass.5C.5D
		// If this is the case, then add it to -entitysets.txt file. Go through comments to check.
		if (database == DB.REACTOME && bpelement instanceof PhysicalEntity) {	
			for(String comment :  ((PhysicalEntity)bpelement).getComment()) 
				if (comment.contains("Converted from EntitySet in Reactome")) {
					isES = true;
					break;
				}			
		}
		
		// In NCIPID, The EntityReference will contain OTHER EntityReferences.
		// According to the bioPAX API, Complexes cannot contain other entity references.
		// Thus, EntitySets may be proteins/small molecules/Rna/Dna/etc.
		// E.g., DVL is a ProteinReference, and DVL1/DVL2/DVL3 are EntityReferences
		// within the ProteinReference.
		// To check this, get the entity reference and then see if there are any
		// members of the entity reference.
		if ((database == DB.NCIPID || database == DB.PC) && bpelement instanceof PhysicalEntity) {
			EntityReference e = null;
			if (bpelement instanceof Protein) {
				e = ((Protein)bpelement).getEntityReference();
				if (e.getMemberEntityReference().size() > 0)  
					isES = true;
			} else if (bpelement instanceof SmallMolecule) {
				e = ((SmallMolecule)bpelement).getEntityReference();
				if (e.getMemberEntityReference().size() > 0)  
					isES = true;
			} else if (bpelement instanceof Dna) {
				e = ((Dna)bpelement).getEntityReference();
				if (e.getMemberEntityReference().size() > 0)  
					isES = true;
			} else if (bpelement instanceof Rna) {
				e = ((Rna)bpelement).getEntityReference();
				if (e.getMemberEntityReference().size() > 0)  
					isES = true;
			} else if (bpelement instanceof DnaRegion) {
				e = ((DnaRegion)bpelement).getEntityReference();
				if (e.getMemberEntityReference().size() > 0)  
					isES = true;
			} else if (bpelement instanceof Rna) {
				e = ((RnaRegion)bpelement).getEntityReference();
				if (e.getMemberEntityReference().size() > 0)  
					isES = true;
			}		
		}
		
		
		
		if (isES) {
			writeEntitySet((PhysicalEntity)bpelement);
			return;
		}
		
		if (bpelement instanceof Pathway) 
			writeSubPathway((Pathway)bpelement);
		else if (bpelement instanceof Complex) 
			writeComplex((Complex)bpelement);
		else if (bpelement instanceof BiochemicalReaction)
			writeBiochemicalReaction((BiochemicalReaction)bpelement);
		else if (bpelement instanceof Catalysis) 
			writeCatalysis((Catalysis)bpelement);
		else if (bpelement instanceof Control) 
			writeControl((Control)bpelement);
		else if (bpelement instanceof MolecularInteraction)
			writeMolecularInteraction((MolecularInteraction)bpelement);
		else if (bpelement instanceof Transport)
			writeTransport((Transport)bpelement);
		else if (bpelement instanceof ComplexAssembly)
			writeComplexAssembly((ComplexAssembly)bpelement);
		else if (bpelement instanceof Conversion)
			writeConversion((Conversion)bpelement);
		else if (bpelement instanceof PhysicalEntity) 
			writeMappableElement((PhysicalEntity)bpelement);
		else if (bpelement instanceof TemplateReaction)
			writeTemplateReaction((TemplateReaction)bpelement);
		else // Increment skippedCounts counter for this class.
			if(!skippedCounts.containsKey(bpelement.getClass()))
				skippedCounts.put(bpelement.getClass(),1);
			else
				skippedCounts.put(bpelement.getClass(),skippedCounts.get(bpelement.getClass())+1);
	}

	/**
	 * Writes PhysicalEntity elements that have IDs in alternate namespaces.
	 * Writes to elements.txt file.
	 * 
	 * @param bpelement - PhysicalEntity element to write.
	 * @throws IOException - issue writing to elements.txt output file
	 */
	private void writeMappableElement(PhysicalEntity bpelement) throws IOException {
		if (verbose)
			System.out.println("  "+RECURSIONSTEP+": "+getID(bpelement)+" --> " + bpelement.getDisplayName());
		// Get all possible IDs for this element
		HashMap<String,String> ids = getAllIDs(bpelement);

		// Write original ID
		elementwriter.write(getID(bpelement));
				
		// Write display name
		// Reactome sometimes has an empty display name; in that case write Name id.
		if (bpelement.getDisplayName() != null) 
			elementwriter.write("\t"+bpelement.getDisplayName().replace(' ','-'));
		else 
			elementwriter.write("\t"+ids.get("Name"));
		
		// Write Alt names
		elementwriter.write("\t"+getName(bpelement,ids));
		
		//Report an error if reactome entry doesn't have reactome ID
		//if(database == DB.REACTOME && !ids.containsKey(REACTOME_RELEASE)) {
		//	System.out.println("ERROR: " + bpelement + " does not have a current Reactome Database ID.");
		//	System.exit(-1);
		//}

		// Write alternate IDs (usually only one is not 'None')
		elementwriter.write("\t"+getAltIDs(ids));
		
		// write element type
		if (bpelement instanceof Protein)
			elementwriter.write("\tprotein");
		else if (bpelement instanceof SmallMolecule)
			elementwriter.write("\tsmall-molecule");
		else if (bpelement instanceof Rna)
			elementwriter.write("\trna");
		else if (bpelement instanceof Dna)
			elementwriter.write("\tdna");
		else if (database==DB.CELLDESIGNER && bpelement instanceof DnaRegion) {
			// in Pavel's model (and perhaps Kitano maps as well), Genes are 
			// reported as DnaRegion elements.  Treat this as dna.
			elementwriter.write("\tdna");
		} else {
			if (database==DB.REACTOME)
				System.out.println("\tNote: bp element " + getID(bpelement) + "("+bpelement.getDisplayName()+") is not a known element type");
			else
				System.out.println("WARNING: bp element " + getID(bpelement) + " is not a known element type");
			elementwriter.write("\tphysical-entity");
		}
		
		// Write features
		elementwriter.write("\t"+getFeatures(bpelement.getFeature()));
		
		// Write cellular location
		if(bpelement.getCellularLocation() == null)
			elementwriter.write("\tNone\n");
		else
			elementwriter.write("\t"+bpelement.getCellularLocation().toString().replace(' ','-')+"\n");
		
	}

	/**
	 * Writes a complex to complexes.txt
	 * 
	 * @param bpelement - Complex to write
	 * @throws IOException - issue with complexes.txt file
	 * 
	 */
	private void writeComplex(Complex bpelement) throws IOException {
		if (verbose)
			System.out.println("  "+RECURSIONSTEP+": "+getID(bpelement)+" --> " + bpelement.getDisplayName());
		// Write original ID
		complexwriter.write(getID(bpelement));
		
		// Write Name
		complexwriter.write("\t"+getName(bpelement,getAllIDs(bpelement)));
				
		// Write features
		complexwriter.write("\t"+getFeatures(bpelement.getFeature()));
		
		// Write cellular location
		if(bpelement.getCellularLocation() == null)
			complexwriter.write("\tNone");
		else
			complexwriter.write("\t"+bpelement.getCellularLocation().toString().replace(' ','-'));

		// Write complex elements
		Set<PhysicalEntity> components = bpelement.getComponent();
		complexwriter.write("\t"+convertList(components));
		
		// Write stoichiometry
		complexwriter.write("\t"+getStoichiometry(bpelement.getComponentStoichiometry())+"\n");

	}

	/**
	 * Writes a biochemical reaction as (id,E,F,R,interactiontype,description)
	 * Writes to reactions.txt
	 * 
	 * @param bpelement - BiochmicalReaction to write
	 * @throws IOException - issue with reactions.txt file
	 */
	private void writeBiochemicalReaction(BiochemicalReaction bpelement) throws IOException {
		if (verbose)
			System.out.println("  "+RECURSIONSTEP+": "+getID(bpelement)+" --> " + bpelement.getDisplayName());
		
		// Write ID
		reactionwriter.write(getID(bpelement));

		// Write E
		reactionwriter.write("\t"+convertList(bpelement.getLeft()));

		// Write F
		reactionwriter.write("\t"+convertList(bpelement.getRight()));

		// Write R
		Set<Control> controllers  = bpelement.getControlledOf();
		if(controllers.size() > 0)
			reactionwriter.write("\t"+convertList(controllers));
		else
			reactionwriter.write("\tNone");

		// Write interactiontype
		if(bpelement.getInteractionType().size() != 0) 
			reactionwriter.write("\t"+convertList(bpelement.getInteractionType()));
		else
			reactionwriter.write("\tBiochemicalReaction");

		// Write stoichiometry
		reactionwriter.write("\t"+getStoichiometry(bpelement.getParticipantStoichiometry()));
		
		// Write spontaneous
		if(bpelement.getSpontaneous() == null)
			reactionwriter.write("\tnull");
		else if(bpelement.getSpontaneous() == true)
			reactionwriter.write("\ttrue");
		else
			reactionwriter.write("\tfalse");
		
		// Write evidence
		reactionwriter.write("\t"+getEvidence(bpelement.getEvidence()));
				
		// Write xrefs:
		reactionwriter.write("\t"+getXrefs(bpelement.getXref())+"\n");
	}

	/**
	 * Writes a catalysis object as (id,controller,controlled,controltype,description)
	 * Writes to controls.txt
	 * 
	 * @param bpelement - Catalysis object to write
	 * @throws IOException - issue with controls.txt
	 * 
	 */
	private void writeCatalysis(Catalysis bpelement) throws IOException {
		if (verbose)
			System.out.println("  "+RECURSIONSTEP+": "+getID(bpelement)+" --> " + bpelement.getDisplayName());
		
		// Write id
		controlwriter.write(getID(bpelement));

		// Write controllers
		if(bpelement.getControlled() != null)
			controlwriter.write("\t"+convertList(bpelement.getController()));
		else
			controlwriter.write("\tNone");

		// Write controlled
		if(bpelement.getControlled() != null)
			controlwriter.write("\t"+convertList(bpelement.getControlled()));
		else
			controlwriter.write("\tNone");

		// Write controltype
		if(bpelement.getControlType() != null) 
			controlwriter.write("\t"+bpelement.getControlType().name());
		else
			controlwriter.write("\tNone");

		// Write catalysisDirection
		if(bpelement.getCatalysisDirection() != null) 
			controlwriter.write("\t"+bpelement.getCatalysisDirection().name());
		else
			controlwriter.write("\tNone");
		
		// Write evidence
		controlwriter.write("\t"+getEvidence(bpelement.getEvidence()));

		// Write xrefs:
		controlwriter.write("\t"+getXrefs(bpelement.getXref())+"\n");
	}

	/**
	 * Writes a Control object as (id,controller,controlled,controltype,description)
	 * Writes to controls.txt
	 * 
	 * @param bpelement - Control object to write
	 * @throws IOException - issue with controls.txt
	 * 
	 */
	private void writeControl(Control bpelement) throws IOException {
		if (verbose)
			System.out.println("  "+RECURSIONSTEP+": "+getID(bpelement)+" --> " + bpelement.getDisplayName());
		
		// Write id
		controlwriter.write(getID(bpelement));

		// Write controllers
		if(bpelement.getControlled() != null)
			controlwriter.write("\t"+convertList(bpelement.getController()));
		else
			controlwriter.write("\tNone");

		// Write controlled
		if(bpelement.getControlled() != null)
			controlwriter.write("\t"+convertList(bpelement.getControlled()));
		else
			controlwriter.write("\tNone");

		// Write controltype
		if(bpelement.getControlType() != null) 
			controlwriter.write("\t"+bpelement.getControlType().name());
		else if (database == DB.CELLDESIGNER) {
			// if no control type is specified, set to ACTIVATION
			// for CellDesigner
			controlwriter.write("\tACTIVATION");
		} else
			controlwriter.write("\tNone");

		// Write catalysisDirection - None if it's a Control Element
		controlwriter.write("\tNone");
		
		// Write evidence
		controlwriter.write("\t"+getEvidence(bpelement.getEvidence()));

		// Write xrefs:
		controlwriter.write("\t"+getXrefs(bpelement.getXref())+"\n");
	}

	/**
	 * Writes a MolecularInteraction object as (id,participants,None,None,interactiontype,description)
	 * Writes to reactions.txt
	 * As far as I can tell, these only appear in NetPath pathways as undirected interactions; 
	 * Thus, the participants are put in E (undirected hyperedge).
	 * 
	 * @param bpelement - MolecularInteractiom object to write
	 * @throws IOException - issue with reactions.txt
	 * 
	 */
	private void writeMolecularInteraction(MolecularInteraction bpelement) throws IOException {
		if (verbose)
			System.out.println("  "+RECURSIONSTEP+": "+getID(bpelement)+" --> " + bpelement.getDisplayName());
		
		if (database != DB.NETPATH) {
			System.err.println("ERROR! Molecular Interactions should only appear in NetPath Pathways.");
			System.exit(-1);
		}
		
		// Write id
		reactionwriter.write(getID(bpelement));

		// Write E
		Set<Entity> components = bpelement.getParticipant();
		reactionwriter.write("\t"+convertList(components));

		// Write F,R = None,None
		reactionwriter.write("\tNone\tNone");

		// Write interactiontype
		if(bpelement.getInteractionType().size() != 0) 
			reactionwriter.write("\t"+convertList(bpelement.getInteractionType()));
		else
			reactionwriter.write("\tMolecularInteraction");

		// Write stoichiometry (MolecularInteraction has no label for stoichiometry)
		reactionwriter.write("\tNone");
		
		// Write spontaneous (MolecularInteraction has no label of spontaneous)
		reactionwriter.write("\tnull");
		
		// Write evidence
		reactionwriter.write("\t"+getEvidence(bpelement.getEvidence()));
		
		// Write xrefs:
		reactionwriter.write("\t"+getXrefs(bpelement.getXref())+"\n");
	}

	/**
	 * Writes a Transport object as (id,left,right,None,None,description)
	 * Writes to reactions.txt
	 * 
	 * @param bpelement - Transport object to write
	 * @throws IOException - issue with reactions.txt
	 * @throws NoSuchFieldException 
	 * 
	 */
	private void writeTransport(Transport bpelement) throws IOException, NoSuchFieldException {
		if (verbose)
			System.out.println("  "+RECURSIONSTEP+": "+getID(bpelement)+" --> " + bpelement.getDisplayName());
		
		// Write id
		reactionwriter.write(getID(bpelement));

		// Write left
		Set<PhysicalEntity> components = bpelement.getLeft();
		reactionwriter.write("\t"+convertList(components));

		// Write right
		components = bpelement.getRight();
		reactionwriter.write("\t"+convertList(components));

		// Write R = None
		reactionwriter.write("\tNone");
		
		// Write interactiontype
		if(bpelement.getInteractionType().size() != 0) 
			reactionwriter.write("\t"+convertList(bpelement.getInteractionType()));
		else
			reactionwriter.write("\tTransport");

		// Write stoichiometry
		reactionwriter.write("\t"+getStoichiometry(bpelement.getParticipantStoichiometry()));

		// Write spontaneous
		if(bpelement.getSpontaneous() == null)
			reactionwriter.write("\tnull");
		else if(bpelement.getSpontaneous() == true)
			reactionwriter.write("\ttrue");
		else
			reactionwriter.write("\tfalse");
	
		// Write evidence
		reactionwriter.write("\t"+getEvidence(bpelement.getEvidence()));
				
		// Write xrefs:
		reactionwriter.write("\t"+getXrefs(bpelement.getXref())+"\n");
	}

	/**
	 * Writes a complex assembly as (id,E,F,R,interactiontype,description)
	 * Writes to reactions.txt
	 * 
	 * @param bpelement - ComplexAssembly to write
	 * @throws IOException - issue with reactions.txt file
	 */
	private void writeComplexAssembly(ComplexAssembly bpelement) throws IOException {
		if (verbose)
			System.out.println("  "+RECURSIONSTEP+": "+getID(bpelement)+" --> " + bpelement.getDisplayName());
		
		// Write ID
		reactionwriter.write(getID(bpelement));

		// Write E
		reactionwriter.write("\t"+convertList(bpelement.getLeft()));

		// Write F
		reactionwriter.write("\t"+convertList(bpelement.getRight()));

		// Write R
		Set<Control> controllers  = bpelement.getControlledOf();
		if(controllers.size() > 0)
			reactionwriter.write("\t"+convertList(controllers));
		else
			reactionwriter.write("\tNone");

		// Write interactiontype
		if(bpelement.getInteractionType().size() != 0) 
			reactionwriter.write("\t"+convertList(bpelement.getInteractionType()));
		else
			reactionwriter.write("\tComplexAssembly");
		
		// Write stoichiometry
		reactionwriter.write("\t"+getStoichiometry(bpelement.getParticipantStoichiometry()));

		// Write spontaneous
		if(bpelement.getSpontaneous() == null)
			reactionwriter.write("\tnull");
		else if(bpelement.getSpontaneous() == true)
			reactionwriter.write("\ttrue");
		else
			reactionwriter.write("\tfalse");
		
		// Write evidence
		reactionwriter.write("\t"+getEvidence(bpelement.getEvidence()));
				
		// Write xrefs:
		reactionwriter.write("\t"+getXrefs(bpelement.getXref())+"\n");
	}

	/**
	 * Writes conversion as (id,E,F,R,interactiontype,description)
	 * Writes to reactions.txt
	 * 
	 * @param bpelement - Conversion to write
	 * @throws IOException - issue with reactions.txt file
	 */
	private void writeConversion(Conversion bpelement) throws IOException {
		if (verbose)
			System.out.println("  "+RECURSIONSTEP+": "+getID(bpelement)+" --> " + bpelement.getDisplayName());
		
		// Write ID
		reactionwriter.write(getID(bpelement));

		// Write E
		reactionwriter.write("\t"+convertList(bpelement.getLeft()));

		// Write F
		reactionwriter.write("\t"+convertList(bpelement.getRight()));

		// Write R
		Set<Control> controllers  = bpelement.getControlledOf();
		if(controllers.size() > 0)
			reactionwriter.write("\t"+convertList(controllers));
		else
			reactionwriter.write("\tNone");

		// Write interactiontype
		if(bpelement.getInteractionType().size() != 0) 
			reactionwriter.write("\t"+convertList(bpelement.getInteractionType()));
		else
			reactionwriter.write("\tConversion");

		// Write stoichiometry
		reactionwriter.write("\t"+getStoichiometry(bpelement.getParticipantStoichiometry()));

		// Write spontaneous
		if(bpelement.getSpontaneous() == null)
			reactionwriter.write("\tnull");
		else if(bpelement.getSpontaneous() == true)
			reactionwriter.write("\ttrue");
		else
			reactionwriter.write("\tfalse");

		// Write evidence
		reactionwriter.write("\t"+getEvidence(bpelement.getEvidence()));
		
		// Write xrefs:
		reactionwriter.write("\t"+getXrefs(bpelement.getXref())+"\n");
	}

	/**
	 * Writes a template reaction as (id-templatereaction,id,products,None,None,description)
	 * Writes to reactions.txt
	 * 
	 * @param bpelement - ComplexAssembly to write
	 * @throws IOException - issue with reactions.txt file
	 */
	private void writeTemplateReaction(TemplateReaction bpelement) throws IOException {
		if (verbose)
			System.out.println("  "+RECURSIONSTEP+": "+getID(bpelement)+" --> " + bpelement.getDisplayName());
		
		// Write ID (in this case, ID = E)
		reactionwriter.write(getID(bpelement));

		// Write E
		reactionwriter.write("\t"+getID(bpelement));

		// Write F
		reactionwriter.write("\t"+convertList(bpelement.getProduct()));

		// Write R
		reactionwriter.write("\tNone");

		// Write interactiontype
		reactionwriter.write("\tTemplateReaction(Transcription)");

		// Write stoichiometry (none for templatereaction)
		reactionwriter.write("\tNone");

		// Write spontaneous (templatereaction is null)
		reactionwriter.write("\tnull");

		// Write evidence
		reactionwriter.write("\t"+getEvidence(bpelement.getEvidence()));
				
		// Write xrefs:
		reactionwriter.write("\t"+getXrefs(bpelement.getXref())+"\n");
	}

	/**
	 * Writes pathway ID and name to a pathway file.
	 * This is necessary when we have pathway "stubs" in BioPAX files (see NCI-PID).
	 * @throws IOException if there's an issue writing to the subpathway file.
	 */
	private void writeSubPathway(Pathway bpelement) throws IOException {
		if (verbose)
			System.out.println("  "+RECURSIONSTEP+": "+getID(bpelement)+" --> " + bpelement.getDisplayName());
		
		String id = getID(bpelement);
		String desc;
		if (bpelement.getDisplayName() != null)
			desc = bpelement.getDisplayName();
		else if (bpelement.getStandardName() != null)
			desc = bpelement.getStandardName();
		else if (bpelement.getName().size() > 0)
			desc = StringUtils.join(bpelement.getName(),DELIM);
		else
			desc = "None";

		//System.out.println(" writing sub-pathway '"+desc+"'");
		subpathwaywriter.write(id+"\t"+desc.replace(' ','-'));
		
		// Write xrefs:
		subpathwaywriter.write("\t"+getXrefs(bpelement.getXref())+"\n");
	}
	
	/**
	 * Writes PhysicalEntity elements that have IDs in alternate namespaces.
	 * Writes to entitysets.txt file.
	 * 
	 * 
	 * @param bpelement - PhysicalEntity element to write.
	 * @throws IOException - issue writing to elements.txt output file
	 */
	private void writeEntitySet(PhysicalEntity bpelement) throws IOException {
		if (verbose)
			System.out.println("  "+RECURSIONSTEP+": "+getID(bpelement)+" --> " + bpelement.getDisplayName());
		
		// ONLY FOR REACTOME and NCI-PID!
		if(database != DB.REACTOME && database != DB.NCIPID && database != DB.PC)  {
			System.out.println("ERROR: EntitySets only belong to Reactome Pathways. Exiting.");
			System.exit(-1);
		}
		
		// Get all possible IDs for this element
		HashMap<String,String> ids = getAllIDs(bpelement);


		//Report an error if reactome entry doesn't have reactome ID
		//if(database == DB.REACTOME && !ids.containsKey(REACTOME_RELEASE)) {
		//	System.out.println("ERROR: " + bpelement + " does not have a current Reactome Database ID.");
		//	System.exit(-1);
		//}
		// Write original ID
		entitysetwriter.write(getID(bpelement));
				
		// Write Display Name
		if (database == DB.REACTOME)
			entitysetwriter.write("\t"+bpelement.getDisplayName().replace(' ','-'));
		else
			entitysetwriter.write("\t"+getName(bpelement,ids));
		
		// write EntitySet Members
		entitysetwriter.write("\t"+convertList(bpelement.getMemberPhysicalEntity()));
		
		// Write features
		entitysetwriter.write("\t"+getFeatures(bpelement.getFeature()));
			
		// Write cellular location
		if(bpelement.getCellularLocation() == null)
			entitysetwriter.write("\tNone\n");
		else
			entitysetwriter.write("\t"+bpelement.getCellularLocation().toString().replace(' ','-')+"\n");
		
	}
	
	/* ************************** */
	/* Concatenation Methods      */
	/* ************************** */

	public String getName(PhysicalEntity bpelement,HashMap<String,String> ids) {
		String namestr= "";

		if (ids.containsKey("Name")) {
			// first check if mapped ids has a "Name" entry (small molecules and some proteins).
			// This is a "cleaner" name than below.
			namestr = ids.get("Name");
		} else if (bpelement.getName().size() > 0) {
			// then check if element has a "name" tag (proteins).
			// replace spaces with dashes. Sometimes alt ids show up here (Uniprot, ChEBI, etc.)
			ArrayList<String> namelist = new ArrayList<String>();
			for(String s : bpelement.getName()) 
				namelist.add(s.replace(' ','-'));
			Collections.sort(namelist,String.CASE_INSENSITIVE_ORDER);
			namestr = StringUtils.join(namelist,DELIM);
		} else if (bpelement.getDisplayName() != null) {
			// then check if the element has a display name tag (complexes).
			// replace spaces with dashes.
			namestr = bpelement.getDisplayName().replace(' ','-');
		} else {
			// No name associated with this element.
			System.out.println("ERROR: " + bpelement + " does not map to any common name.");
			namestr = "None";
			System.out.println(getID(bpelement));
			System.out.println(namestr);
			System.out.println(bpelement.getName());
			System.out.println(bpelement.getName().size());
			listProperties(bpelement);
			System.exit(-1);
		}		
		
		
		
		return namestr;
	}
	
	/**
	 * Return a concatenated string of AltID references.
	 * @param ids
	 * @return
	 */
	public String getAltIDs(HashMap<String,String> ids) {
		ArrayList<String> idlist = new ArrayList<String>();
		for(String altid : ALT_IDS) 
			if(ids.containsKey(altid)) 
				idlist.add(altid.replace(" ", "-")+":"+ids.get(altid));
		Collections.sort(idlist,String.CASE_INSENSITIVE_ORDER);
		String idstr = StringUtils.join(idlist,DELIM);
		if (idstr.isEmpty())
			idstr = "None";
		return idstr;
	}
	
	/**
	 * Return a concatenated sting of Feature references.
	 * @param features
	 * @return
	 */
	public String getFeatures(Set<EntityFeature> features) {
		ArrayList<String> featurelist = new ArrayList<String>();
		for(EntityFeature f : features) {
			if (f.toString().contains("residue modification, active"))
				featurelist.add("activating-residue-modification");
			else if (f.toString().contains("residue modification, inactive"))
				featurelist.add("inactivating-residue-modification");
			else
				featurelist.add(f.toString().replace(' ','-'));
		}
		Collections.sort(featurelist,String.CASE_INSENSITIVE_ORDER);
		String featurestr = StringUtils.join(featurelist,DELIM);
		if (featurestr.isEmpty())
			featurestr = "None";
		return featurestr;
	}
	
	/**
	 * Return a concatenated string of Evidnece references
	 * @param evidence
	 * @return
	 */
	public String getEvidence(Set<Evidence> evidence) {
		ArrayList<String> evlist = new ArrayList<String>();
		String evname,tmp;
		for(Evidence e: evidence) { // set of Evidence Objects
			for(EvidenceCodeVocabulary ec : e.getEvidenceCode()) { // set of evidence code vocabs
				evname = ec.toString().replace(' ','-');
				if (ec.getXref().size() == 0) {
					evlist.add(evname);
				} else {
					tmp = evname;
					for (Xref xref: ec.getXref())
						tmp+=":"+xref.getId();
					evlist.add(tmp);
				}
			}
		}
		Collections.sort(evlist,String.CASE_INSENSITIVE_ORDER);
		String evstr = StringUtils.join(evlist,DELIM);
		if (evstr.isEmpty())
			evstr = "None";
		return evstr;
	}
	
	/**
	 * Return a concatenated string of the X-references.
	 * @param xrefs 
	 * @return xrefstring or "none" if there are no x references.
	 */
	public String getXrefs(Set<Xref> xrefs) {
		ArrayList<String> xreflist = new ArrayList<String>();
		for(Xref xref : xrefs) {
			if (xref.getDb() == null || xref.getId()==null)
				continue;
			if (xref.getDb().equals("Reactome") && xref.getId().contains("REACT_")) {
				xreflist.add((xref.getDb()+"Link:http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID="+xref.getId()).replace(' ','-'));
			} else
				xreflist.add((xref.getDb()+":"+xref.getId()).replace(' ','-'));
		}
		Collections.sort(xreflist,String.CASE_INSENSITIVE_ORDER);
		String xrefstr = StringUtils.join(xreflist,DELIM);
		if (xrefstr.equals(""))
			xrefstr = "None";
		return xrefstr;
	}
	
	/**
	 * Return a concatenated string of stoichiometry of the form "pid:float".
	 * This may be called with either participantStoichiometry() or componentStoichiometry().
	 * @param stoic
	 * @return stroicstr or "none" if there is no stoichiometry information.
	 */
	public String getStoichiometry(Set<Stoichiometry> stoic) {
		ArrayList<String> stoicarray = new ArrayList<String>();
		for(Stoichiometry s: stoic) {
			stoicarray.add(getID(s.getPhysicalEntity())+":"+s.getStoichiometricCoefficient()+"");
		}
		Collections.sort(stoicarray,String.CASE_INSENSITIVE_ORDER);
		String stoicstr = StringUtils.join(stoicarray,DELIM);
		if (stoicstr.isEmpty())
			stoicstr = "None";
		return stoicstr;
	}
	/**
	 * Converts a set of components to a string of comma-separated IDs.
	 * Sorts this list alphabetically.
	 *  
	 * @param components - components to convert
	 * @return - String of a comma-separated, sorted IDs
	 * 
	 */
	public String convertList(Set<?> components) {

		// If the list is empty, return 'None'.
		if (components.size()==0)
			return "None";

		// Make a list of the IDs (strings)
		ArrayList<String> newlist = new ArrayList<String>(components.size());
		for (BioPAXElement entity : (Set<BioPAXElement>) components) 
			newlist.add(getID(entity));

		// Sort the list
		Collections.sort(newlist,String.CASE_INSENSITIVE_ORDER);

		// Join on ';'
		String strToReturn = StringUtils.join(newlist,DELIM);

		return strToReturn;
	}

	/* ************************** */
	/* Access and Utility Methods */
	/* ************************** */

	/**
	 * For a PhysicalEntity, search Xref and Reference objects to get all
	 * possible mappings in different namespaces.
	 * 
	 * @param entity - element to get mappings for
	 * @return HashMap of NameSpace to ID for this element 
	 */
	public HashMap<String,String> getAllIDs(PhysicalEntity entity) {

		// Initialized HashMap object for storing mapped IDs
		HashMap<String,String> mappedIDs = new HashMap<String,String>();

		// For all Xref objects, store namespace (DB) and ID (Id).
		// This is usually the REACTOME_DB and stable Reactome Identifier.
		for(Xref xref : entity.getXref()) 
			mappedIDs.put(xref.getDb(),xref.getId());

		// Other namespaces are stored in EntityReferences.
		// To get the entity reference, need to determine the type
		// of PhysicalEntity.  For example, the entity might have a 
		// ProteinReference, SmallMoleculeReference, etc.  These 
		// are all implementations of the EntityReference interface.
		//
		// If the entity reference is null, Reactome might list ANOTHER 
		// PhysicalEntity its a member of - this parent element might contain
		// the proper EntityReference.  For example, different phosphorylated
		// forms of beta-catenin are all members of 'beta-catenin', which 
		// contains the proper entity reference for UniProt mapping.
		EntityReference e = null;
		if (entity instanceof Protein) {
			e = ((Protein)entity).getEntityReference();
			if (database == DB.REACTOME && e == null) 
				for(PhysicalEntity e2 : ((Protein)entity).getMemberPhysicalEntity())
					mappedIDs = addEntityReference(((Protein)e2).getEntityReference(),mappedIDs);
		} else if (entity instanceof SmallMolecule) {
			e = ((SmallMolecule)entity).getEntityReference();
			if (database == DB.REACTOME && e == null) 
				for(PhysicalEntity e2 : ((SmallMolecule)entity).getMemberPhysicalEntity())
					mappedIDs = addEntityReference(((SmallMolecule)e2).getEntityReference(),mappedIDs);
		} else if (entity instanceof Dna) {
			e = ((Dna)entity).getEntityReference();
			if (database == DB.REACTOME && e == null) 
				for(PhysicalEntity e2 : ((Dna)entity).getMemberPhysicalEntity())
					mappedIDs = addEntityReference(((Dna)e2).getEntityReference(),mappedIDs);
		} else if (entity instanceof Rna) {
			e = ((Rna)entity).getEntityReference();
			if (database == DB.REACTOME && e == null) 
				for(PhysicalEntity e2 : ((Rna)entity).getMemberPhysicalEntity())
					mappedIDs = addEntityReference(((Rna)e2).getEntityReference(),mappedIDs);
		}

		// if the Entity Reference is not null, then one of the getEntityReference()
		// method calls returned a proper EntityReference.  Add this to mappedIDs.
		if (e!=null)
			mappedIDs = addEntityReference(e,mappedIDs);
		
		return mappedIDs;
	}

	/**
	 * This takes a HashMap of Namespace to IDs and an entity reference, and returns the
	 * HashMap with the new references added.
	 * 
	 * @param e - EntityReference to add
	 * @param mappedIDs - HashMap of Namespace to IDs
	 * @return HashMap with newly-added mapping information
	 * 
	 */
	public HashMap<String,String> addEntityReference(EntityReference e, HashMap<String,String> mappedIDs) {

		// Iterate through "Names"
		for(String n: e.getName()) {
			// If the name contains no colons, it is assumed to be a 
			// standard (common) name.  Add it to the key "Name", which
			// combines all of the common names in a comma-separated list.
			if (!n.contains(":"))
				if(!mappedIDs.containsKey("Name"))
					mappedIDs.put("Name",n.replace(' ','-'));
				else if(!mappedIDs.get("Name").contains(n.replace(' ','-')))
					mappedIDs.put("Name",mappedIDs.get("Name")+DELIM+n.replace(' ','-'));
		}
		
		// For all Xref objects for this EntityReference, 
		// store namespace (DB) and ID (Id).
		for(Xref xref : e.getXref()) {
			mappedIDs.put(xref.getDb(),xref.getId());
		}
		
		// NOTE: NCI-PID may have additional Entity References within an Entity Reference;
		// this is an Entity Set.  This condition is checked in the writeBioPAXElement() 
		// function and handled in the writeEntitySet() function.
		
		return mappedIDs;
	}

	/**
	 * Strip the prefix of the RDFId() for the element.
	 * 
	 * @param bpelement
	 * @return the suffix of the RDFId split on '#'
	 * 
	 */
	public String getID(BioPAXElement bpelement) {
		String idstr = "";
		if (database == DB.REACTOME) 
			idstr = bpelement.getRDFId().split("#")[1];
		else if (database == DB.NETPATH) 
			idstr = bpelement.getRDFId();
		else if (database == DB.CELLDESIGNER) {
			String id = bpelement.getRDFId();
			if (id.contains("#"))
				idstr = id.split("#")[1];
			else
				idstr = id;
		}
		else if (database == DB.NCIPID) {
			idstr = bpelement.getRDFId().split("biopax")[1];
			if (idstr.contains("#"))
				idstr = idstr.split("#")[1];
		} else if (database == DB.PC) {
			idstr = bpelement.getRDFId();
		}
		if (idstr.isEmpty()) {
			System.out.println("ERROR: " + bpelement + " does not conform to RDFID syntax (pid #).");
			System.out.println(getID(bpelement));
			System.out.println(bpelement.getRDFId());
			System.exit(-1);
		}
		
		return idstr;
	}

	/**
	 * Return the properties of a BioPAXElement.  
	 * Not used in this parser, but potentially useful; see paxtools.pdf 
	 * 
	 * @param bpe - BioPAXElement
	 * @return array of properties for that BioPAXElement
	 * 
	 */
	public String[] listProperties(BioPAXElement bpe) {
		EditorMap editorMap = SimpleEditorMap.L3;
		Set<PropertyEditor> editors = editorMap.getEditorsOf(bpe);
		String[] properties = new String[editors.size()];
		int i=0;
		for (PropertyEditor editor : editors) {
			properties[i]=editor.getProperty();
			i++;
		}
		return properties;
	}	 
}
