package main;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.biopax.paxtools.controller.EditorMap;
import org.biopax.paxtools.controller.PropertyEditor;
import org.biopax.paxtools.controller.SimpleEditorMap;
import org.biopax.paxtools.controller.Traverser;
import org.biopax.paxtools.controller.Visitor;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.level3.Process;

/**
 * This visitor class is used to traverse BioPAX files.
 * See the Visitor class in PAXTools.
 *
 */
public class PathwayMemberVisitor implements Visitor {

	/* Class Variables */
	
	// Current pathway being traversed
	public static Pathway pathway = null;

	// origModel is the original BioPAX object that is traversed.
	private Model origModel;

	// Objects for the Visitor interface. See biopax.pdf for an example.
	private Traverser traverser;
	private EditorMap editorMap;

	// Keep track of the number of properties and objects we visit.
	public HashSet<String> visitedObjects;

	/* ***************************** */
	/* Constructor and Reset Methods */
	/* ***************************** */

	/**
	 * Constructor.
	 * @param em - EditorMap
	 * @param om - Original Model
	 * @param prefix - output file prefix
	 */
	public PathwayMemberVisitor(EditorMap em, Model om) {
		editorMap = em;
		origModel = om;
		traverser = new Traverser(editorMap,this);

		// Initialize variables for tracking objects and properties.
		visitedObjects = new HashSet<String>();
	}

	/**
	 * Resets the variables for tracking objects and properties.
	 */
	public void reset() {
		visitedObjects.clear();
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
		if (visitedObjects.contains(bpelement.getUri()))
			return;
			

		// Add this element to the list of visited objects.
		visitedObjects.add(bpelement.getUri());

		// Recurse. Note this is ONLY if we haven't seen the element yet.
		traverser.traverse(bpelement, model);

	}

	
	/**
	 * Traverses a BioPAX object, calling the visit() method 
	 * for each traversed element.
	 * 
	 * @param pathway - BioPAX Pathway to traverse
	 * @throws IOException - if output files cannot be written.
	 * @throws NoSuchFieldException 
	 * 
	 */
	public void traverseAndCount(Pathway pathway) throws IOException, NoSuchFieldException {
		
		// Reset counts (allows us to call traverse() on different pathways)
		reset();
		System.out.println("---------------------------------------------------------------------------------");
		System.out.println(pathway.getDisplayName());
		traverser.traverse(pathway,origModel);
	}
}
