package main;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.biopax.paxtools.controller.PathAccessor;
import org.biopax.paxtools.controller.SimpleEditorMap;
import org.biopax.paxtools.io.BioPAXIOHandler;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level2.pathway;
import org.biopax.paxtools.model.level3.BiochemicalReaction;
import org.biopax.paxtools.model.level3.Pathway;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.model.level3.Process;
import org.biopax.paxtools.query.QueryExecuter;
import org.biopax.paxtools.query.algorithm.Direction;

// API: http://biopax.github.io/Paxtools/4.3.1/apidocs/


public class GetPathwayMembers {
	String infile,outfile,entityFile;
	
	public GetPathwayMembers(String[] args) {
		if(args.length != 2 && args.length != 3) {
			System.out.println("USAGE: java StreamSurvey <OWL File> <Output File> <Entities_OPTIONAL>");
			System.exit(0);
		}
		
		infile = args[0];
		outfile = args[1];
		if (args.length == 3) 
			entityFile = args[2];
		else
			entityFile = null;
		System.out.println("Infile: "+infile);
		System.out.println("Outfile: "+outfile);
		System.out.println("Entity File: " + entityFile);
	}

	public static void main(String[] args) throws IOException, NoSuchFieldException {
		
		GetPathwayMembers SS = new GetPathwayMembers(args);
				
		BioPAXIOHandler handler = new SimpleIOHandler(BioPAXLevel.L3);
		Model model = handler.convertFromOWL(new FileInputStream(SS.infile));
		SS.getStats(model);
		
		// SIGNAL TRANSDUCTION: http://identifiers.org/reactome/R-HSA-162582
		// OTHER PATHWYAS
		HashSet<String> pathwaysToExpand = new HashSet<String>();
		pathwaysToExpand.add("http://identifiers.org/reactome/R-HSA-9006934"); // Signaling by Receptor Tyrosine Kinases
		pathwaysToExpand.add("http://identifiers.org/reactome/R-HSA-9006936"); // Signaling by TGF Beta Family Members
		pathwaysToExpand.add("http://identifiers.org/reactome/R-HSA-73887"); // Death Receptor Signaling
		pathwaysToExpand.add("http://identifiers.org/reactome/R-HSA-9006925"); //Intracellular signaling by second messengers
		pathwaysToExpand.add("http://identifiers.org/reactome/R-HSA-9006927"); //	Signaling by Non-Receptor Tyrosine Kinases 
		pathwaysToExpand.add("http://identifiers.org/reactome/R-HSA-5683057");	// MAPK family signaling cascades
		
		Pathway p = (Pathway) model.getByID("http://identifiers.org/reactome/R-HSA-162582");
		
		HashSet<Pathway> pathways = new HashSet<Pathway>();
		for (Process process : p.getPathwayComponent()) {
			if (process instanceof Pathway){
				//System.out.println(process.getUri());
				if (pathwaysToExpand.contains(process.getUri())) {
					for (Process subprocess: ((Pathway)process).getPathwayComponent()) 
						if (subprocess instanceof Pathway) 
							pathways.add((Pathway)subprocess);
				} else { 
					pathways.add((Pathway)process);
				}
			}
		}
		System.out.println(pathways.size()+" pathways");
		
		BufferedWriter out = new BufferedWriter(new FileWriter(SS.outfile));
		out.write("#Pathway\tName\tNumMembers\tMembers\n");
		HashSet<String> entityFilter = SS.getFilteredEntities(SS.entityFile);
		
		PathwayMemberVisitor visitor = new PathwayMemberVisitor(SimpleEditorMap.L3,model);
				
		// Parse each pathway
		int i = 1;
		for (Pathway pathway : pathways) {
			System.out.println(i + " of " + pathways.size());
			i++;
			visitor.reset();
			visitor.traverseAndCount(pathway);
			System.out.println(visitor.visitedObjects.size());
			HashSet<String> result = new HashSet<String>(visitor.visitedObjects);
			result.retainAll(entityFilter);
				
			out.write(pathway.getUri()+"\t"+pathway.getDisplayName()+"\t"+result.size()+"\t"+String.join(";",result)+"\n");
		}
		
		out.close();
		System.out.println("Wrote  pathways to " + SS.outfile);
	}

	public HashSet<String> getFilteredEntities(String entityFile) throws IOException {
		if(entityFile == null)
			return null;
		HashSet<String> entityFilter = new HashSet<String>();
		try (BufferedReader br = new BufferedReader(new FileReader(entityFile))) {
		    String line;
		    while ((line = br.readLine()) != null) {
		       if (line.startsWith("#"))
		    	   continue;
		       entityFilter.add(line.split("\t")[0]);
		       //System.out.println(line.split("\t")[0]);
		    }
		}
		System.out.println("Read " + entityFilter.size() + " entities");
		return entityFilter;	
	}
	
	
	public void getStats(Model model) {
		System.out.println("BioPAX Model ");
		System.out.println("  "+ model.getObjects(PhysicalEntity.class).size() + " Physical Entities");
		System.out.println("  "+ model.getObjects(BiochemicalReaction.class).size() + " Biochemical Reactions");
		System.out.println("  "+ model.getObjects(Pathway.class).size() + " Pathways");
	}
}