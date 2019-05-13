package main;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.biopax.paxtools.io.BioPAXIOHandler;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.BiochemicalReaction;
import org.biopax.paxtools.model.level3.Pathway;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.query.QueryExecuter;
import org.biopax.paxtools.query.algorithm.Direction;

// API: http://biopax.github.io/Paxtools/4.3.1/apidocs/


public class StreamSurvey2 {
	String infile,outfile,entityFile;
	int limit;
	
	public StreamSurvey2(String[] args) {
		if(args.length != 3 && args.length != 4) {
			System.out.println("USAGE: java StreamSurvey <OWL File> <Output File> <LIMIT> <Entities_OPTIONAL>");
			System.exit(0);
		}
		
		infile = args[0];
		outfile = args[1];
		limit = Integer.parseInt(args[2]);
		if (args.length == 4) 
			entityFile = args[3];
		else
			entityFile = null;
		System.out.println("Infile: "+infile);
		System.out.println("Outfile: "+outfile);
		System.out.println("Limit: " + limit);
		System.out.println("Entity File: " + entityFile);
	}

	public static void main(String[] args) throws IOException {
		
		StreamSurvey2 SS = new StreamSurvey2(args);
				
		BioPAXIOHandler handler = new SimpleIOHandler(BioPAXLevel.L3);
		Model model = handler.convertFromOWL(new FileInputStream(SS.infile));
		SS.getStats(model);
		
		HashSet<String> entityFilter = SS.getFilteredEntities(SS.entityFile);
		HashSet<String> intersection = null;
		if (entityFilter != null) {
			intersection = SS.findMissingElements(model, entityFilter);
			if(intersection.size() != 0) {
				System.out.println(intersection.size() + " elements are IN filter but NOT a physical element of BioPAX:");
				for (String s: intersection) 
					System.out.println("  "+s);
				System.exit(0);
			} else
				System.out.println("All elements IN filter are IN OWL file. Continuing.");
		}
		
		HashMap<String,String> id2name = SS.getDownstream(model,SS.limit,entityFilter);
		String namefile = SS.outfile+".names";
		BufferedWriter out = new BufferedWriter(new FileWriter(namefile));
		for (String key: id2name.keySet()) {
			out.write(key+"\t"+id2name.get(key)+"\n");
		}
		out.close();
		System.out.println("Wrote "+id2name.size()+" names to " + namefile);
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
	
	public HashMap<String,String> getDownstream(Model model, int limit, HashSet<String> entityFilter) throws IOException {
		Set<PhysicalEntity> entitySet = model.getObjects(PhysicalEntity.class);
		HashMap<String,String> id2name = new HashMap<String,String>();
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
		out.write("#Entity\tNumDownstream\tDownstreamElements\n");
	
		Direction downstream = Direction.DOWNSTREAM;
		Set<BioPAXElement> source;
		int i = 0;
		int size = entitySet.size();
		if (entityFilter != null) 
			size = entityFilter.size();
		int skipped = 0;
		Set<String> found = new HashSet<String>();
		for (PhysicalEntity entity: entitySet) {
			if(entityFilter != null && !entityFilter.contains(entity.getUri())) {
				skipped++;
				continue;
			} else if (found.contains(entity.getUri()))
				System.out.println("Skipping already-added entity " + entity.getUri());
			i++;
			found.add(entity.getUri());
			if (i % 100 == 0){
				System.out.println(" "+i+" of "+size + "; " + skipped + " skipped.");
			}
			id2name.put(entity.getUri(),entity.getDisplayName());
			source = new HashSet<BioPAXElement>();
			source.add(entity);
			Set<BioPAXElement> result = QueryExecuter.runNeighborhood(source, model, limit, downstream);
			Set<String> members = new HashSet<String>();
			for (BioPAXElement e : result) {
				if(entityFilter == null || entityFilter.contains(e.getUri())) 
					members.add(e.getUri());
			}
			//if (result.size() != members.size()) 
			//	System.out.println("Entity: "+entity.getDisplayName()+": "+result.size()+" downstream elements, "+members.size() + " in filtered set");
			
			out.write(entity.getUri()+"\t"+members.size()+"\n");
			
			//System.exit(0);
		}
		
		out.close();
		System.out.println("Wrote to " + outfile);
		return id2name;
	}
	
	public HashSet<String> findMissingElements(Model model, HashSet<String> entityFilter) {
		Set<PhysicalEntity> entitySet = model.getObjects(PhysicalEntity.class);
		HashSet<String> entityNames = new HashSet<String>();
		for (PhysicalEntity e : entitySet) {
			entityNames.add(e.getUri());
		}
		HashSet<String> intersection = new HashSet<String>(entityFilter);
		intersection.removeAll(entityNames);
		return intersection;
	}
	
	public void getStats(Model model) {
		System.out.println("BioPAX Model ");
		System.out.println("  "+ model.getObjects(PhysicalEntity.class).size() + " Physical Entities");
		System.out.println("  "+ model.getObjects(BiochemicalReaction.class).size() + " Biochemical Reactions");
		System.out.println("  "+ model.getObjects(Pathway.class).size() + " Pathways");
	}
}