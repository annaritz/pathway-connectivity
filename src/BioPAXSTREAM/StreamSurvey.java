import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
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
import org.biopax.paxtools.model.level2.pathway;
import org.biopax.paxtools.model.level3.BiochemicalReaction;
import org.biopax.paxtools.model.level3.Pathway;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.query.QueryExecuter;
import org.biopax.paxtools.query.algorithm.Direction;

// API: http://biopax.github.io/Paxtools/4.3.1/apidocs/


public class StreamSurvey {
	String infile,outfile;
	
	public StreamSurvey(String[] args) {
		if(args.length != 2) {
			System.out.println("USAGE: java StreamSurvey <OWL File> <Output File>");
			System.exit(0);
		}
		
		infile = args[0];
		outfile = args[1];
		System.out.println("Infile: "+infile+"\nOutfile: "+outfile);
	}

	public static void main(String[] args) throws IOException {
		
		StreamSurvey SS = new StreamSurvey(args);
				
		BioPAXIOHandler handler = new SimpleIOHandler(BioPAXLevel.L3);
		Model model = handler.convertFromOWL(new FileInputStream(SS.infile));
		SS.getStats(model);
		HashMap<String,String> id2name = SS.getDownstream(model);
		String namefile = SS.outfile+".names";
		BufferedWriter out = new BufferedWriter(new FileWriter(namefile));
		for (String key: id2name.keySet()) {
			out.write(key+"\t"+id2name.get(key)+"\n");
		}
		out.close();
		System.out.println("Wrote "+id2name.size()+" names to " + namefile);
	}

	public HashMap<String,String> getDownstream(Model model) throws IOException {
		Set<PhysicalEntity> entitySet = model.getObjects(PhysicalEntity.class);
		System.out.println(entitySet.size() + " entity sets.");
		int size = entitySet.size();
		int i = 0;
		HashMap<String,String> id2name = new HashMap<String,String>();
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
		out.write("#Entity\tNumDownstream\n");
		
		int limit = 30;
		Direction downstream = Direction.DOWNSTREAM;
		Set<BioPAXElement> source;
		for (PhysicalEntity entity: entitySet) {
			i++;
			if (i % 100 == 0){
				System.out.println(" "+i+" of "+size);
			}
			id2name.put(entity.getRDFId(),entity.getDisplayName());
			source = new HashSet<BioPAXElement>();
			source.add(entity);
			Set<BioPAXElement> result = QueryExecuter.runNeighborhood(source, model, limit, downstream);
			//System.out.println("Entity: "+entity.getDisplayName()+": "+result.size()+" downstream elements");
			Set<String> members = new HashSet<String>();
			for (BioPAXElement e : result) {
				members.add(entity.getRDFId());
				/**
				if (e instanceof PhysicalEntity) 
					System.out.println("  "+((PhysicalEntity)e).getDisplayName());
				else
					System.out.println("  "+e.getUri());
				**/
			}
			//out.write(entity.getUri()+"\t"+result.size()+"\t"+String.join(",", members)+"\n");
			out.write(entity.getRDFId()+"\t"+result.size()+"\n");
			
			//System.exit(0);
		}
		
		out.close();
		System.out.println("Wrote to " + outfile);
		return id2name;
	}
	
	public void getStats(Model model) {
		System.out.println("BioPAX Model ");
		System.out.println("  "+ model.getObjects(PhysicalEntity.class).size() + " Physical Entities");
		System.out.println("  "+ model.getObjects(BiochemicalReaction.class).size() + " Biochemical Reactions");
		System.out.println("  "+ model.getObjects(Pathway.class).size() + " Pathways");
	}
}
