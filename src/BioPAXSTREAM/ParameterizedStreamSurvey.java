package main;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import org.biopax.paxtools.io.BioPAXIOHandler;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;

public class ParameterizedStreamSurvey {

	public String infile = "/Users/annaritz/Documents/github/pathway-connectivity/data/OWL/PathwayCommons10.reactome.BIOPAX.owl";
	public String outdir = "/Users/annaritz/Documents/github/pathway-connectivity/src/BioPAXSTREAM/output/";
	public String filterfile = "/Users/annaritz/Documents/github/pathway-connectivity/src/hypergraph_code/output/reactome.txt";
	public int min_limit = 1;
	public int max_limit = 30;
	
	// Constructor
	public ParameterizedStreamSurvey() throws IOException {
		//java StreamSurvey2 <OWL File> <Output File> <LIMIT> <Entities_OPTIONAL>
		//EXAMPLE ARGS: 
		// /Users/annaritz/Documents/github/pathway-connectivity/data/OWL/PathwayCommons10.reactome.BIOPAX.owl 
		// /Users/annaritz/Documents/github/pathway-connectivity/src/BioPAXSTREAM/output/reactome_limit30_filtered.txt 
		// 30 
		// /Users/annaritz/Documents/github/pathway-connectivity/src/hypergraph_code/output/reactome.txt
		//END EXAMPLE ARGS
		
		// get StreamSurvey2 object
		String[] args = new String[4];
		args[0] = infile;
		args[1] = null; // override this
		args[2] = "0"; // override this
		args[3] = filterfile;
		
		StreamSurvey2 SS = new StreamSurvey2(args);
		BioPAXIOHandler handler = new SimpleIOHandler(BioPAXLevel.L3);
		Model model = handler.convertFromOWL(new FileInputStream(SS.infile));
		SS.getStats(model);
		
		HashSet<String> entityFilter = getEntityFilter(SS,model);
		System.out.println("RUNNING");
		run(SS,model,entityFilter);
		System.out.println("Done!");
	}
	
	public HashSet<String> getEntityFilter(StreamSurvey2 SS, Model model) throws IOException {
		HashSet<String> entityFilter = SS.getFilteredEntities(SS.entityFile);
		HashSet<String> intersection = null;
		if (entityFilter != null) {
			intersection = SS.findMissingElements(model, entityFilter);
			if(intersection.size() != 0) {
				// if running ENTIRE Reactome BioPAX, all elements from filter file should be in OWL file.
				System.out.println(intersection.size() + " elements are IN filter but NOT a physical element of BioPAX:");
				for (String s: intersection) 
					System.out.println("  "+s);
				System.exit(0);
			} else
				System.out.println("All elements IN filter are IN OWL file. Continuing.");
		}
		return entityFilter;
	}
	
	public void run(StreamSurvey2 SS, Model model, HashSet<String> entityFilter) throws IOException {
		for (int i=min_limit; i<= max_limit; i++) {
			// set outfile & limit
			// example outfile: /Users/annaritz/Documents/github/pathway-connectivity/src/BioPAXSTREAM/output/reactome-signaling-by-wnt_limit30_filtered.txt
			SS.limit = i;
			System.out.println("------------ " + SS.limit + " -------------");
			SS.outfile = outdir+"reactome_limit"+SS.limit+"_filtered.txt";
			System.out.println(SS.outfile);
			File tmp = new File(SS.outfile);
			if (tmp.exists()) {
				System.out.println("FILE EXISTS - not re-running.");
				continue; // skip rest of the loop.
			}
			HashMap<String,String> id2name = SS.getDownstream(model,SS.limit,entityFilter);
			String namefile = SS.outfile+".names";
			BufferedWriter out = new BufferedWriter(new FileWriter(namefile));
			for (String key: id2name.keySet()) {
				out.write(key+"\t"+id2name.get(key)+"\n");
			}
			out.close();
			System.out.println("Wrote "+id2name.size()+" names to " + namefile);
			System.out.println("------------ DONE WITH " + SS.limit + " -------------");
		}
	}
	
	public static void main(String[] args) throws IOException {
		new ParameterizedStreamSurvey();
	}
	
	
}
