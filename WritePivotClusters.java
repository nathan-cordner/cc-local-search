import java.util.*;
import java.io.*;

public class WritePivotClusters {

    private static void writeConsensusGraph(ArrayList<HashMap<Integer,Integer>> cluster_labels, int num_nodes) {

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter("clusterings_graph.txt"));
            

            for (int i = 0; i < num_nodes; i++) {
                StringBuilder cur_line = new StringBuilder(2*cluster_labels.size());
                cur_line.append(cluster_labels.get(0).get(i));

                for (int j = 1; j < cluster_labels.size(); j++) {
                    cur_line.append("," + cluster_labels.get(j).get(i));
                }
                writer.write(cur_line.toString());
                writer.newLine();

                if ((i + 1) % 1000 == 0)
                    System.out.println("Finished writing " + (i + 1) + " lines");
            }            

            writer.close();
        } catch (Exception e) {
            System.out.println("Something went wrong. :(");
        }
    }

    private static ArrayList<ArrayList<Integer>> getGroundTruth(String data_set, String delimiter) {

        String file_name = "Data/"+ data_set + "/graph.txt";
        HashMap<Integer, Integer> my_mapping = new HashMap<Integer, Integer>();

        try {
            BufferedReader myReader = new BufferedReader(new FileReader(file_name));
            
            String data = myReader.readLine();
            String[] split = data.split(" ");
            int num_pts = Integer.parseInt(split[0]);

            int cur_index = 0;
            while ((data = myReader.readLine()) != null) {

                split = data.split(delimiter);
                int x = Integer.parseInt(split[0]);
                int y = Integer.parseInt(split[1]);

                // RELABEL
                if (!my_mapping.containsKey(x)) {
                    my_mapping.put(x, cur_index);
                    cur_index += 1;
                }    
                if (!my_mapping.containsKey(y)) {
                    my_mapping.put(y, cur_index);
                    cur_index += 1; 
                }

            }
            myReader.close();

            BufferedReader myReader2 = new BufferedReader(new FileReader("Data/" + data_set + "/ground_truth.txt"));

            //ArrayList<ArrayList<Integer>> my_clustering = new ArrayList<ArrayList<Integer>>();
            HashMap<Integer, ArrayList<Integer>> map_clustering = new HashMap<Integer, ArrayList<Integer>>();
			    
            cur_index = 0;
            while ((data = myReader2.readLine()) != null) {
		// answers are ordered by question!!
                split = data.split(delimiter);
                int x = Integer.parseInt(split[0]);
		if (!map_clustering.containsKey(x)) {
                    ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
		    cur_cluster.add(my_mapping.get(cur_index));
		    map_clustering.put(x, cur_cluster);
		} else {
                    map_clustering.get(x).add(my_mapping.get(cur_index));
		}
                
		cur_index += 1;

            }

            myReader2.close();
	    ArrayList<ArrayList<Integer>> my_clustering = new ArrayList<ArrayList<Integer>>();
	    for (Integer x : map_clustering.keySet()) {
                my_clustering.add(map_clustering.get(x));
	    }

            return my_clustering;

          } catch (Exception e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            return null;
          }


        // return null;
    }
    

    public static void main(String args[]){

        // ----- PARAMETERS -----
        String data_set = args[0];
        String delimiter = ","; 
        int ROUNDS = 100;
        // ----------------------

        ArrayList<ArrayList<Integer>> prob_matrix = Helper.read_large_network_relabel("Data/"+ data_set + "/graph.txt", delimiter);
        ArrayList<ArrayList<ArrayList<Integer>>> pivot_clusterings = new ArrayList<ArrayList<ArrayList<Integer>>>();
	ArrayList<ArrayList<Integer>> ground_truth = getGroundTruth(data_set, delimiter);

        long pivot_disagreement = 0;
	long hybrid_disagreement = 0;


        System.out.println("Start");
        for (int j = 0; j < ROUNDS; j++) {

            ArrayList<ArrayList<Integer>> result = PKwik.pKwikClustering(prob_matrix);
            pivot_clusterings.add(result);
	    // long pCost = Consensus.quickEditDist(result, ground_truth);
	    // pivot_disagreement += pCost;
	    // System.out.print(pCost);

	    // ArrayList<ArrayList<Integer>> fix = Hybrid.large_graph_fix_clusters(result, prob_matrix);
	    // long hCost = Consensus.quickEditDist(fix, ground_truth);
	    // hybrid_disagreement += hCost;
	    // System.out.print(" " + hCost);
	    // System.out.println();

        }
        System.out.println("Finish");
        System.out.println();

        double pivotAvg = pivot_disagreement / ((double) ROUNDS);
	double hybridAvg = hybrid_disagreement / ((double) ROUNDS);

	System.out.println("Average Pivot: " + pivotAvg);
	System.out.println("Average Hybrid: " + hybridAvg);

        // ArrayList<ArrayList<Integer>> ground_truth = getGroundTruth(data_set, delimiter);
	// for (ArrayList<ArrayList<Integer>> pivot_clustering : pivot_clusterings) {
        //     System.out.println(Consensus.quickEditDist(pivot_clustering, ground_truth));
        // }


	    // ----------------------
        Consensus myCon = Consensus.consensusLabels(pivot_clusterings);
	writeConsensusGraph(myCon.getClusterLabels(), myCon.getNumNodes());


    }
    
}
