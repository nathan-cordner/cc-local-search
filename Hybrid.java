import java.util.*;

public class Hybrid {

    public static ArrayList<ArrayList<Integer>> large_graph_fix_clusters_local(ArrayList<ArrayList<Integer>> cur_clustering, ArrayList<ArrayList<Integer>> prob_matrix, boolean ECS) {
    /*
        Apply local search to clusters resulting from Pivot

    */

    ArrayList<ArrayList<Integer>> new_clustering = new ArrayList<ArrayList<Integer>>();
    int negative_examples = 0;
    int total_examples = 0;
    for (int c = 0; c < cur_clustering.size(); c++) {
        ArrayList<Integer> cluster = cur_clustering.get(c);
        int c_size = cluster.size();

        

        if (c_size <= 3) {
            new_clustering.add(cluster);
        } else { 

            HashMap<Integer, Integer> relabelling = new HashMap<Integer, Integer>();
            for (int i = 0; i < c_size; i++)
                relabelling.put(cluster.get(i), i);

            HashSet<Integer> cluster_set = new HashSet<Integer>();
            cluster_set.addAll(cluster);
            ArrayList<ArrayList<Integer>> cur_prob_matrix = new ArrayList<ArrayList<Integer>>();
	    ArrayList<Integer> cluster_relabel = new ArrayList<Integer>();

            //int half = 0;
            //int total = 0;
            for (int i = 0; i < c_size; i++) {
                ArrayList<Integer> new_edges = new ArrayList<Integer>();
                ArrayList<Integer> cur_edges = prob_matrix.get(cluster.get(i));
                int edge_count = 0;
                for (int j = 0; j < cur_edges.size(); j++){
                    int cur_edge = cur_edges.get(j);
                    if (cluster_set.contains(cur_edge)) {
                        new_edges.add(relabelling.get(cur_edge));
                        edge_count += 1;
                    }
                }
                //if (edge_count < c_size / 2)
                //    half += 1;
                //total += edge_count;
                cur_prob_matrix.add(new_edges);
		cluster_relabel.add(i);
            }
            /*
            double precision = ((double) total) / (c_size * (c_size -1));
            double half_ratio = ((double) half) / c_size;
            
            if (c_size >= 10 && precision > 0.33)
                total_examples += 1;               

            if (c_size >= 10 && precision > 0.33 && half_ratio < 1 - Math.sqrt(precision)) {
                System.out.println("Cluster size: " + c_size);
                System.out.println("Cluster precision: " + precision);
                System.out.println("Number of nodes with degree < n/2: " + half_ratio); 
                negative_examples += 1;
            }
            */

            ArrayList<ArrayList<Integer>> input_clustering = new ArrayList<ArrayList<Integer>>();
            input_clustering.add(cluster_relabel);
	    // System.out.println(input_clustering.size());
            ArrayList<ArrayList<Integer>> adjusted_clustering = DNode.local_search_network(input_clustering, cur_prob_matrix, ECS);

            // replace node labels
            for (int j = 0; j < adjusted_clustering.size(); j++) {
                ArrayList<Integer> adj_cluster = adjusted_clustering.get(j);
                for (int i = 0; i < adj_cluster.size(); i++)
                    adj_cluster.set(i, cluster.get(adj_cluster.get(i)));
                new_clustering.add(adj_cluster);
            }

        } 
    }
    // System.out.println("Probability of negative example: " + ((double) negative_examples) / total_examples);


    return new_clustering;

    }


        // INPUT: cluser labels, will need to compute |C|^2 probabilities for each cluster
        public static ArrayList<ArrayList<Integer>> large_graph_fix_clusters_consensus(ArrayList<ArrayList<Integer>> cur_clustering, ArrayList<HashMap<Integer, Integer>> cluster_labels, int MAX_ATTRIBUTES, boolean useRandom) {
            /*
                Apply DeterministicNode to clusters resulting from pKwik algorithm    
            */

            MAX_ATTRIBUTES = Math.min(MAX_ATTRIBUTES, cluster_labels.size()); // don't want max attributes to exceed num input clusterings

            ArrayList<ArrayList<Integer>> new_clustering = new ArrayList<ArrayList<Integer>>();
            for (int c = 0; c < cur_clustering.size(); c++) {
                ArrayList<Integer> cluster = cur_clustering.get(c);
                int c_size = cluster.size();
        
                if (c_size <= 2) { // avoid hybrid on clusters that are too small or too large
                    new_clustering.add(cluster);
                } else { // deal with the memory issues by computing probabilites on the fly 
                    ArrayList<ArrayList<Integer>> adjusted_clustering = DNode.local_search_consensus(cluster_labels, cluster, MAX_ATTRIBUTES, useRandom);
        
                    // CHECK IF NEW CLUSTERING GIVES BETTER DISAGREEMENT DISTANCE:
		    ArrayList<ArrayList<Integer>> old_clustering = new ArrayList<ArrayList<Integer>>();
		    old_clustering.add(cluster);
		    long old_cost = Consensus.getSubsetEditDist(old_clustering, cluster, cluster_labels);
		    long new_cost = Consensus.getSubsetEditDist(adjusted_clustering, cluster, cluster_labels);
		    if (new_cost < old_cost) {
                        // replace node labels
                        for (int j = 0; j < adjusted_clustering.size(); j++) {
                            ArrayList<Integer> adj_cluster = adjusted_clustering.get(j);
                            new_clustering.add(adj_cluster);
                        }
		   } else {
                       new_clustering.add(cluster);
                   }

                }
            }
            return new_clustering;
        
        }

    
}
