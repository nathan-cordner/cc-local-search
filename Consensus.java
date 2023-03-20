import java.util.*;

public class Consensus {

    private ArrayList<HashMap<Integer, Integer>> cluster_labels;
    private int num_nodes;

    public Consensus(ArrayList<HashMap<Integer, Integer>> cluster_labels, int num_nodes) {
        this.cluster_labels = cluster_labels;
        this.num_nodes = num_nodes;
    }

    public ArrayList<HashMap<Integer, Integer>> getClusterLabels() {
        return cluster_labels;
    }

    public int getNumNodes() {
        return num_nodes;
    }


    private static String getLabel (ArrayList<HashMap<Integer, Integer>> cluster_labels, int i) {
        String cur_label = "";
        for (int j = 0; j < cluster_labels.size(); j++)
            cur_label += "c" + cluster_labels.get(j).get(i);
        return cur_label;
    }

    public static Consensus consensusLabels (ArrayList<ArrayList<ArrayList<Integer>>> input_clusterings) {

        int num_clusterings = input_clusterings.size();
        int num_nodes = 0;

        // step 1: compute cluster labels 
        ArrayList<HashMap<Integer, Integer>> cluster_labels = new ArrayList<HashMap<Integer, Integer>>();
        for (int i = 0; i < num_clusterings; i++)
            cluster_labels.add(new HashMap<Integer, Integer>());

        for (int i = 0; i < num_clusterings; i++) {
            ArrayList<ArrayList<Integer>> cur_clustering = input_clusterings.get(i);
            HashMap<Integer, Integer> cur_map = cluster_labels.get(i);
            for (int j = 0; j < cur_clustering.size(); j++) {
                ArrayList<Integer> cur_cluster = cur_clustering.get(j);
                for (int k = 0; k < cur_cluster.size(); k++) {
                    // if (cur_map.containsKey(cur_cluster.get(k)))
                    //    System.out.println("Already have: " + cur_cluster.get(k) + " in clustering " + i);
                    cur_map.put(cur_cluster.get(k), j);
                    if (i == 0) {
                        num_nodes +=1; // COMPUTE THIS ONCE
                    }
                    
                }
            }
        }

        return new Consensus(cluster_labels, num_nodes);

    }

    public static void attributeShuffle(ArrayList<HashMap<Integer, Integer>> cluster_labels, int MAX_ATTRIBUTES) {

        int num_attributes = cluster_labels.size();
        Random r = new Random();

        if (MAX_ATTRIBUTES < num_attributes / 2.0) {
            for (int ii = 0; ii < MAX_ATTRIBUTES; ii++) // independent samples each time
                Collections.swap(cluster_labels, ii , ii + r.nextInt(num_attributes - ii));
        } else {
            for (int ii = num_attributes - 1; ii >= num_attributes - MAX_ATTRIBUTES; ii--)
                Collections.swap(cluster_labels, ii, r.nextInt(ii + 1));
        }
    }


    public static Helper consensusMatrix (ArrayList<ArrayList<ArrayList<Integer>>> input_clusterings) {

        int num_clusterings = input_clusterings.size();
        int num_nodes = 0;

        // step 1: compute cluster labels 
        ArrayList<HashMap<Integer, Integer>> cluster_labels = new ArrayList<HashMap<Integer, Integer>>();
        for (int i = 0; i < num_clusterings; i++)
            cluster_labels.add(new HashMap<Integer, Integer>());

        for (int i = 0; i < num_clusterings; i++) {
            ArrayList<ArrayList<Integer>> cur_clustering = input_clusterings.get(i);
            HashMap<Integer, Integer> cur_map = cluster_labels.get(i);
            for (int j = 0; j < cur_clustering.size(); j++) {
                ArrayList<Integer> cur_cluster = cur_clustering.get(j);
                for (int k = 0; k < cur_cluster.size(); k++) {
                    cur_map.put(cur_cluster.get(k), j);
                    if (i == 0)
                        num_nodes +=1; // COMPUTE THIS ONCE
                }
            }
        }

        // figure out how to return which nodes are already clustered together
        ArrayList<Integer> my_nodes = new ArrayList<Integer>();
        // HashSet<<ArrayList<Integer>>  ignore_list = new HashSet<Integer>();
        HashMap<String, Integer> unique_labels = new HashMap<String, Integer>();
        ArrayList<ArrayList<Integer>> index_mapping = new ArrayList<ArrayList<Integer>>();
        int cur_index = 0;
        for (int i = 0; i < num_nodes; i++) {
            String cur_label = getLabel(cluster_labels, i);
            if (!unique_labels.containsKey(cur_label)) {
                my_nodes.add(i);
                unique_labels.put(cur_label, cur_index);
                ArrayList<Integer> bucket = new ArrayList<Integer>();
                bucket.add(i);
                index_mapping.add(bucket);
                cur_index += 1;
            } else {
                // put i in an existing bucket
                index_mapping.get(unique_labels.get(cur_label)).add(i);
            }
        }

     
        int new_nodes = my_nodes.size();
        System.out.println("Num Mergeable Nodes: " + (num_nodes - new_nodes));

        // TEST
        /*
        for (int i = 0; i < index_mapping.size(); i++) {
            if (index_mapping.get(i).size() > 1)
                System.out.println(index_mapping.get(i));
        }
        */

        /* OLD METHOD
        double [][] cur_prob_matrix = new double[num_nodes][num_nodes];
        for (int i = 0; i < num_nodes; i++) {
            // int cur_node = my_nodes.get(i);
            cur_prob_matrix[i][i] = 1;

            for (int j = i + 1; j < num_nodes; j++) {
                int num_matches = 0;
                // int new_node = my_nodes.get(j);
                for (int k = 0; k < num_clusterings; k++) {
                    if (cluster_labels.get(k).get(i).equals(cluster_labels.get(k).get(j))) {
                        num_matches += 1;
                    }
                }
                double cur_prob = ((double) num_matches) / num_clusterings;
                cur_prob_matrix[i][j] = cur_prob;
                cur_prob_matrix[j][i] = cur_prob;
            }
        }
        */
        
        // double [][] cur_prob_matrix = new double[new_nodes][new_nodes];
        // int pos_edges = 0;
        // form new probability matrix
        ArrayList<ArrayList<Double[]>> cur_prob_matrix = new ArrayList<ArrayList<Double[]>>();
        for (int i = 0; i < new_nodes; i++)
            cur_prob_matrix.add(new ArrayList<Double[]>());

        for (int i = 0; i < new_nodes; i++) {
            int cur_node = my_nodes.get(i);

            for (int j = i + 1; j < new_nodes; j++) {
                int num_matches = 0;
                int new_node = my_nodes.get(j);
                for (int k = 0; k < num_clusterings; k++) {
                    if (cluster_labels.get(k).get(cur_node).equals(cluster_labels.get(k).get(new_node))) {
                        num_matches += 1;
                    }
                }                
                if (num_matches > 0) {
                    double cur_prob = ((double) num_matches) / num_clusterings;
                    Double[] j_to_i = new Double[2]; // {i, cur_prob};
                    Double[] i_to_j = new Double[2];
                    i_to_j[0] = (double) j; i_to_j[1] = cur_prob;
                    j_to_i[0] = (double) i; j_to_i[1] = cur_prob;

                    cur_prob_matrix.get(i).add(i_to_j);
                    cur_prob_matrix.get(j).add(j_to_i);

                }
            }
        }

        // System.out.println("potential edges: " + (num_nodes * (num_nodes - 1) / 2));
        // System.out.println("pos edges: " + pos_edges);
        
        // Pair<double[][], ArrayList<ArrayList<Integer>>> my_pair = new Pair();

        Helper myHelper = new Helper(cur_prob_matrix, index_mapping);
        return myHelper;
    }

    public static ArrayList<ArrayList<Integer>> relabelClustering(ArrayList<ArrayList<Integer>> clustering, ArrayList<ArrayList<Integer>> index_mapping) {
        // unmerge nodes for purposes of computing edit distance 
        ArrayList<ArrayList<Integer>> new_clustering = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < clustering.size(); i++) {
            ArrayList<Integer> cur_cluster = clustering.get(i);
            ArrayList<Integer> new_cluster = new ArrayList<Integer>();
            for (int j = 0; j < cur_cluster.size(); j++) {
                int cur_index = cur_cluster.get(j);
                ArrayList<Integer> bucket = index_mapping.get(cur_index);
                for (int k = 0; k < bucket.size(); k++)
                    new_cluster.add(bucket.get(k));

            }
            new_clustering.add(new_cluster);
        }

        return new_clustering;

    }

    public static long quickEditDist(ArrayList<ArrayList<Integer>> clustering1, ArrayList<ArrayList<Integer>> clustering2 ){
    
        int num_nodes = 0;
        long c1_energy = 0;
        for (int i = 0; i < clustering1.size(); i++) {
            long cur_size = (long) clustering1.get(i).size();
            c1_energy += cur_size * (cur_size - 1) / 2;
        }

        long c2_energy = 0;
        for (int i = 0; i < clustering2.size(); i++) {
            long cur_size = (long) clustering2.get(i).size();
            c2_energy += cur_size * (cur_size - 1) / 2;
        }

        HashMap<Integer, Integer> c1_labels = new HashMap<Integer, Integer>();
        for (int i = 0; i < clustering1.size(); i++) {
            ArrayList<Integer> cur_cluster = clustering1.get(i);
            for (int j = 0; j < cur_cluster.size(); j++) {
                c1_labels.put(cur_cluster.get(j), i);
                num_nodes +=1; // COMPUTE THIS ONCE
            }
        }

        HashMap<Integer, Integer> c2_labels = new HashMap<Integer, Integer>();
        for (int i = 0; i < clustering2.size(); i++) {
            ArrayList<Integer> cur_cluster = clustering2.get(i);
            for (int j = 0; j < cur_cluster.size(); j++) {
                c2_labels.put(cur_cluster.get(j), i);
            }
        }
        
        long union_energy = 0; 
        /*	
        String cur_label = "c" + c1_labels.get(0) + "c" + c2_labels.get(0);
        long cur_size = 1;
        for (int i = 1; i < num_nodes; i++) {
            
            String new_label = "c" + c1_labels.get(i) + "c" + c2_labels.get(i);
            if (new_label.equals(cur_label)) {
                cur_size += 1;
            } else {
                union_energy += cur_size * (cur_size - 1) / 2;
                cur_size = 1;
                cur_label = new_label;
            }               

        }
        union_energy += cur_size * (cur_size - 1) / 2;
        */

	long default_val = 0;

        HashMap<String, Long> union_labels = new HashMap<String, Long>();
        for (int i = 0; i < num_nodes; i++) {

            String cur_label = "c" + c1_labels.get(i) + "c" + c2_labels.get(i);
            union_labels.put(cur_label, union_labels.getOrDefault(cur_label, default_val) + 1);

        }

        long cur_size = 0;
        for (String my_label : union_labels.keySet()) {
            cur_size = union_labels.get(my_label);
            union_energy += cur_size * (cur_size - 1) / 2;
        }


        long dist = c1_energy + c2_energy - (2*union_energy);
        return dist;
    }

    public static long getEditDist(ArrayList<ArrayList<Integer>> consensusClustering, ArrayList<ArrayList<ArrayList<Integer>>> clusterings) {
        long dist = 0;
        for (int i = 0; i < clusterings.size(); i++) {
            ArrayList<ArrayList<Integer>> cur_clustering = clusterings.get(i);
            dist += quickEditDist(consensusClustering, cur_clustering);
        }

        return dist;

    }

    public static long getSubsetEditDist(ArrayList<ArrayList<Integer>> consensusClustering, ArrayList<Integer> nodeSet, ArrayList<HashMap<Integer, Integer>> cluster_labels) {
        long dist = 0;
        for (int i = 0; i < cluster_labels.size(); i++) {
            HashMap<Integer, Integer> cur_clustering = cluster_labels.get(i);
            dist += subsetEditDist(consensusClustering, nodeSet, cur_clustering);
            // System.out.println();
        }

        return dist;

    }

    public static long subsetEditDist(ArrayList<ArrayList<Integer>> clustering1, ArrayList<Integer> nodeSet, HashMap<Integer, Integer> clustering2){
    
        // find disagreement distance quickly given a subset of nodes to search

        int num_nodes = nodeSet.size();
        long c1_energy = 0;
        for (int i = 0; i < clustering1.size(); i++) {
            long cur_size = (long) clustering1.get(i).size();
            c1_energy += cur_size * (cur_size - 1) / 2;
        }

        HashMap<Integer, Integer> c1_labels = new HashMap<Integer, Integer>();
        for (int i = 0; i < clustering1.size(); i++) {
            ArrayList<Integer> cur_cluster = clustering1.get(i);
            // System.out.print("[");
            for (int j = 0; j < cur_cluster.size(); j++) {
                c1_labels.put(cur_cluster.get(j), i);
            }
            // System.out.print("] ");
        }
       // System.out.println();

        // HashMap<Integer, Integer> c2_labels = new HashMap<Integer, Integer>();
        HashMap<Integer, Long> c2_counts = new HashMap<Integer, Long>();
        long default_val = 0;       

        for (int i = 0; i < num_nodes; i++) {
            int cur_cluster = clustering2.get(nodeSet.get(i));
            // c2_labels.put(nodeSet.get(i), cur_cluster);            
            c2_counts.put(cur_cluster, c2_counts.getOrDefault(cur_cluster, default_val) + 1);;
        }

        long c2_energy = 0;
        for (int i : c2_counts.keySet()){
            long cur_size = (long) c2_counts.get(i);
            c2_energy += cur_size * (cur_size - 1) / 2;
        }
        
        long union_energy = 0;        
        
        HashMap<String, Long> union_labels = new HashMap<String, Long>();
        for (int i = 0; i < num_nodes; i++) {
            int cur_node = nodeSet.get(i);
            
            String cur_label = "c" + c1_labels.get(cur_node) + "c" + clustering2.get(cur_node);
            union_labels.put(cur_label, union_labels.getOrDefault(cur_label, default_val) + 1);             

        }

        long cur_size = 0;
        for (String my_label : union_labels.keySet()) {
            cur_size = union_labels.get(my_label);
            union_energy += cur_size * (cur_size - 1) / 2;
        }

        long dist = c1_energy + c2_energy - (2*union_energy);
        return dist;
    }


    // --- OLD CODE ---


    
    public ArrayList<ArrayList<Integer>> unionClustering(ArrayList<ArrayList<ArrayList<Integer>>> input_clusters, int num_nodes) {

        // step 1: compute cluster labels 
        ArrayList<HashMap<Integer, Integer>> cluster_labels = new ArrayList<HashMap<Integer, Integer>>();
        for (int i = 0; i < input_clusters.size(); i++)
            cluster_labels.add(new HashMap<Integer, Integer>());
        for (int i = 0; i < input_clusters.size(); i++) {
            ArrayList<ArrayList<Integer>> cur_clustering = input_clusters.get(i);
            HashMap<Integer, Integer> cur_map = cluster_labels.get(i);
            for (int j = 0; j < cur_clustering.size(); j++) {
                ArrayList<Integer> cur_cluster = cur_clustering.get(j);
                for (int k = 0; k < cur_cluster.size(); k++)
                    cur_map.put(cur_cluster.get(k), j);
            }
        }

        // step 2: form union clusters 
        HashMap<String, ArrayList<Integer>> union_map = new HashMap<String, ArrayList<Integer>>();
        ArrayList<String> unique_labels = new ArrayList<String>();

        for (int i = 0; i < num_nodes; i++) {
            
            String cur_label = getLabel (cluster_labels, i);
            if (union_map.containsKey(cur_label))
                union_map.get(cur_label).add(i);
            else {
                ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
                cur_cluster.add(i);
                union_map.put(cur_label, cur_cluster);
                unique_labels.add(cur_label);
            }            
        
        }

        // convert data structure 
        ArrayList<ArrayList<Integer>> union = new ArrayList<ArrayList<Integer>>();

        for (int i = 0; i < unique_labels.size(); i++) {
            String cur_label = unique_labels.get(i);
            union.add(union_map.get(cur_label));
        }

        return union; // would it be better just to output nodes to ignore?

    }


    
}
