import java.util.*;

public class DNode {

    // ECS Reverse Compare
    private static class ECSComparator implements Comparator<Double[]> {  
        public int compare (Double[] a1, Double[] a2) {  

            if (a1[0] == a2[0]) {   
                return 0;
            } else if (a1[0] > a2[0]) {  
                return -1;
            } else if (a1[0] < a2[0]) {
                return 1;
            }
            return 0;  
        }  
    }  

    // ECS Forward Compare (for local search)
    private static class ECSComparatorLS implements Comparator<Double[]> {
        public int compare (Double[] a1, Double[] a2) {

            if (a1[0] == a2[0]) {
                return 0;
            } else if (a1[0] > a2[0]) {
                return 1;
            } else if (a1[0] < a2[0]) {
                return -1;
            }
            return 0;
        }
    }


    // Deterministic Node
    private static double sum(ArrayList<Double> input_list) {
        double my_sum = 0;
        for (double i : input_list)
            my_sum += i;
        return my_sum;

    }

    private static double sum(double[] input_list) {
        double my_sum = 0;
        for (double i : input_list)
            my_sum += i;
        return my_sum;

    }

    private static Double[] find_min(ArrayList<Double[]> input_list) {

        double cur_min = input_list.get(0)[0];
        int cur_min_index = 0;

        for (int i = 0; i < input_list.size(); i++) {
            double cand = input_list.get(i)[0];
            if (cand < cur_min) {
                cur_min = cand;
                cur_min_index = i;
            }
        }
        return input_list.get(cur_min_index);

    }

    private static double intercluster_weight(int[] cluster1, ArrayList<Integer> cluster2, double[][] prob_matrix) { 
        // find average value of edges between clusters
        ArrayList<Double> cluster_weights = new ArrayList<Double>();
        for (int i : cluster1) {
            for (int j : cluster2) {
                cluster_weights.add(prob_matrix[i][j]);
            }
        }        

        return sum(cluster_weights) / ( (double) cluster_weights.size()); 
    }
    
    public static ArrayList<Double[]> node_expected_cluster_size(double[][] prob_matrix, int num_nodes) {
        // rank nodes by probability sums 
        ArrayList<Double[]> ecs_list = new ArrayList<Double[]>();

        for (int i = 0; i < num_nodes; i++) {

            double[] cur_row = prob_matrix[i];

            Double[] cur_val = {sum(cur_row), (double) i};
            ecs_list.add(cur_val);
            // System.out.println(cur_val[0]);

        }

        ecs_list.sort(new ECSComparator());
        
        return ecs_list;
    }

    public static ArrayList<Double[]> node_expected_cluster_size_LS(double[][] prob_matrix, int num_nodes) {
        // rank nodes by probability sums
        ArrayList<Double[]> ecs_list = new ArrayList<Double[]>();

        for (int i = 0; i < num_nodes; i++) {

            double[] cur_row = prob_matrix[i];

            Double[] cur_val = {sum(cur_row), (double) i};
            ecs_list.add(cur_val);
            // System.out.println(cur_val[0]);

        }

        ecs_list.sort(new ECSComparatorLS());

        return ecs_list;
    }

    
    private static double[] min_prob_edit_dist_choice(int cur_node, ArrayList<ArrayList<Integer>> cur_clustering, double[][] prob_matrix) {
        // compute clustering choice that minimizes increase of prob edit dist
        int num_clusters = cur_clustering.size(); 
        
        ArrayList<Double[]> avg_weights = new ArrayList<Double[]>();

        for (int i = 0; i < num_clusters; i++){

            int[] cluster1 = {cur_node};
            double cur_icweight = intercluster_weight(cluster1, cur_clustering.get(i), prob_matrix);
            Double[] cur_val = {cur_icweight, (double) i};
            avg_weights.add(cur_val);

        }
        
        double[] total_weights = new double[num_clusters];
        for (int i = 0; i < num_clusters; i++) {
            total_weights[i] = avg_weights.get(i)[0] * cur_clustering.get(i).size(); 
        }
        
        double all_weights = sum(total_weights);
    
        // figure out which cluster descreases total score the least, or if we should create a new cluster
        ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
        for (int i = 0; i < num_clusters; i++) {

            Double[] cur_val = {all_weights + cur_clustering.get(i).size() -  2 * total_weights[i], (double) i};
            cost_increase.add(cur_val);

        }

        Double[] last_val = {all_weights, (double) num_clusters};

        cost_increase.add(last_val); // cost of opening a new cluster
        Double[] least_cost = find_min(cost_increase);

        if (least_cost[1] < num_clusters) {
            double[] r_val = {1, least_cost[1]};
            return r_val;
        }
        else {
            double[] r_val = {0, least_cost[1]};
            return r_val; // force to open new cluster
        }
    }
    
    public static ArrayList<ArrayList<Integer>> adapted_node(double[][] prob_matrix) {
        /*
            My approach: order nodes by expected cluster size
            -- settle nodes in that order
            -- iterate over all edges between current node and the settled set
            -- only query those edges that we can't yet infer
            -- break out of the loop if we find a match for the current node
        */
        int num_nodes = prob_matrix.length; 
    
        // top-v():
        ArrayList<Double[]> ecs_list = node_expected_cluster_size(prob_matrix, num_nodes);
        ArrayList<Integer> settled_nodes = new ArrayList<Integer>();
        Integer cur_node = (int) Math.floor(ecs_list.get(0)[1]); // should be integer-valued anyway...
        settled_nodes.add(cur_node);

        // clustering is incomplete until algorithm finishes
        ArrayList<ArrayList<Integer> > cur_clustering = new ArrayList<ArrayList<Integer> >();
        ArrayList<Integer> first_cluster = new ArrayList<Integer>();
        first_cluster.add(cur_node);
     
        cur_clustering.add(first_cluster);
    
        int node_counter = 1;
        
        while (settled_nodes.size() < num_nodes) {
            cur_node = (int) Math.floor(ecs_list.get(node_counter)[1]);
            node_counter += 1;
            
            // get all cluster weights between cur_node and cur_clusters
            // cluster_weights = [[intercluster_weight([cur_node], cur_clustering[i], prob_matrix), i] for i in range(len(cur_clustering))]
            // max_cluster = max(cluster_weights)
            double[] max_cluster = min_prob_edit_dist_choice(cur_node, cur_clustering, prob_matrix);
    
            if (max_cluster[0] >= 0.5) {
                // add cur_node to the clustering
                int new_index = (int) max_cluster[1];
                cur_clustering.get(new_index).add(cur_node);
            }
            else {
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cur_node);
                cur_clustering.add(new_cluster);
            }
    
            settled_nodes.add(cur_node);
        }
            
        return cur_clustering;

    }    

    // --- 0/1 NETWORK ---

    public static ArrayList<Double[]> node_expected_cluster_size_network(ArrayList<ArrayList<Integer>> prob_matrix, int num_nodes) {
        // rank nodes by probability sums 

        ArrayList<Double[]> ecs_list = new ArrayList<Double[]>(num_nodes);
        for (int i = 0; i < num_nodes; i++) {
            Double[] cur_val = {(double) prob_matrix.get(i).size(), (double) i};
            ecs_list.add(cur_val);
        }

        ecs_list.sort(new ECSComparator());
        
        return ecs_list;
    }

    public static ArrayList<Double[]> node_expected_cluster_size_network_LS(ArrayList<ArrayList<Integer>> prob_matrix, int num_nodes) {
        // rank nodes by probability sums

        ArrayList<Double[]> ecs_list = new ArrayList<Double[]>(num_nodes);
        for (int i = 0; i < num_nodes; i++) {
            Double[] cur_val = {(double) prob_matrix.get(i).size(), (double) i};
            ecs_list.add(cur_val);
        }

        ecs_list.sort(new ECSComparatorLS());

        return ecs_list;
    }

    private static double intercluster_weight(int[] cluster1, ArrayList<Integer> cluster2, ArrayList<ArrayList<Integer>> prob_matrix, double REPS_PROB) { 
        // find average value of edges between clusters
        ArrayList<Double> cluster_weights = new ArrayList<Double>();
        int num_reps = (int) Math.ceil(REPS_PROB * cluster2.size());
        for (int i : cluster1) {
            for (int k = 0; k < num_reps; k++) {
                int j = cluster2.get(k);
                double edge_weight = 0.0;
                if (prob_matrix.get(i).contains(j))
                    edge_weight = 1.0;
                cluster_weights.add(edge_weight);
            }
        }      

        return sum(cluster_weights) / ( (double) num_reps); 

        // return sum(cluster_weights) / ( (double) cluster_weights.size()); 
    }
    
    private static double[] min_prob_edit_dist_choice(int cur_node, ArrayList<ArrayList<Integer>> cur_clustering, HashMap<Integer, Integer> settled_nodes, ArrayList<ArrayList<Integer>> prob_matrix, double REPS_PROB) {
        // compute clustering choice that minimizes increase of prob edit dist
        int num_clusters = cur_clustering.size(); 
        ArrayList<Integer> cur_edges = prob_matrix.get(cur_node);
        ArrayList<Double[]> avg_weights = new ArrayList<Double[]>();

        if (settled_nodes.size() < cur_edges.size()) { // DO IT THE OLD WAY

        for (int i = 0; i < num_clusters; i++){

            int[] cluster1 = {cur_node};
            double cur_icweight = intercluster_weight(cluster1, cur_clustering.get(i), prob_matrix, REPS_PROB);
            Double[] cur_val = {cur_icweight, (double) i};
            avg_weights.add(cur_val);

        }
        
        double[] total_weights = new double[num_clusters];
        for (int i = 0; i < num_clusters; i++) {
            total_weights[i] = avg_weights.get(i)[0] * cur_clustering.get(i).size(); 
        }
        
        double all_weights = sum(total_weights);
    
        // figure out which cluster descreases total score the least, or if we should create a new cluster
        ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
        for (int i = 0; i < num_clusters; i++) {

            Double[] cur_val = {all_weights + cur_clustering.get(i).size() -  2 * total_weights[i], (double) i};
            cost_increase.add(cur_val);

        }

        Double[] last_val = {all_weights, (double) num_clusters};

        cost_increase.add(last_val); // cost of opening a new cluster
        Double[] least_cost = find_min(cost_increase);

        if (least_cost[1] < num_clusters) {
            double[] r_val = {1, least_cost[1]};
            return r_val;
        }
        else {
            double[] r_val = {0, least_cost[1]};
            return r_val; // force to open new cluster
        }
        } else { // TAKE ADVANTAGE OF NEIGHBORHOOD ORACLE TO REDUCE "QUERIES" 

            // Calculate positive connections to established clusters 
            HashMap<Integer, Integer> pos_connections = new HashMap<Integer, Integer>();
            int total_connections = 0;
            for (int i = 0; i < cur_edges.size(); i++) {
                int cur_edge = cur_edges.get(i);
                if (settled_nodes.keySet().contains(cur_edge)) {
                    int cluster_label = settled_nodes.get(cur_edge);
                    if (pos_connections.keySet().contains(cluster_label)) 
                        pos_connections.put(cluster_label, pos_connections.get(cluster_label) + 1);
                    else
                        pos_connections.put(cluster_label, 1);
                    total_connections += 1;
                }
            }

            // Figure out best clustering choice 
            ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
            for (int i : pos_connections.keySet()) { // i is the cluster label now              

                Double[] cur_val = {total_connections + cur_clustering.get(i).size() -  (2.0 * pos_connections.get(i)), (double) i};
                cost_increase.add(cur_val);

            }

            Double[] last_val = {(double) total_connections, (double) num_clusters};
            cost_increase.add(last_val); // cost of opening a new cluster
            Double[] least_cost = find_min(cost_increase);

            if (least_cost[1] < num_clusters) {
                double[] r_val = {1, least_cost[1]};
                return r_val;
            }
            else {
                double[] r_val = {0, least_cost[1]};
                return r_val; // force to open new cluster
            }
        }
    }

    public static ArrayList<ArrayList<Integer>> adapted_node_network(ArrayList<ArrayList<Integer>> prob_matrix, int MAX_REPS) {
        /*
            My approach: order nodes by expected cluster size
            -- settle nodes in that order
            -- iterate over all edges between current node and the settled set
            -- only query those edges that we can't yet infer
            -- break out of the loop if we find a match for the current node
        */
        int num_nodes = prob_matrix.size();
        MAX_REPS = Math.min(MAX_REPS, num_nodes); // don't want max representatives to exceed cur cluster size
    
        // top-v():
        ArrayList<Double[]> ecs_list = node_expected_cluster_size_network(prob_matrix, num_nodes);
        // ArrayList<Integer> settled_nodes = new ArrayList<Integer>();
        HashMap<Integer, Integer> settled_nodes = new HashMap<Integer, Integer>();
        Integer cur_node = (int) Math.floor(ecs_list.get(0)[1]); // should be integer-valued anyway...
        settled_nodes.put(cur_node, 0);

        // clustering is incomplete until algorithm finishes
        ArrayList<ArrayList<Integer> > cur_clustering = new ArrayList<ArrayList<Integer> >();
        ArrayList<Integer> first_cluster = new ArrayList<Integer>();
        first_cluster.add(cur_node);
     
        cur_clustering.add(first_cluster);
    
        int node_counter = 1;
        double REPS_PROB = 1.0;
        
        while (settled_nodes.size() < num_nodes) {
            cur_node = (int) Math.floor(ecs_list.get(node_counter)[1]);
            node_counter += 1;

            if (MAX_REPS < settled_nodes.size())
                REPS_PROB = ((double) MAX_REPS) / settled_nodes.size();
            
            // get all cluster weights between cur_node and cur_clusters
            // cluster_weights = [[intercluster_weight([cur_node], cur_clustering[i], prob_matrix), i] for i in range(len(cur_clustering))]
            // max_cluster = max(cluster_weights)
            double[] max_cluster = min_prob_edit_dist_choice(cur_node, cur_clustering, settled_nodes, prob_matrix, REPS_PROB);
    
            if (max_cluster[0] >= 0.5) {
                // add cur_node to the clustering
                int new_index = (int) max_cluster[1];
                cur_clustering.get(new_index).add(cur_node);
                settled_nodes.put(cur_node, new_index);
            }
            else {
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cur_node);
                settled_nodes.put(cur_node, cur_clustering.size());
                cur_clustering.add(new_cluster);
            }
    
            // settled_nodes.add(cur_node);
        }
            
        return cur_clustering;

    }   

    // --- RANDOM NODE WITH NEIGHBORHOOD ORACLE --- // 


    private static double intercluster_weight_rn(int[] cluster1, ArrayList<Integer> cluster2, ArrayList<ArrayList<Integer>> prob_matrix) { 
        // find average value of edges between clusters
        ArrayList<Double> cluster_weights = new ArrayList<Double>();
        for (int i : cluster1) {
            for (int k = 0; k < cluster2.size(); k++) {
                int j = cluster2.get(k);
                double edge_weight = 0.0;
                if (prob_matrix.get(i).contains(j))
                    edge_weight = 1.0;
                cluster_weights.add(edge_weight);
            }
        }      

        // return sum(cluster_weights) / ( (double) cluster2.size()); 

        return sum(cluster_weights) / ( (double) cluster_weights.size()); 
    }
    
    private static double[] min_prob_edit_dist_choice_rn(int cur_node, ArrayList<ArrayList<Integer>> cur_clustering, HashMap<Integer, Integer> settled_nodes, ArrayList<ArrayList<Integer>> prob_matrix) {
        // compute clustering choice that minimizes increase of prob edit dist
        int num_clusters = cur_clustering.size();         
        ArrayList<Double[]> avg_weights = new ArrayList<Double[]>();
        ArrayList<Integer> cur_edges = prob_matrix.get(cur_node);
        if (settled_nodes.size() < cur_edges.size()) { // DO IT THE OLD WAY

        for (int i = 0; i < num_clusters; i++){

            int[] cluster1 = {cur_node};
            double cur_icweight = intercluster_weight_rn(cluster1, cur_clustering.get(i), prob_matrix);
            Double[] cur_val = {cur_icweight, (double) i};
            avg_weights.add(cur_val);

        }
        
        double[] total_weights = new double[num_clusters];
        for (int i = 0; i < num_clusters; i++) {
            total_weights[i] = avg_weights.get(i)[0] * cur_clustering.get(i).size(); 
        }
        
        double all_weights = sum(total_weights);
    
        // figure out which cluster decreases total score the least, or if we should create a new cluster
        ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
        for (int i = 0; i < num_clusters; i++) {

            Double[] cur_val = {all_weights + cur_clustering.get(i).size() -  2 * total_weights[i], (double) i};
            cost_increase.add(cur_val);

        }

        Double[] last_val = {all_weights, (double) num_clusters};

        cost_increase.add(last_val); // cost of opening a new cluster
        Double[] least_cost = find_min(cost_increase);

        if (least_cost[1] < num_clusters) {
            double[] r_val = {1, least_cost[1]};
            return r_val;
        }
        else {
            double[] r_val = {0, least_cost[1]};
            return r_val; // force to open new cluster
        }
        } else { // TAKE ADVANTAGE OF NEIGHBORHOOD ORACLE TO REDUCE "QUERIES" 

            // Calculate positive connections to established clusters 
            HashMap<Integer, Integer> pos_connections = new HashMap<Integer, Integer>();
            int total_connections = 0;
            for (int i = 0; i < cur_edges.size(); i++) {
                int cur_edge = cur_edges.get(i);
                if (settled_nodes.keySet().contains(cur_edge)) {
                    int cluster_label = settled_nodes.get(cur_edge);
                    if (pos_connections.keySet().contains(cluster_label)) 
                        pos_connections.put(cluster_label, pos_connections.get(cluster_label) + 1);
                    else
                        pos_connections.put(cluster_label, 1);
                    total_connections += 1;
                }
            }

            // Figure out best clustering choice 
            ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
            for (int i : pos_connections.keySet()) { // i is the cluster label now              

                Double[] cur_val = {total_connections + cur_clustering.get(i).size() -  (2.0 * pos_connections.get(i)), (double) i};
                cost_increase.add(cur_val);

            }

            Double[] last_val = {(double) total_connections, (double) num_clusters};
            cost_increase.add(last_val); // cost of opening a new cluster
            Double[] least_cost = find_min(cost_increase);

            if (least_cost[1] < num_clusters) {
                double[] r_val = {1, least_cost[1]};
                return r_val;
            }
            else {
                double[] r_val = {0, least_cost[1]};
                return r_val; // force to open new cluster
            }

        }
    }

    public static ArrayList<ArrayList<Integer>> random_node_network(ArrayList<ArrayList<Integer>> prob_matrix) {
        /*
            RandomNode ordering // TODO: use same order as corresponding Pivot clustering?
            Assume 0/1 graph; Make use of Neighborhood oracle provided  
        */
        int num_nodes = prob_matrix.size();
    
        // top-v():
        ArrayList<Integer> ecs_list = new ArrayList<Integer>(num_nodes);
        for (int i = 0; i < num_nodes; i++)
            ecs_list.add(i);
        Collections.shuffle(ecs_list);

        HashMap<Integer, Integer> settled_nodes = new HashMap<Integer, Integer>();
        Integer cur_node = (int) ecs_list.get(0); 
        settled_nodes.put(cur_node, 0); // first node is assigned to cluster 0 

        // clustering is incomplete until algorithm finishes
        ArrayList<ArrayList<Integer> > cur_clustering = new ArrayList<ArrayList<Integer> >();
        ArrayList<Integer> first_cluster = new ArrayList<Integer>();
        first_cluster.add(cur_node);
     
        cur_clustering.add(first_cluster);
    
        int node_counter = 1;
        
        while (settled_nodes.size() < num_nodes) {
            cur_node = (int) ecs_list.get(node_counter);
            node_counter += 1;
            
            // get all cluster weights between cur_node and cur_clusters
            // cluster_weights = [[intercluster_weight([cur_node], cur_clustering[i], prob_matrix), i] for i in range(len(cur_clustering))]
            // max_cluster = max(cluster_weights)
            double[] max_cluster = min_prob_edit_dist_choice_rn(cur_node, cur_clustering, settled_nodes, prob_matrix);
    
            if (max_cluster[0] >= 0.5) {
                // add cur_node to the clustering
                int new_index = (int) max_cluster[1];
                cur_clustering.get(new_index).add(cur_node);
                settled_nodes.put(cur_node, new_index);
            }
            else {
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cur_node);
                settled_nodes.put(cur_node, cur_clustering.size());
                cur_clustering.add(new_cluster);
            }
    
            
        }
            
        return cur_clustering;

    }  


    // LOCAL SEARCH IMPROVEMENTS

    private static double[] min_cost_choice_ls(int cur_node, ArrayList<Integer> cluster_sizes, HashMap<Integer, Integer> cluster_map, ArrayList<ArrayList<Integer>> prob_matrix) {
        // compute clustering choice that minimizes increase of prob edit dist
        int num_clusters = cluster_sizes.size();         
        ArrayList<Double[]> avg_weights = new ArrayList<Double[]>();
        ArrayList<Integer> cur_edges = prob_matrix.get(cur_node);

        // Calculate positive connections to established clusters 
        HashMap<Integer, Integer> pos_connections = new HashMap<Integer, Integer>();
        int total_connections = 0;
        for (int i = 0; i < cur_edges.size(); i++) {
            int cur_edge = cur_edges.get(i);
            int cluster_label = cluster_map.get(cur_edge);
            if (pos_connections.keySet().contains(cluster_label)) 
                pos_connections.put(cluster_label, pos_connections.get(cluster_label) + 1);
            else
                pos_connections.put(cluster_label, 1);
            total_connections += 1;

        }

        // Figure out best clustering choice 
        ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
        for (int i : pos_connections.keySet()) { // i is the cluster label now              

            Double[] cur_val = {total_connections + cluster_sizes.get(i) -  (2.0 * pos_connections.get(i)), (double) i};
            cost_increase.add(cur_val);

        }

        Double[] last_val = {(double) total_connections, (double) num_clusters};
        cost_increase.add(last_val); // cost of opening a new cluster
        Double[] least_cost = find_min(cost_increase);

        if (least_cost[1] < num_clusters) {
            double[] r_val = {1, least_cost[1]};
            return r_val;
        }
        else {
            double[] r_val = {0, least_cost[1]};
            return r_val; // force to open new cluster
        }

    }


    public static ArrayList<ArrayList<Integer>> local_search_network(ArrayList<ArrayList<Integer>> cur_clustering, ArrayList<ArrayList<Integer>> prob_matrix, boolean ECS) {
        /*
            Make one pass through given clustering, using "best of one-element moves" strategy
            Assume 0/1 graph; Make use of Neighborhood oracle provided  
        */
        int num_nodes = prob_matrix.size();
        ArrayList<Double[]> ecs_list;

        if (ECS)
            ecs_list = node_expected_cluster_size_network_LS(prob_matrix, num_nodes);
        else {    
            // Get random node order
            ecs_list = new ArrayList<Double[]>(num_nodes);
            for (int i = 0; i < num_nodes; i++) {
                Double[] cur_entry = {0.0, (double) i};
                ecs_list.add(cur_entry);
            }
            Collections.shuffle(ecs_list);
        }

        // get map of <node, cluster_label> pairs
        HashMap<Integer, Integer> cluster_map = Helper.get_clustering_map(cur_clustering);
        
        ArrayList<Integer> cluster_sizes = new ArrayList<Integer>();
        for (int i = 0; i < cur_clustering.size(); i++) 
            cluster_sizes.add(cur_clustering.get(i).size());

  
        long cur_score = Helper.quick_edit_dist(cur_clustering, prob_matrix);
        HashMap<Integer, Integer> cluster_relabel = new HashMap<Integer, Integer>();
        ArrayList<ArrayList<Integer>> new_clustering = new ArrayList<ArrayList<Integer>>();  

        while (true) {
       
            for (int i = 0; i < num_nodes; i++) {
                Integer cur_node = (int) Math.floor(ecs_list.get(i)[1]);            
                double[] max_cluster = min_cost_choice_ls(cur_node, cluster_sizes, cluster_map, prob_matrix);

                int old_index = cluster_map.get(cur_node);

                // Improvements
                if (max_cluster[0] >= 0.5) {    
                    // add to existing cluster
                    int new_index = (int) max_cluster[1];
                    if (!cluster_relabel.keySet().contains(new_index)) {
                        ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                        new_cluster.add(cur_node);
                        cluster_relabel.put(new_index, new_clustering.size());
                        new_clustering.add(new_cluster);
                    } else {
                        new_clustering.get(cluster_relabel.get(new_index)).add(cur_node);
                    }

                    cluster_map.put(cur_node, new_index);
                    cluster_sizes.set(old_index, cluster_sizes.get(old_index) - 1);
                    cluster_sizes.set(new_index, cluster_sizes.get(new_index) + 1);

                }
                else {
                    // create new cluster 
                    ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                    new_cluster.add(cur_node);
                    cluster_map.put(cur_node, cluster_sizes.size());
                    cluster_relabel.put(cluster_sizes.size(), new_clustering.size());
                    cluster_sizes.add(new_cluster.size());
                    new_clustering.add(new_cluster);
                }
            }

            long new_score = Helper.quick_edit_dist(new_clustering, prob_matrix);
            if (new_score < cur_score) {
                cur_score = new_score; 
                cluster_relabel = new HashMap<Integer, Integer>();
                new_clustering = new ArrayList<ArrayList<Integer>>(); 
            } else {
                break;
            }
        }

            
        return new_clustering;

    }    

    public static ArrayList<ArrayList<Integer>> local_search_network_score(ArrayList<ArrayList<Integer>> cur_clustering, ArrayList<ArrayList<Integer>> prob_matrix, boolean ECS, long target_score) {
        /*
            Make one pass through given clustering, using "best of one-element moves" strategy
            Assume 0/1 graph; Make use of Neighborhood oracle provided  
        */
        int num_nodes = prob_matrix.size();
        ArrayList<Double[]> ecs_list;

        if (ECS)
            ecs_list = node_expected_cluster_size_network_LS(prob_matrix, num_nodes);
        else {    
            // Get random node order
            ecs_list = new ArrayList<Double[]>(num_nodes);
            for (int i = 0; i < num_nodes; i++) {
                Double[] cur_entry = {0.0, (double) i};
                ecs_list.add(cur_entry);
            }
            Collections.shuffle(ecs_list);
        }

        // get map of <node, cluster_label> pairs
        HashMap<Integer, Integer> cluster_map = Helper.get_clustering_map(cur_clustering);
        
        ArrayList<Integer> cluster_sizes = new ArrayList<Integer>();
        for (int i = 0; i < cur_clustering.size(); i++) 
            cluster_sizes.add(cur_clustering.get(i).size());

  
        // long cur_score = Helper.quick_edit_dist(cur_clustering, prob_matrix);
        HashMap<Integer, Integer> cluster_relabel = new HashMap<Integer, Integer>();
        ArrayList<ArrayList<Integer>> new_clustering = new ArrayList<ArrayList<Integer>>();  

        while (true) {
       
            for (int i = 0; i < num_nodes; i++) {
                Integer cur_node = (int) Math.floor(ecs_list.get(i)[1]);            
                double[] max_cluster = min_cost_choice_ls(cur_node, cluster_sizes, cluster_map, prob_matrix);

                int old_index = cluster_map.get(cur_node);

                // Improvements
                if (max_cluster[0] >= 0.5) {    
                    // add to existing cluster
                    int new_index = (int) max_cluster[1];
                    if (!cluster_relabel.keySet().contains(new_index)) {
                        ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                        new_cluster.add(cur_node);
                        cluster_relabel.put(new_index, new_clustering.size());
                        new_clustering.add(new_cluster);
                    } else {
                        new_clustering.get(cluster_relabel.get(new_index)).add(cur_node);
                    }

                    cluster_map.put(cur_node, new_index);
                    cluster_sizes.set(old_index, cluster_sizes.get(old_index) - 1);
                    cluster_sizes.set(new_index, cluster_sizes.get(new_index) + 1);

                }
                else {
                    // create new cluster 
                    ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                    new_cluster.add(cur_node);
                    cluster_map.put(cur_node, cluster_sizes.size());
                    cluster_relabel.put(cluster_sizes.size(), new_clustering.size());
                    cluster_sizes.add(new_cluster.size());
                    new_clustering.add(new_cluster);
                }
            }

            long new_score = Helper.quick_edit_dist(new_clustering, prob_matrix);
            if (new_score > target_score) {
                // cur_score = new_score; 
                cluster_relabel = new HashMap<Integer, Integer>();
                new_clustering = new ArrayList<ArrayList<Integer>>(); 
            } else {
                break;
            }
        }

            
        return new_clustering;

    } 

    public static ArrayList<ArrayList<Integer>> local_search_network_timed(ArrayList<ArrayList<Integer>> cur_clustering, ArrayList<ArrayList<Integer>> prob_matrix, boolean ECS, long startTime, long limit) {
        /*
            Make one pass through given clustering, using "best of one-element moves" strategy
            Assume 0/1 graph; Make use of Neighborhood oracle provided  
        */
        int num_nodes = prob_matrix.size();
        ArrayList<Double[]> ecs_list;
        long currentTime = 0; 

        if (ECS)
            ecs_list = node_expected_cluster_size_network_LS(prob_matrix, num_nodes);
        else {    
            // Get random node order
            ecs_list = new ArrayList<Double[]>(num_nodes);
            for (int i = 0; i < num_nodes; i++) {
                Double[] cur_entry = {0.0, (double) i};
                ecs_list.add(cur_entry);
            }
            Collections.shuffle(ecs_list);
        }

        // get map of <node, cluster_label> pairs
        HashMap<Integer, Integer> cluster_map = Helper.get_clustering_map(cur_clustering);
        
        ArrayList<Integer> cluster_sizes = new ArrayList<Integer>();
        for (int i = 0; i < cur_clustering.size(); i++) 
            cluster_sizes.add(cur_clustering.get(i).size());

  
        long cur_score = Helper.quick_edit_dist(cur_clustering, prob_matrix);
        HashMap<Integer, Integer> cluster_relabel = new HashMap<Integer, Integer>();
        ArrayList<ArrayList<Integer>> new_clustering = new ArrayList<ArrayList<Integer>>();
        
        boolean flag = true;

        while (flag) {
       
            for (int i = 0; i < num_nodes; i++) {

                currentTime = System.currentTimeMillis() - startTime;
                if (currentTime > limit) {
                    flag = false;
                }

                Integer cur_node = (int) Math.floor(ecs_list.get(i)[1]);
                int old_index = cluster_map.get(cur_node);

                double[] max_cluster = {1.0, (double) old_index};

                if (flag) 
                    max_cluster = min_cost_choice_ls(cur_node, cluster_sizes, cluster_map, prob_matrix);


                // Improvements
                if (max_cluster[0] >= 0.5) {    
                    // add to existing cluster
                    int new_index = (int) max_cluster[1];
                    if (!cluster_relabel.keySet().contains(new_index)) {
                        ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                        new_cluster.add(cur_node);
                        cluster_relabel.put(new_index, new_clustering.size());
                        new_clustering.add(new_cluster);
                    } else {
                        new_clustering.get(cluster_relabel.get(new_index)).add(cur_node);
                    }

                    cluster_map.put(cur_node, new_index);
                    cluster_sizes.set(old_index, cluster_sizes.get(old_index) - 1);
                    cluster_sizes.set(new_index, cluster_sizes.get(new_index) + 1);

                }
                else {
                    // create new cluster 
                    ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                    new_cluster.add(cur_node);
                    cluster_map.put(cur_node, cluster_sizes.size());
                    cluster_relabel.put(cluster_sizes.size(), new_clustering.size());
                    cluster_sizes.add(new_cluster.size());
                    new_clustering.add(new_cluster);
                }
                
            }
         
            if (flag) {
                // prep for next round
                cluster_relabel = new HashMap<Integer, Integer>();
                new_clustering = new ArrayList<ArrayList<Integer>>();
            }
        }

            
        return new_clustering;

    } 


    // --- CONSENSUS ON THE FLY ---

    private static double getSimilarity(ArrayList<HashMap<Integer, Integer>> cluster_labels, int node1, int node2, int MAX_ATTRIBUTES) {

        // Random r = new Random();

        //for (int ii = 0; ii < MAX_ATTRIBUTES; ii++) // independent samples each time
         //   Collections.swap(cluster_labels, ii , ii + r.nextInt(cluster_labels.size() - ii));
        
        // Consensus.attributeShuffle(cluster_labels, MAX_ATTRIBUTES);

        int num_matches = 0;
        for (int k = 0; k < MAX_ATTRIBUTES; k++) {

            if (cluster_labels.get(k).get(node1).equals(cluster_labels.get(k).get(node2))) 
                num_matches += 1;
        }
        double sim = ((double) num_matches) / MAX_ATTRIBUTES;

        return sim;

    }


    public static ArrayList<Double[]> node_expected_cluster_size_consensus(ArrayList<HashMap<Integer, Integer>> cluster_labels, ArrayList<Integer> cluster_to_fix, int MAX_ATTRIBUTES, int MAX_REPS, boolean useRandom) {
        // rank nodes by probability sums 
        int num_nodes = cluster_to_fix.size();

        ArrayList<Double[]> ecs_list = new ArrayList<Double[]>();
        for (int i = 0; i < num_nodes; i++) {
            Double[] cur_val = {0.0, (double) i};
            ecs_list.add(cur_val);
        }

        // EXPERIMENT: when to use random order instead of ECS 
	if (useRandom == true)
	    return ecs_list;

        
        for (int i = 0; i < MAX_REPS; i++) {
            for (int j = i+1; j < MAX_REPS; j++) {                
                double sim = getSimilarity(cluster_labels, cluster_to_fix.get(i), cluster_to_fix.get(j), MAX_ATTRIBUTES);
                ecs_list.get(i)[0] += sim;
                ecs_list.get(j)[0] += sim;                
            }
        }

        for (int i = MAX_REPS; i < num_nodes; i++) {
            for (int j = 0; j < MAX_REPS; j++) {
                double sim = getSimilarity(cluster_labels, cluster_to_fix.get(i), cluster_to_fix.get(j), MAX_ATTRIBUTES);
                ecs_list.get(i)[0] += sim;
            }

            /*
            int j = 0;
            int count = 0;
            while (count < MAX_REPS - 1) { // -1 to balance i == j case from before
                if (i != j) {
                    double sim = getSimilarity(cluster_labels, cluster_to_fix.get(i), cluster_to_fix.get(j));
                    ecs_list.get(i)[0] += sim;
                    // ecs_list.get(j)[0] += sim;
                    count++;
                }
                j++;               
            } */
        }

        ecs_list.sort(new ECSComparator());
        
        return ecs_list;
        
    }

    private static double intercluster_weight_consensus(int[] cluster1, ArrayList<Integer> cluster2, ArrayList<HashMap<Integer, Integer>> cluster_labels, int MAX_ATTRIBUTES, double REPS_PROB) { 
        // find average value of edges between clusters
        ArrayList<Double> cluster_weights = new ArrayList<Double>();
        int num_reps = (int) Math.ceil(REPS_PROB * cluster2.size()); // Math.min(MAX_REPS, cluster2.size()); 

        for (int i : cluster1) {
            for (int k = 0 ; k < num_reps; k++) {
                int j = cluster2.get(k);
                cluster_weights.add(getSimilarity(cluster_labels, i, j, MAX_ATTRIBUTES));
            }
        }        

        return sum(cluster_weights) / ( (double) num_reps); 
    }

    private static double[] min_prob_edit_dist_choice_consensus(int cur_node, ArrayList<ArrayList<Integer>> cur_clustering, ArrayList<HashMap<Integer, Integer>> cluster_labels, int MAX_ATTRIBUTES, double REPS_PROB) {
        // compute clustering choice that minimizes increase of prob edit dist
        int num_clusters = cur_clustering.size(); 
        
        ArrayList<Double[]> avg_weights = new ArrayList<Double[]>();

        for (int i = 0; i < num_clusters; i++){

            int[] cluster1 = {cur_node};
            double cur_icweight = intercluster_weight_consensus(cluster1, cur_clustering.get(i), cluster_labels, MAX_ATTRIBUTES, REPS_PROB);
            Double[] cur_val = {cur_icweight, (double) i};
            avg_weights.add(cur_val);

        }
        
        double[] total_weights = new double[num_clusters];
        for (int i = 0; i < num_clusters; i++) {
            total_weights[i] = avg_weights.get(i)[0] * cur_clustering.get(i).size(); 
        }
        
        double all_weights = sum(total_weights);
    
        // figure out which cluster descreases total score the least, or if we should create a new cluster
        ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
        for (int i = 0; i < num_clusters; i++) {

            Double[] cur_val = {all_weights + cur_clustering.get(i).size() -  2 * total_weights[i], (double) i};
            cost_increase.add(cur_val);

        }

        Double[] last_val = {all_weights, (double) num_clusters};

        cost_increase.add(last_val); // cost of opening a new cluster
        Double[] least_cost = find_min(cost_increase);

        if (least_cost[1] < num_clusters) {
            double[] r_val = {1, least_cost[1]};
            return r_val;
        }
        else {
            double[] r_val = {0, least_cost[1]};
            return r_val; // force to open new cluster
        }
    }

    public static ArrayList<ArrayList<Integer>> adapted_node_consensus(ArrayList<HashMap<Integer, Integer>> cluster_labels, ArrayList<Integer> cluster_to_fix, int MAX_ATTRIBUTES, int MAX_REPS, boolean useRandom) {
        /*
            My approach: order nodes by expected cluster size
            -- settle nodes in that order
            -- iterate over all edges between current node and the settled set
            -- only query those edges that we can't yet infer
            -- break out of the loop if we find a match for the current node
        */
        int num_nodes = cluster_to_fix.size();
        MAX_REPS = Math.min(MAX_REPS, num_nodes); // don't want max representatives to exceed cur cluster size
        
    
        // top-v():
        ArrayList<Double[]> ecs_list = node_expected_cluster_size_consensus(cluster_labels, cluster_to_fix, MAX_ATTRIBUTES, MAX_REPS, useRandom);
        ArrayList<Integer> settled_nodes = new ArrayList<Integer>();
        Integer cur_node = (int) Math.floor(ecs_list.get(0)[1]); // should be integer-valued anyway...
        settled_nodes.add(cur_node);

        // clustering is incomplete until algorithm finishes
        ArrayList<ArrayList<Integer> > cur_clustering = new ArrayList<ArrayList<Integer> >();
        ArrayList<Integer> first_cluster = new ArrayList<Integer>();
        first_cluster.add(cur_node);
     
        cur_clustering.add(first_cluster);
    
        int node_counter = 1;
        double REPS_PROB = 1.0;
        
        while (settled_nodes.size() < num_nodes) {
            cur_node = (int) Math.floor(ecs_list.get(node_counter)[1]);
            node_counter += 1;
            if (MAX_REPS < settled_nodes.size())
                REPS_PROB = ((double) MAX_REPS) / settled_nodes.size();
            double[] max_cluster = min_prob_edit_dist_choice_consensus(cluster_to_fix.get(cur_node), cur_clustering, cluster_labels, MAX_ATTRIBUTES, REPS_PROB); //MAX_REPS); 
    
            if (max_cluster[0] >= 0.5) {
                // add cur_node to the clustering
                int new_index = (int) max_cluster[1];
                cur_clustering.get(new_index).add(cluster_to_fix.get(cur_node));
            }
            else {
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cluster_to_fix.get(cur_node));
                cur_clustering.add(new_cluster);
            }
    
            settled_nodes.add(cur_node);
        }
            
        return cur_clustering;

    }

    // --- LOCAL SEARCH CONSENSUS ---

    public static ArrayList<Double[]> node_expected_cluster_size_consensus_ls(ArrayList<HashMap<Integer, Integer>> cluster_labels, ArrayList<Integer> cluster_to_fix, int MAX_ATTRIBUTES, boolean useRandom) {
        // rank nodes by probability sums 
        int num_nodes = cluster_to_fix.size();

        ArrayList<Double[]> ecs_list = new ArrayList<Double[]>();
        for (int i = num_nodes - 1; i >= 0; i--) { // add in reverse order so that pivot node is last
            Double[] cur_val = {0.0, (double) i};
            ecs_list.add(cur_val);
        }

        // EXPERIMENT: when to use random order instead of ECS 
	if (useRandom == true) {
            Collections.shuffle(ecs_list); // going back to random order instead of reverse order
	    return ecs_list;    
        }    
        
        for (int i = 0; i < num_nodes; i++) {
            for (int j = i+1; j < num_nodes; j++) {                
                double sim = getSimilarity(cluster_labels, cluster_to_fix.get(i), cluster_to_fix.get(j), MAX_ATTRIBUTES);
                ecs_list.get(i)[0] += sim;
                ecs_list.get(j)[0] += sim;                
            }
        }

        ecs_list.sort(new ECSComparatorLS());
        
        return ecs_list;
        
    }

    private static double intercluster_weight_consensus_ls(int[] cluster1, ArrayList<Integer> cluster2, ArrayList<HashMap<Integer, Integer>> cluster_labels, int MAX_ATTRIBUTES) { 
        // find average value of edges between clusters
        // ArrayList<Double> cluster_weights = new ArrayList<Double>();
        int num_reps = cluster2.size(); // Math.min(MAX_REPS, cluster2.size()); 
        double weight_total = 0;

        for (int i : cluster1) {
            for (int k = 0 ; k < num_reps; k++) {
                int j = cluster2.get(k);
                if (i != j)
                    weight_total += getSimilarity(cluster_labels, i, j, MAX_ATTRIBUTES);
            }
        }        

        return weight_total; 
    }

    private static double[] min_prob_edit_dist_choice_consensus_ls(int cur_node, int cur_index, ArrayList<ArrayList<Integer>> cur_clustering, ArrayList<HashMap<Integer, Integer>> cluster_labels, int MAX_ATTRIBUTES) {
        // compute clustering choice that minimizes increase of prob edit dist
        int num_clusters = cur_clustering.size(); 
        
        // ArrayList<Double[]> avg_weights = new ArrayList<Double[]>();
        double[] total_weights = new double[num_clusters];
        for (int i = 0; i < num_clusters; i++){

            int[] cluster1 = {cur_node};
            total_weights[i] = intercluster_weight_consensus_ls(cluster1, cur_clustering.get(i), cluster_labels, MAX_ATTRIBUTES);
            // Double[] cur_val = {cur_icweight, (double) i};
            // avg_weights.add(cur_val);

        }        
       
        // for (int i = 0; i < num_clusters; i++) {
        //     total_weights[i] = avg_weights.get(i)[0] * cur_clustering.get(i).size(); 
        // }
        
        double all_weights = sum(total_weights);
    
        // figure out which cluster descreases total score the least, or if we should create a new cluster
        ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
        for (int i = 0; i < num_clusters; i++) {

            if (i != cur_index) {
                Double[] cur_val = {all_weights + cur_clustering.get(i).size() -  (2 * total_weights[i]), (double) i};
                cost_increase.add(cur_val);
            } else {
                Double[] cur_val = {all_weights + (cur_clustering.get(i).size() - 1) -  (2 * total_weights[i]), (double) i};
                cost_increase.add(cur_val);
            }

        }

        Double[] last_val = {all_weights, (double) num_clusters};

        cost_increase.add(last_val); // cost of opening a new cluster
        Double[] least_cost = find_min(cost_increase);

        if (least_cost[1] < num_clusters) {
            double[] r_val = {1, least_cost[1]};
            return r_val;
        }
        else {
            double[] r_val = {0, least_cost[1]};
            return r_val; // force to open new cluster
        }
    }


    public static ArrayList<ArrayList<Integer>> local_search_consensus(ArrayList<HashMap<Integer, Integer>> cluster_labels, ArrayList<Integer> cluster_to_fix, int MAX_ATTRIBUTES, boolean useRandom) {
        /*
            My approach: order nodes by expected cluster size
            -- settle nodes in that order
            -- iterate over all edges between current node and the settled set
            -- only query those edges that we can't yet infer
            -- break out of the loop if we find a match for the current node
        */
        int num_nodes = cluster_to_fix.size();
    
        // top-v():
        ArrayList<Double[]> ecs_list = node_expected_cluster_size_consensus_ls(cluster_labels, cluster_to_fix, MAX_ATTRIBUTES, useRandom);
        ArrayList<Integer> settled_nodes = new ArrayList<Integer>();

        // clustering is incomplete until algorithm finishes
        ArrayList<ArrayList<Integer> > cur_clustering = new ArrayList<ArrayList<Integer> >();
        ArrayList<Integer> cluster_copy = new ArrayList<Integer>();
        cluster_copy.addAll(cluster_to_fix); 
        cur_clustering.add(cluster_copy);
    
        int node_counter = 0;
        
        while (settled_nodes.size() < num_nodes) {
            int cur_node = (int) Math.floor(ecs_list.get(node_counter)[1]);
            node_counter += 1;
            double[] max_cluster = min_prob_edit_dist_choice_consensus_ls(cluster_to_fix.get(cur_node), 0, cur_clustering, cluster_labels, MAX_ATTRIBUTES); //MAX_REPS); 
    
            if (max_cluster[0] >= 0.5) {
                // add cur_node to the clustering
                int new_index = (int) max_cluster[1];
                if (new_index != 0) { // moves out of original cluster
                    cur_clustering.get(0).remove(cluster_to_fix.get(cur_node));
                    cur_clustering.get(new_index).add(cluster_to_fix.get(cur_node));
                }
            }
            else {
                cur_clustering.get(0).remove(cluster_to_fix.get(cur_node));
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cluster_to_fix.get(cur_node));
                cur_clustering.add(new_cluster);
            }    
            settled_nodes.add(cur_node);
        }
            
        return cur_clustering;

    }

    // FULL LOCAL SEARCH
    public static ArrayList<ArrayList<Integer>> full_local_search_consensus(int num_nodes, ArrayList<HashMap<Integer, Integer>> cluster_labels, ArrayList<ArrayList<Integer>> clustering_to_fix, int MAX_ATTRIBUTES) {
        /*
            My approach: order nodes by expected cluster size
            -- settle nodes in that order
            -- iterate over all edges between current node and the settled set
            -- only query those edges that we can't yet infer
            -- break out of the loop if we find a match for the current node
        */
        // int num_nodes = cluster_labels.size();
    
        // top-v():
        ArrayList<Integer> ecs_list = new ArrayList<Integer>(num_nodes);
        for (int i = 0; i < clustering_to_fix.size(); i++) {
            ArrayList<Integer> cur_cluster = clustering_to_fix.get(i);
            for (int j = 0; j < cur_cluster.size(); j++) {
                ecs_list.add(cur_cluster.get(j));               
            }
        }
        Collections.shuffle(ecs_list); // random order

        HashMap<Integer, Integer> cluster_map = Helper.get_clustering_map(clustering_to_fix);
        ArrayList<Integer> settled_nodes = new ArrayList<Integer>();

        // clustering is incomplete until algorithm finishes
        // ArrayList<ArrayList<Integer> > cur_clustering = new ArrayList<ArrayList<Integer> >();
        // ArrayList<Integer> cluster_copy = new ArrayList<Integer>();
        // cluster_copy.addAll(cluster_to_fix); 
        // cur_clustering.add(cluster_copy);
    
        int node_counter = 0;
        
        while (settled_nodes.size() < num_nodes) {
            Integer cur_node = ecs_list.get(node_counter);
            // System.out.println(cur_node);
            node_counter += 1;
            double[] max_cluster = min_prob_edit_dist_choice_consensus_ls(cur_node, cluster_map.get(cur_node), clustering_to_fix, cluster_labels, MAX_ATTRIBUTES); //MAX_REPS); 
    
            if (max_cluster[0] >= 0.5) {
                // add cur_node to the clustering
                int new_index = (int) max_cluster[1];
                if (new_index != cluster_map.get(cur_node)) { // moves out of original cluster
                    clustering_to_fix.get(cluster_map.get(cur_node)).remove(cur_node);
                    clustering_to_fix.get(new_index).add(cur_node);
                    cluster_map.put(cur_node, new_index);
                }
            }
            else {
                
                clustering_to_fix.get(cluster_map.get(cur_node)).remove(cur_node);
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cur_node);
                cluster_map.put(cur_node, clustering_to_fix.size());
                clustering_to_fix.add(new_cluster);
            }    
            settled_nodes.add(cur_node);
        }
            
        return clustering_to_fix;

    } 

    
}
