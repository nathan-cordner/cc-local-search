import java.util.*;
import java.io.*;

public class PKwik {

    // LARGE NETWORK: READ AS ADJACENCY LIST
    public static ArrayList<ArrayList<Integer>> pKwikClustering(ArrayList<ArrayList<Integer>> prob_matrix) {
        /*
            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 
    
        // int loop1 = 0;
        // int loop2 = 0;
        // long start = 0;

        int cur_index = 0;

        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);
            // available_nodes.remove(cur_index);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);
      

            // start = System.currentTimeMillis();
            for (int neighbor : prob_matrix.get(cur_node)) {
                if (available_nodes.contains(neighbor)) {
                    cur_cluster.add(neighbor);
                }
            }
            // loop1 += System.currentTimeMillis() - start;

            // start = System.currentTimeMillis();
            available_nodes.removeAll(cur_cluster);
            // for (int i = 1; i < cur_cluster.size(); i++) 
            //     available_nodes.remove(cur_cluster.get(i));
            // loop2 += System.currentTimeMillis() - start;
            clusters.add(cur_cluster);            
        }
        // System.out.println("Loop 1 finished in: " + loop1 / 1000.0 + " s");
        // System.out.println("Loop 2 finished in: " + loop2 / 1000.0 + " s");
            
        return clusters;
    }

    // INPUT: probability matrix
    public static ArrayList<ArrayList<Integer>> pKwikClusteringProb(ArrayList<ArrayList<Double[]>> prob_matrix) {
        /*
            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 

        int cur_index = 0;

        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);

            for (Double[] values : prob_matrix.get(cur_node)) {
                double temp = values[0];
                int neighbor = (int) temp;
                double prob = values[1];

                if (prob >= 0.5 && available_nodes.contains(neighbor))
                    cur_cluster.add(neighbor);
            }

            available_nodes.removeAll(cur_cluster);
            clusters.add(cur_cluster);            
        }
            
        return clusters;
    }

    public static ArrayList<ArrayList<Integer>> pKwikClusteringConsensus(ArrayList<HashMap<Integer, Integer>> cluster_labels, int num_nodes) {
        /*
            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node
        */

        double threshold = cluster_labels.size() / 2.0; // number of attributes
        
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 

        int cur_index = 0;

        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);

            Iterator<Integer> myIterator = available_nodes.iterator();
            while(myIterator.hasNext()) {
                int other_node = myIterator.next();
                if (cur_node != other_node) {
                    // calculate probability of matching cur_pivot
                    int num_matches = 0;
                    for (int k = 0; k < cluster_labels.size(); k++) {
                        // System.out.println("cur node: " + cur_node + ", other node: " + other_node);
                        try {
                            cluster_labels.get(k).get(cur_node).equals(cluster_labels.get(k).get(other_node));
                        } catch (NullPointerException e) {
                            System.out.println(k);
                            System.out.println("Total nodes: " + num_nodes);
                            System.out.println("cur node: " + cur_node + ", other node: " + other_node);
                            System.out.println(cluster_labels.get(k).get(cur_node));
                            System.out.println(cluster_labels.get(k).get(other_node));
                        }

                        if (cluster_labels.get(k).get(cur_node).equals(cluster_labels.get(k).get(other_node))) {
                            num_matches += 1;
                        }
                        
                        if (num_matches >= threshold) {
                            cur_cluster.add(other_node);
                            break;
                        }
                    }


                }
            }           


            available_nodes.removeAll(cur_cluster);
            clusters.add(cur_cluster);            
        }
        return clusters;
    }

    public static ArrayList<ArrayList<Integer>> permutationPKwikConsensus(ArrayList<HashMap<Integer, Integer>> cluster_labels, int num_nodes, int MAX_PIVOTS) {
        /*
            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node
        */
        
        double threshold = cluster_labels.size() / 2.0; // MAX_PIVOTS / 2.0; //cluster_labels.size() / 2.0; // number of attributes
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i); 

        Collections.shuffle(permutation);

	// EXPERIMENT: now using MAX_PIVOTS to specify max number of attributes
	int MAX_R = cluster_labels.size(); // MAX_PIVOTS;
	MAX_PIVOTS = num_nodes;
	// Random r = new Random();


        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 
        ArrayList<Integer> pivot_nodes = new ArrayList<Integer>();
        pivot_nodes.add(permutation.get(0));

        clusters.add(new ArrayList<Integer>());
        clusters.get(0).add(permutation.get(0));
        for (int i = 1; i < num_nodes; i++) {   

            int cur_node = permutation.get(i);
            boolean added = false;
            
            // Collections.shuffle(cluster_labels); // EXPERIMENT (CAN MAKE THIS MORE EFFICIENT LATER)
            // for (int ii = 0; ii < MAX_R; ii++)
            //    Collections.swap(cluster_labels, ii , ii + r.nextInt(cluster_labels.size() - ii));


            for(int j = 0; j < pivot_nodes.size(); j++) {
                int cur_pivot = pivot_nodes.get(j);

                // calculate probability of matching cur_pivot
                int num_matches = 0;
                for (int k = 0; k < MAX_R; k++) { // EXPERIMENT
                    if (cluster_labels.get(k).get(cur_node).equals(cluster_labels.get(k).get(cur_pivot))) {
                        num_matches += 1;
                    }
                    if (num_matches >= threshold) {
                        clusters.get(j).add(cur_node);
                        added = true;
                        break;
                    }
                }
                if (added)
                    break; 
            }
            if (!added) {                
                if (pivot_nodes.size() < MAX_PIVOTS) // limiting size of pivot set
                    pivot_nodes.add(cur_node);
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cur_node);
                clusters.add(new_cluster);
            }

        }
            
        return clusters;
    }





    // subroutine to shuffle pivots
    // how do I preserve original order?

    // https://stackoverflow.com/questions/4702036/take-n-random-elements-from-a-liste

    private static void pickNRandomElements(ArrayList<Integer> list, int n, Random r) {
        int length = list.size();

        if (length < n) return;

        //We don't need to shuffle the whole list
        for (int i = 0; i < n; i++)
        {
            int randIndex = r.nextInt(length - i) + i;
            Collections.swap(list, i , randIndex);            
        }
        
    }


    public static ArrayList<ArrayList<Integer>> permutationPKwik(ArrayList<ArrayList<Integer>> prob_matrix, int MAX_PIVOTS) {
        /*
            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node

            Shuffle pivot set at each step, and only query up to MAX_PIVOT nodes
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i); 

        Collections.shuffle(permutation);

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 
        ArrayList<Integer> pivot_nodes = new ArrayList<Integer>();
        pivot_nodes.add(permutation.get(0));

        clusters.add(new ArrayList<Integer>());
        clusters.get(0).add(permutation.get(0));

        Random myRand = new Random();
        for (int i = 1; i < num_nodes; i++) {   

            int cur_node = permutation.get(i);        
            boolean added = false;

            int LIMIT = Math.min(pivot_nodes.size(), MAX_PIVOTS);
            pickNRandomElements(pivot_nodes, LIMIT, myRand);
        
            for(int j = 0; j < LIMIT; j++) {
                int cur_pivot = pivot_nodes.get(j);
                if (prob_matrix.get(cur_node).contains(cur_pivot)) {
                    clusters.get(j).add(cur_node);
                    added = true;
                    break;
                }
            }
            if (!added) {
                pivot_nodes.add(cur_node);
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cur_node);
                clusters.add(new_cluster);
            }

        }
            
        return clusters;
    }


    // --- OLD CODE ---

    // read pKwikClusters from file

    public static ArrayList<ArrayList<Integer>> pKwikClusteringFromFile (String file_name, String path, int step) {

        try {
        BufferedReader myReader = new BufferedReader(new FileReader(file_name));
        String data = myReader.readLine();
        String[] split = data.split(" ");
        int num_nodes = Integer.parseInt(split[0]);
        myReader.close();

        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 
    
        // int loop1 = 0;
        // int loop2 = 0;
        // long start = 0;

        int cur_index = 0;

        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);
            // available_nodes.remove(cur_index);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);

            // get neighbors list of cur_node 
            int pivot_mod = cur_node % step;
            int pivot_div = cur_node / step;
            String file_path = path + "/" + pivot_div + ".txt";

            BufferedReader myReader2 = new BufferedReader(new FileReader(file_path));
            for (int k = 0; k < pivot_mod; k++)
                myReader2.readLine(); // is there a way to speed this up?
            String neighbors[] = myReader2.readLine().split(",");
            myReader2.close();

            // start = System.currentTimeMillis();
            int cur_neighbor;
            for (String neighbor : neighbors) {
                cur_neighbor = Integer.parseInt(neighbor);
                if (available_nodes.contains(cur_neighbor))
                    cur_cluster.add(cur_neighbor);
            }
            // loop1 += System.currentTimeMillis() - start;

            // start = System.currentTimeMillis();
            available_nodes.removeAll(cur_cluster);
            // for (int i = 1; i < cur_cluster.size(); i++) 
            //     available_nodes.remove(cur_cluster.get(i));
            // loop2 += System.currentTimeMillis() - start;
            clusters.add(cur_cluster);
        }
        // System.out.println("Loop 1 finished in: " + loop1 / 1000.0 + " s");
        // System.out.println("Loop 2 finished in: " + loop2 / 1000.0 + " s");
            
        return clusters;
        } catch (Exception e) {
            System.out.println("Oops! It broke. :(");
            return null;
        }

    }


    public static ArrayList<ArrayList<Integer>> permutationPKwikFromFile (String file_name, String path) {

        try {
        BufferedReader myReader = new BufferedReader(new FileReader(file_name));
        String data = myReader.readLine();
        String[] split = data.split(" ");
        int num_nodes = Integer.parseInt(split[0]);
        myReader.close();
        // System.out.println("okay");

        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i); 

        Collections.shuffle(permutation);

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 
        ArrayList<Integer> pivot_nodes = new ArrayList<Integer>();
        pivot_nodes.add(permutation.get(0));

        clusters.add(new ArrayList<Integer>());
        clusters.get(0).add(permutation.get(0));
        for (int i = 1; i < num_nodes; i++) {   

            int cur_node = permutation.get(i);        
            boolean added = false;

            // read neighbors of cur_node from file
        
            for(int j = 0; j < pivot_nodes.size(); j++) {
                int cur_pivot = pivot_nodes.get(j);
                int pivot_mod = cur_pivot % 10000;
                int pivot_div = cur_pivot / 10000;

                String file_path = path + "/" + pivot_div + ".txt";
                
                // System.out.println(cur_pivot);
                    BufferedReader myReader2 = new BufferedReader(new FileReader(file_path));
                    for (int k = 0; k < pivot_mod; k++)
                        myReader2.readLine(); // is there a way to speed this up?
                    String neighbors[] = myReader2.readLine().split(",");
                    myReader2.close();
                    for (String num : neighbors) {
                        if (Integer.parseInt(num) == cur_pivot) {
                            clusters.get(j).add(cur_node);
                            added = true;
                            break;
                        }
                    }
                

                /*
                if (prob_matrix.get(cur_node).contains(cur_pivot)) {
                    clusters.get(j).add(cur_node);
                    added = true;
                    break;
                }
                */
            }
            if (!added) {
                pivot_nodes.add(cur_node);
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cur_node);
                clusters.add(new_cluster);
            }

        }
            
        return clusters;
    } catch (Exception e) {
        System.out.println("FAIL! :(");
        return null;

    }

    }


    
}
