import java.util.*;

public class RunFastCategoricalConsensusHybridRepresentatives {

    public static void main(String args[]){


        // ----- PARAMETERS -----
        String data_set = args[0];
        String delimiter = ",";
        int ROUNDS = 100;
        // ----------------------

        ArrayList<ArrayList<ArrayList<Integer>>> inputClusterings = Helper.cluster_categorical_data("Data/"+ data_set + "/graph.txt", delimiter);
        // ArrayList<ArrayList<Integer>> ground_truth1 = RunCorrelation.getGroundTruth("cor_stack_overflow", "\\s", 1);
        // ArrayList<ArrayList<ArrayList<Integer>>> gt_clustering = Helper.cluster_categorical_data("Data/" + data_set + "/ground_truth1.txt", "\\s");
        // ArrayList<ArrayList<Integer>> ground_truth1 = gt_clustering.get(0);


        // --- CONSENSUS CLUSTERING ---
        System.out.println("Creating Cluster Labels...");
        long consensusStart = System.currentTimeMillis();
        Consensus myConsensus = Consensus.consensusLabels(inputClusterings);
        long consensusTime = System.currentTimeMillis() - consensusStart;
        System.out.println("Cluster Labels finished in: " + consensusTime / 1000.0 + " s");

        ArrayList<HashMap<Integer, Integer>> cluster_labels = myConsensus.getClusterLabels();
        int num_nodes = myConsensus.getNumNodes();
        int num_attributes = cluster_labels.size();

        int x = num_attributes / 11;
        int[] attributes = {x, 2*x, 3*x, 4*x, 5*x, 6*x, 7*x, 8*x, 9*x, 10*x, num_attributes};
        boolean[] useRandom = {true}; // , false};

        System.out.println("Start: " + data_set + " with " + num_nodes + " entries and " + inputClusterings.size() + " attributes");
        System.out.println("Attribute sizes: ");
        for (int i = 0; i < inputClusterings.size(); i++)
            System.out.print(inputClusterings.get(i).size() + " ");
        System.out.println();
        System.out.println();

        for (int a = 0; a < attributes.length; a++) {

        

        // Collect numbers here
        double[] pivotTimes = new double[ROUNDS];
        long[] pivotScores = new long[ROUNDS];
        int[] pivotNumClusters = new int[ROUNDS];
        int[] pivotLargestCluster = new int[ROUNDS];
	long[] pivotDisagreement = new long[ROUNDS];

        double pKwikTimeTotal = 0;  
        int pKwikNumClusters = 0;
        int largestPKwikCluster = 0;
	long pKwikScores = 0;
	long pivot_disagreement1 = 0;

	double[] localTimes = new double[ROUNDS];
        long[] localScores = new long[ROUNDS];
        int[] localNumClusters = new int[ROUNDS];
        int[] localLargestCluster = new int[ROUNDS];

        double localTimeTotal = 0;
        int localNumClustersTotal = 0;
        int localLargestClusterTotal = 0;
        long localScoresTotal = 0;

        double[] voteTimes = new double[ROUNDS];
        long[] voteScores = new long[ROUNDS];
        int[] voteNumClusters = new int[ROUNDS];
        int[] voteLargestCluster = new int[ROUNDS];

        double voteTimeTotal = 0;
        int voteNumClustersTotal = 0;
        int voteLargestClusterTotal = 0;
        long voteScoresTotal = 0;

	long pCost = 0;
        long hCost = 0;

        int[] representatives = {num_nodes}; // , num_nodes}; // threshold for max cluster size
        int rep_rounds = representatives.length;
        double[] hybridTimeTotal = new double[rep_rounds];
        int[] hybridNumClusters = new int[rep_rounds];
        int[] largestHybridCluster = new int[rep_rounds];
	long[] hybridScores = new long[rep_rounds];
	long[] hybridDisagreement = new long[rep_rounds];

        double[][] hybridTimes = new double[rep_rounds][ROUNDS];
        long[][] hybridAllScores = new long[rep_rounds][ROUNDS];
        int[][] hybridAllNumClusters = new int[rep_rounds][ROUNDS];
        int[][] hybridLargestCluster = new int[rep_rounds][ROUNDS];
	long[][] hybridDisagreements = new long[rep_rounds][ROUNDS];

    ArrayList<Integer> all_nodes = new ArrayList<Integer>();
    for (int b = 0; b < num_nodes; b++)
        all_nodes.add(b);

        System.out.println("Start: " + attributes[a] + " attributes");
        for (int j = 0; j < ROUNDS; j++) {

        // RANDOM SAMPLE
        Collections.shuffle(cluster_labels);
        int cur_attributes = Math.min(cluster_labels.size(), attributes[a]);

        ArrayList<HashMap<Integer, Integer>> new_labels = new ArrayList<HashMap<Integer, Integer>>();
        for (int i = 0; i < cur_attributes; i++)
            new_labels.add(cluster_labels.get(i));

        System.out.print(j + ": ");

        long pKwikStart = System.currentTimeMillis();
        
	ArrayList<ArrayList<Integer>> result = PKwik.permutationPKwikConsensus(new_labels, num_nodes, num_nodes); 
        long pKwikTime = System.currentTimeMillis() - pKwikStart;
        pKwikTimeTotal += (pKwikTime / 1000.0);
        pivotTimes[j] = pKwikTime / 1000.0;
        // System.out.println("Pivot finished in: " + pKwikTime / 1000.0 + " s");
	
	pCost = 0; // Consensus.quickEditDist(result, ground_truth1);
        pivot_disagreement1 += pCost;
	pivotDisagreement[j] = pCost;

        pKwikNumClusters += result.size();
        pivotNumClusters[j] = result.size();

        int max_pkwik = 0;
        for (int i = 0; i < result.size(); i++) {
            if (result.get(i).size() > max_pkwik)
                max_pkwik = result.get(i).size();
        }
        largestPKwikCluster += max_pkwik;
        pivotLargestCluster[j] = max_pkwik;
        // System.out.println("Largest Pivot cluster: " + max_pkwik);

        // pKwikStart = System.currentTimeMillis();

        long p_score = Consensus.getEditDist(result, inputClusterings);
        // pKwikTime = System.currentTimeMillis() - pKwikStart;
        // System.out.println("pKwik Score finished in: " + pKwikTime / 1000.0 + " s");        
        // int score = Helper.large_graph_prob_edit_dist(result, prob_matrix);
        pKwikScores += p_score;
        pivotScores[j] = p_score;
	// System.out.println("Pivot Score: " + p_score);


        System.out.print("P ");


        // ------------- HYBRID ----------------

        for (int b = 0; b < rep_rounds; b++) {

            // MAKE A COPY OF PIVOT RESULT (not needed?)
	    ArrayList<ArrayList<Integer>> pivotCopy = new ArrayList<ArrayList<Integer>>();
	    for (int i = 0; i < result.size(); i++) {
                ArrayList<Integer> pivotClusterCopy = new ArrayList<Integer>();
		pivotClusterCopy.addAll(result.get(i));
		pivotCopy.add(pivotClusterCopy);
	    }

            long hybridStart = System.currentTimeMillis();
            ArrayList<ArrayList<Integer>> fix = Hybrid.large_graph_fix_clusters_consensus(pivotCopy, cluster_labels, cur_attributes, useRandom[b]);
            long hybridTime = System.currentTimeMillis() - hybridStart;
            hybridTimeTotal[b] += (hybridTime / 1000.0);
            hybridTimes[b][j] = hybridTime / 1000.0;
            // System.out.println("Hybrid finished in: " + hybridTime / 1000.0 + " s");
            // System.out.println();
	    
	    hCost = 0; // Consensus.quickEditDist(fix, ground_truth1);
            hybridDisagreement[b] += hCost;
	    hybridDisagreements[b][j] = hCost;

            // System.out.println("Num Pivot Clusters: " + result.size());
            hybridNumClusters[b] += fix.size();
            hybridAllNumClusters[b][j] = fix.size();
            // System.out.println("Num Hybrid Clusters: " + fix.size());
            // System.out.println();

        
            int max_hybrid = 0;
            for (int i = 0; i < fix.size(); i++) {
                if (fix.get(i).size() > max_hybrid)
                    max_hybrid = fix.get(i).size();
            }
            largestHybridCluster[b] += max_hybrid;
            hybridLargestCluster[b][j] = max_hybrid;
        
            // System.out.println("Largest Hybrid cluster: " + max_hybrid);
        
            // hybridStart = System.currentTimeMillis();
            long h_score = Consensus.getEditDist(fix, inputClusterings);
            // hybridTime = System.currentTimeMillis() - hybridStart;
            // System.out.println("Hybrid Score finished in: " + hybridTime / 1000.0 + " s");

            if (h_score > p_score) 
                h_score = p_score; // retain original clustering
    
            hybridScores[b] += h_score;
            hybridAllScores[b][j] = h_score;
	    // System.out.println("Hybrid Score: " + h_score);

            System.out.print("H" + b + " ");
        }


        //System.out.println(num_nodes);
        long localStart = System.currentTimeMillis();
       	ArrayList<ArrayList<Integer>> local_result = DNode.full_local_search_consensus(num_nodes, cluster_labels, result, cur_attributes); 
        long localTime = System.currentTimeMillis() - localStart;
        localTimeTotal += (localTime / 1000.0);
        localTimes[j] = localTime / 1000.0;

        localNumClustersTotal += local_result.size();
        localNumClusters[j] = local_result.size();

        int max_local = 0;
        for (int i = 0; i < local_result.size(); i++) {
            if (local_result.get(i).size() > max_local)
                max_local = local_result.get(i).size();
        }
        localLargestClusterTotal += max_local;
        localLargestCluster[j] = max_local;

        long l_score = Consensus.getEditDist(local_result, inputClusterings);
        if (l_score < p_score) {
            localScoresTotal += l_score;
            localScores[j] = l_score; 
        } else { // retain original clustering
            localScoresTotal += p_score;
            localScores[j] = p_score; 
        }

        System.out.print("L ");

        long voteStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> vote_result = DNode.adapted_node_consensus(cluster_labels, all_nodes, cur_attributes, num_nodes, true);
        long voteTime = System.currentTimeMillis() - voteStart;
        voteTimeTotal += (voteTime / 1000.0);
        voteTimes[j] = voteTime / 1000.0;

        voteNumClustersTotal += vote_result.size();
        voteNumClusters[j] = vote_result.size();

        int max_vote = 0;
        for (int i = 0; i < vote_result.size(); i++) {
            if (vote_result.get(i).size() > max_vote)
                max_vote = vote_result.get(i).size();
        }
        voteLargestClusterTotal += max_vote;
        voteLargestCluster[j] = max_vote;

        long v_score = Consensus.getEditDist(vote_result, inputClusterings);
        voteScoresTotal += v_score;
        voteScores[j] = v_score; 


        System.out.print("V");
        System.out.println();



        }

        
        System.out.println();
        System.out.println("Pivot times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Pivot scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Pivot num clusters: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotNumClusters[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Pivot max cluster size: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotLargestCluster[i] + " ");
        System.out.println();
                System.out.println();
        //System.out.println("Pivot disagreements: ");
        //for (int i = 0; i < ROUNDS; i++)
        //    System.out.print(pivotDisagreement[i] + " ");
        //System.out.println();

    
        //System.out.println("Finish");
        //System.out.println();
        //System.out.println("Num Attributes: " + new_labels.size());
        //System.out.println();
        System.out.println("Average Pivot time: " + pKwikTimeTotal / ((double) ROUNDS));
	System.out.println("Average Pivot score: " + pKwikScores / ((double) ROUNDS));
        System.out.println("Average Pivot num clusters: " + pKwikNumClusters / ((double) ROUNDS));
        System.out.println("Average Pivot max cluster size: " + largestPKwikCluster / ((double) ROUNDS));
        //System.out.println("Average Pivot disagreement: " + pivot_disagreement1 / ((double) ROUNDS));
	System.out.println();

        for (int b = 0; b < rep_rounds; b++) {
            System.out.println("ILS with up to " + representatives[b] + " representatives per cluster:");

        System.out.println();
        System.out.println("ILS"+b+" times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(hybridTimes[b][i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("ILS"+b+" scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(hybridAllScores[b][i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("ILS"+b+" num clusters: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(hybridAllNumClusters[b][i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("ILS"+b+" max cluster size: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(hybridLargestCluster[b][i] + " ");
        System.out.println();
            System.out.println();
        //System.out.println("Hybrid disagreements: ");
        //for (int i = 0; i < ROUNDS; i++)
        //    System.out.print(hybridDisagreements[b][i] + " ");
        //System.out.println();
        //System.out.println();


            System.out.println("Average ILS"+b+" time: " + hybridTimeTotal[b] / ((double) ROUNDS));
	    System.out.println("Average ILS"+b+" score: " + hybridScores[b] / ((double) ROUNDS));
            System.out.println("Percent Improvement: " + (pKwikScores - hybridScores[b]) / ((double) pKwikScores) * 100.0);
            System.out.println("Average ILS"+b+" num clusters: " + hybridNumClusters[b] / ((double) ROUNDS));
            System.out.println("Average ILS"+b+" max cluster size: " + largestHybridCluster[b] / ((double) ROUNDS));
            //System.out.println("Average hybrid disagreement: " + hybridDisagreement[b] / ((double) ROUNDS));
            //System.out.println("Percent Improvement: " + (pivot_disagreement1 - hybridDisagreement[b]) / ((double) pivot_disagreement1) * 100.0);
	    System.out.println();
        }

        System.out.println("Local times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(localTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Local scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(localScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Local num clusters: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(localNumClusters[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Local max cluster size: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(localLargestCluster[i] + " ");
        System.out.println();
        System.out.println();

        System.out.println("Average Local time: " + localTimeTotal / ((double) ROUNDS));
	System.out.println("Average Local score: " + localScoresTotal / ((double) ROUNDS));
        System.out.println("Percent Improvement: " + (pKwikScores - localScoresTotal) / ((double) pKwikScores) * 100.0);
        System.out.println("Average Local num clusters: " + localNumClustersTotal / ((double) ROUNDS));
        System.out.println("Average Local max cluster size: " + localLargestClusterTotal / ((double) ROUNDS));
	System.out.println();


    System.out.println("Vote times: ");
    for (int i = 0; i < ROUNDS; i++)
        System.out.print(voteTimes[i] + " ");
    System.out.println();
    System.out.println();
    System.out.println("Vote scores: ");
    for (int i = 0; i < ROUNDS; i++)
        System.out.print(voteScores[i] + " ");
    System.out.println();
    System.out.println();
    System.out.println("Vote num clusters: ");
    for (int i = 0; i < ROUNDS; i++)
        System.out.print(voteNumClusters[i] + " ");
    System.out.println();
    System.out.println();
    System.out.println("Vote max cluster size: ");
    for (int i = 0; i < ROUNDS; i++)
        System.out.print(voteLargestCluster[i] + " ");
    System.out.println();
    System.out.println();

    System.out.println("Average Vote time: " + voteTimeTotal / ((double) ROUNDS));
System.out.println("Average Vote score: " + voteScoresTotal / ((double) ROUNDS));
    System.out.println("Percent Improvement: " + (pKwikScores - voteScoresTotal) / ((double) pKwikScores) * 100.0);
    System.out.println("Average Vote num clusters: " + voteNumClustersTotal / ((double) ROUNDS));
    System.out.println("Average Vote max cluster size: " + voteLargestClusterTotal / ((double) ROUNDS));
System.out.println();


        }


    }
    
}
