import java.util.*;

public class RunCategoricalConsensus {

    public static void main(String args[]){

        String data_set = args[0];
        String delimiter = ","; 

        ArrayList<ArrayList<ArrayList<Integer>>> pKwikClusterings = Helper.cluster_categorical_data("Data/"+ data_set + "/graph.txt", delimiter);

        for (int i = 0; i < pKwikClusterings.size(); i++)
            System.out.print(pKwikClusterings.get(i).size() + " ");
        System.out.println();

        // --- CONSENSUS CLUSTERING ---
        System.out.println("Creating Consensus Matrix...");
        long consensusStart = System.currentTimeMillis();
        Helper myHelper = Consensus.consensusMatrix(pKwikClusterings);
        long consensusTime = System.currentTimeMillis() - consensusStart;
        System.out.println("Consensus Matrix finished in: " + consensusTime / 1000.0 + " s");

        ArrayList<ArrayList<Double[]>> match = myHelper.getMatches();
        ArrayList<ArrayList<Integer>> index_mapping = myHelper.getMapping();
        
        int ROUNDS = 10;
        double pKwikTimeTotal = 0;
        double hybridTimeTotal = 0;
  
        int pKwikNumClusters = 0;
        int hybridNumClusters = 0;

        int largestPKwikCluster = 0;
        int largestHybridCluster = 0;

	long pKwikScores = 0;
	long hybridScores = 0;

        System.out.println("Start");
        for (int j = 0; j < ROUNDS; j++) {

        long pKwikStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> result = PKwik.pKwikClusteringProb(match);
        long pKwikTime = System.currentTimeMillis() - pKwikStart;
        pKwikTimeTotal += (pKwikTime / 1000.0); 
        System.out.println("pKwik finished in: " + pKwikTime / 1000.0 + " s");

        long hybridStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> fix = result; // Hybrid.large_graph_fix_clusters_prob(result, match);
        long hybridTime = System.currentTimeMillis() - hybridStart;
        hybridTimeTotal += (hybridTime / 1000.0);
        System.out.println("Hybrid finished in: " + hybridTime / 1000.0 + " s");

        System.out.println();
        pKwikNumClusters += result.size();
        System.out.println("Num pKwik Clusters: " + result.size());
        hybridNumClusters += fix.size();
        System.out.println("Num Hybrid Clusters: " + fix.size());
        System.out.println();

        int max_pkwik = 0;
        for (int i = 0; i < result.size(); i++) {
            if (result.get(i).size() > max_pkwik)
                max_pkwik = result.get(i).size();
        }
        largestPKwikCluster += max_pkwik;
        System.out.println("Largest pKwik cluster: " + max_pkwik);
        
        int max_hybrid = 0;
        for (int i = 0; i < fix.size(); i++) {
            if (fix.get(i).size() > max_hybrid)
                max_hybrid = fix.get(i).size();
        }
        largestHybridCluster += max_hybrid;
        System.out.println("Largest Hybrid cluster: " + max_hybrid);
        

        
        pKwikStart = System.currentTimeMillis();

        long p_score = Consensus.getEditDist(Consensus.relabelClustering(result, index_mapping), pKwikClusterings);
        pKwikTime = System.currentTimeMillis() - pKwikStart;
        System.out.println("pKwik Score finished in: " + pKwikTime / 1000.0 + " s");        
        // int score = Helper.large_graph_prob_edit_dist(result, prob_matrix);
        pKwikScores += p_score;
	System.out.println("PKwik Score: " + p_score);
	
        
        hybridStart = System.currentTimeMillis();
        // int h_score = Helper.large_graph_prob_edit_dist_hash(fix, prob_matrix);

        long h_score = Consensus.getEditDist(Consensus.relabelClustering(fix, index_mapping), pKwikClusterings);
        hybridTime = System.currentTimeMillis() - hybridStart;
        System.out.println("Hybrid Score finished in: " + hybridTime / 1000.0 + " s");
        hybridScores += h_score;
	System.out.println("Hybrid Score: " + h_score);
        


        }
        System.out.println("Finish");
        System.out.println();
        System.out.println("Average pKwik time: " + pKwikTimeTotal / ((double) ROUNDS));
        System.out.println("Average hybrid time: " + hybridTimeTotal / ((double) ROUNDS));
        System.out.println();
        System.out.println("Average pKwik num clusters: " + pKwikNumClusters / ((double) ROUNDS));
        System.out.println("Average hybrid num clusters: " + hybridNumClusters / ((double) ROUNDS));
        System.out.println();
        System.out.println("Average pKwik max cluster size: " + largestPKwikCluster / ((double) ROUNDS));
        System.out.println("Average hybrid max cluster size: " + largestHybridCluster / ((double) ROUNDS));
     
        System.out.println();
	    System.out.println("Average pKwik score: " + pKwikScores / ((double) ROUNDS));
	    System.out.println("Average hybrid score: " + hybridScores / ((double) ROUNDS));
        System.out.println();




    }
    
}
