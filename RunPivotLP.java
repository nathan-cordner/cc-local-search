// Copyright 2010-2022 Google LLC
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Minimal example to call the GLOP solver.
// [START program]
package com.google.ortools;
// [START import]
import com.google.ortools.Loader;
import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPObjective;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.linearsolver.MPVariable;

import java.util.*;

// [END import]

public final class RunPivotLP {
  public static void main(String[] args) {

    String data_set = args[0]; // "cor_landmarks";
    String delimiter = "\\s"; 
    ArrayList<ArrayList<Integer>> prob_matrix = Helper.read_large_network_relabel("Data/"+ data_set + "/graph.txt", delimiter);
    int num_nodes = prob_matrix.size();
    System.out.println("Data set with " + num_nodes + " nodes");

    Loader.loadNativeLibraries();
    // [START solver]
    // Create the linear solver with the GLOP backend.

    long lpStart = System.currentTimeMillis();

    MPSolver solver = MPSolver.createSolver("CLP");
    // [END solver]

    // [START variables]
    // Create the variables x and y.

    MPVariable[][] var_array = new MPVariable[num_nodes][num_nodes];

    for (int i = 0; i < num_nodes; i++) {
        for (int j = i+1; j < num_nodes; j++) {
            // String name = i + "_" + j;
            var_array[i][j] = solver.makeNumVar(0.0, 1.0, i + "_" + j); // + if i < j
            var_array[j][i] = solver.makeNumVar(0.0, 1.0, j + "_" + i); // - if i > j
        }
    }

    System.out.println("Number of variables = " + solver.numVariables());
    // [END variables]

    // [START constraints]

    // probability constraints
    for (int i = 0; i < num_nodes; i++) {
        for (int j = i + 1; j < num_nodes; j++) {
          MPConstraint ct = solver.makeConstraint(1.0, 1.0, i + "_" + j);
          ct.setCoefficient(var_array[i][j], 1); // assume all other variables are 0?
          ct.setCoefficient(var_array[j][i], 1);
        }
    }

    System.out.println("Constraints so far = " + solver.numConstraints());


    // triangle inequality constraints
    for (int i = 0; i < num_nodes; i++) {
        for (int j = i + 1; j < num_nodes; j++) {
            for (int k = j + 1; k < num_nodes; k++) {
                MPConstraint ct = solver.makeConstraint(0.0, 2.0, i + "_" + j + "_" + k);
                ct.setCoefficient(var_array[j][i], 1); // assume all other variables are 0?
                ct.setCoefficient(var_array[k][j], 1);     
                ct.setCoefficient(var_array[k][i], -1); 
            }

        }
    }


    System.out.println("Number of constraints = " + solver.numConstraints());
    // [END constraints]

    // [START objective]
    // Create the objective function, 3 * x + y.
    MPObjective objective = solver.objective();
    for (int i = 0; i < num_nodes; i++) {
        for (int j = i + 1; j < num_nodes; j++) {
            int prob = 0;
            if (prob_matrix.get(i).contains(j))
                prob = 1; // edge exists!
            objective.setCoefficient(var_array[i][j], 1 - prob);
            objective.setCoefficient(var_array[j][i], prob);
        }
    }

    objective.setMinimization();
    solver.setTimeLimit(600000); // 10 minute limit
    solver.solve();
    // [END solve]

    long lpTime = System.currentTimeMillis() - lpStart;

    System.out.println("LP Objective value = " + objective.value());
    System.out.println("Total time: " + (lpTime / 1000.0));
    System.out.println();

    int ROUNDS = 50;


        // Collect numbers here
        double[] pivotTimes = new double[ROUNDS];
        long[] pivotScores = new long[ROUNDS];

        double[] blendTimes = new double[ROUNDS];
        long[] blendScores = new long[ROUNDS];

        double pivotTimeTotal = 0;
        double blendTimeTotal = 0; // using "blend" for constrained RNode

        int pivotNumClusters = 0;
        int blendNumClusters = 0;

        int largestPivotCluster = 0;
        int largestBlendCluster = 0;

        long pivotScoreTotal = 0;
    	long blendScoreTotal = 0;

        for (int j = 0; j < ROUNDS; j++) {

            long pivotStart = System.currentTimeMillis();
            ArrayList<ArrayList<Integer>> pivot_result = PivotLP(var_array);
            long pivotTime = System.currentTimeMillis() - pivotStart;
            pivotTimeTotal += ((pivotTime + lpTime) / 1000.0);
            pivotTimes[j] = pivotTime + lpTime;
            // System.out.println("finish pivot " + j);
    
            long blendStart = System.currentTimeMillis();
            ArrayList<ArrayList<Integer>> blend_result = ChawlaLP(var_array, prob_matrix); // RoundingImproved(var_array, K);
            long blendTime = System.currentTimeMillis() - blendStart;
            blendTimeTotal += ((blendTime + lpTime) / 1000.0); 
            blendTimes[j] = blendTime + lpTime;
            // System.out.println("finish rnode " + j);
    
            pivotNumClusters += pivot_result.size();
            blendNumClusters += blend_result.size();
    
            int max_pivot = 0;
            for (int i = 0; i < pivot_result.size(); i++) {
                if (pivot_result.get(i).size() > max_pivot)
                    max_pivot = pivot_result.get(i).size();
            }
            largestPivotCluster += max_pivot;
    
            int max_blend = 0;
            for (int i = 0; i < blend_result.size(); i++) {
                if (blend_result.get(i).size() > max_blend)
                    max_blend = blend_result.get(i).size();
            }
            largestBlendCluster += max_blend;
    
            long pivot_score = Helper.quick_edit_dist(pivot_result, prob_matrix);
            pivotScoreTotal += pivot_score;
            pivotScores[j] = pivot_score;
    
            long blend_score = Helper.quick_edit_dist(blend_result, prob_matrix);
            blendScoreTotal += blend_score;
            blendScores[j] = blend_score;
    
            }
    
            System.out.println("Finish");
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
            
            System.out.println("Chawla times: ");
            for (int i = 0; i < ROUNDS; i++)
                System.out.print(blendTimes[i] + " ");
            System.out.println();
            System.out.println();
            System.out.println("Chawla scores: ");
            for (int i = 0; i < ROUNDS; i++)
                System.out.print(blendScores[i] + " ");
            System.out.println();
            System.out.println();
            
            System.out.println("Average Pivot time: " + pivotTimeTotal / ((double) ROUNDS));
            System.out.println("Average Chawla time: " + blendTimeTotal / ((double) ROUNDS));
            System.out.println();
            System.out.println("Average Pivot score: " + pivotScoreTotal / ((double) ROUNDS));
            System.out.println("Average Chawla score: " + blendScoreTotal / ((double) ROUNDS));
    
            System.out.println();
            System.out.println("Average Pivot num clusters: " + pivotNumClusters / ((double) ROUNDS));
            System.out.println("Average Chawla num clusters: " + blendNumClusters / ((double) ROUNDS));
    
            System.out.println();
            System.out.println("Average Pivot max cluster size: " + largestPivotCluster / ((double) ROUNDS));
            System.out.println("Average Chawla max cluster size: " + largestBlendCluster / ((double) ROUNDS));
    
            System.out.println();
    
    
    }
    


 
  


  public static ArrayList<ArrayList<Integer>> PivotLP(MPVariable[][] var_array) {
    /*
        prob_matrix is now a dictionary of lists 
        index of graph is a node, list contains all neighbors of node
    */
    
    int num_nodes = var_array.length;
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

        for (int i = 0; i < num_nodes; i++) {
            if (available_nodes.contains(i) && i != cur_node) {
                double prob = var_array[cur_node][i].solutionValue();
                if (i < cur_node)
                    prob = var_array[i][cur_node].solutionValue(); // get x+ value
                if (Math.random() < prob)
                    cur_cluster.add(i); 
            }
        }

        available_nodes.removeAll(cur_cluster);
        clusters.add(cur_cluster);            
    }
        
    return clusters;
}

private static double roundPositive (double x) {

    double a = 0.19;
    double b = 0.5095; // parameters set by paper 
    if (x < a)
        return 0.0;
    else if (x < b) {
        double result = (x - a) / (b - a);
        return Math.pow(result, 2);
    } else 
        return 1.0; 

}

public static ArrayList<ArrayList<Integer>> ChawlaLP(MPVariable[][] var_array, ArrayList<ArrayList<Integer>> prob_matrix) {
    /*
        prob_matrix is now a dictionary of lists 
        index of graph is a node, list contains all neighbors of node
    */


    
    int num_nodes = var_array.length;
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

        for (int i = 0; i < num_nodes; i++) {
            if (available_nodes.contains(i) && i != cur_node) {
                double prob = var_array[cur_node][i].solutionValue();
                if (i > cur_node)
                    prob = var_array[i][cur_node].solutionValue(); // get x+ value
                // check if edge between cur_node and i 
                if (prob_matrix.get(cur_node).contains(i))
                    prob = roundPositive(prob);
                // decide whether to add i to cur_cluster 
                if (Math.random() < (1 - prob))
                    cur_cluster.add(i); 
            }
        }

        available_nodes.removeAll(cur_cluster);
        clusters.add(cur_cluster);            
    }
        
    return clusters;
}





  // private BasicExample() {}
}
// [END program]
