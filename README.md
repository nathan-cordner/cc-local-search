Data Format:
-- Data files are stored in Data/[name of data set]/graph.txt
  -- All data sets are either publicly available, provided here, or the code used to generate them is provided
  -- We provide data set samples to show the format for the various experiments 
-- Delimiter is hard-coded in driver files
-- For 0/1 Correlation Clustering: 
  -- The first line of the file must contain the total number of nodes (anything that follows it on the line will be ignored)
  -- Rest of file lists positive edges as [node1] [node2]
-- For Consensus Clustering, no header line is needed
  -- File contains categorical data, which will be interpreted into input clusterings

To Compile: javac *.java
To Run: java [DriverName] [data set folder name]
  -- additional heap space may be needed for some experiments; increase the max heap size with the -Xmx flag

Drivers
-------

Hard coded parameters:
-- delimiter for data set
-- number of Pivot rounds, number of attributes used, etc. 

RunCorrelation.java
-- Input: positive edge list
-- Method: "neighborhood oracle" for Pivot algorithm; also runs LS and Vote
-- Sample data sets: cor_*

RunCorrelationTimed.java
-- Input: positive edge list
-- Method: "neighborhood oracle" for Pivot algorithm;
-- also runs ILS, "timed" LS (same amount of time as ILS), and full LS
-- Sample data sets: fb_*

RunCategoricalConsensus.java
-- Input: categorical data
-- Method: precomputed edges for Pivot algorithm
-- Sample data set: consensus_mushroom, consensus_fb_government

RunFastCategoricalConsensusHybridRepresentatives.java
-- Input: categorical data
-- Method: "Pivot on the fly" for consensus clustering, using attribute sampling
  -- Also runs multiple iterations of ILS, LS, and Vote per Pivot, wtih the same attributes that Pivot used 
-- Sample data sets: consensus_mushroom, consensus_fb_government

Code Files
----------

Consensus.java
-- Helper functions for Consensus clustering

DNode.java
-- Implementations of LocalSearch

Helper.java
-- Helper functions for reading data sets

Hybrid.java
-- ILS algorithm implementation (calls LocalSearch on input clusters)

PKwik.java
-- Pivot algorithm implementations

Data Generation
---------------

WritePivotClusters.java
-- runs the Pivot correlation algorithm a large number (~100) of times, writes resulting clusterings as categorical data
-- one attribute generated per clustering 
-- Sample data set: fb_government

Plots
-----

cc_ls_plots.ipynb
-- Jupyter Notebook for general CC data plots

local_search_consensus.ipynb, syn_cor_plots_new.ipynb
-- Jupyter Notebooks for consensus clustering plots
