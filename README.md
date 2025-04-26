# CURE-Clustering
 An implementation of the CURE (Clustering Using REpresentatives) clustering algorithm in R. 
 
 The implementation of CURE follows the design described in section 7.4.1 (pg 275) of [2] and [3]. 
 
 Folder testdata/ contains csv files to test CURE implementation. Some of the files in testdata/ were downloaded from [4].

# How to execute
Make sure to edit the configuration file CURE.config (written in the yaml serialization language) to set the appropriate settings (e.g. number of clusters, number of representatives etc). Each block can refer to specific file. To execute:

source CURE.R

By default, CURE.yaml file wil be loaded when execution starts.

# Configuration settings in yaml file
``targetColumns`` : LIST OF STRINGS. Which columns (names/header) from the csv file to use for clustering. An empty list means all columns

``randomDataSampleSize`` : INTEGER. Number of sampled lines to read from the csv file to start CURE. Lines are randomly selected.

``chunkSize`` : INTEGER. Number of lines to read each time from the csv file during the final step of point assignment to CURE clusters.

``nClusters`` : INTEGER. Number of initial clusters. Final number of clusters can be different due to merge operation. 

``nRepresentatives`` : INTEGER. Number of representative points to select from each cluster. 

``alpha`` : DOUBLE. Representative shrinking factor. Percetage of distance to shrink representative points towards cluster center. 

``clusterMergeDistance`` : DOUBLE. Representatives from different clusters with a distance smaller than clusterMergeDistance will have their clusters merged.



# References
1) Guha S, Rastogi R, Shim K. CURE: An efficient clustering algorithm for large databases[J]. ACM Sigmod record, 1998, 27(2): 73-84. doi: 10.1145/276305.276312

2) Leskovec, J., Rajaraman, A. and Ullman, J. D.: Mining of Massive Datasets. Cambridge University Press, 2014. Available from http://infolab.stanford.edu/~ullman/mmds/book0n.pdf

3) Lecture 62 — The CURE Algorithm (Advanced) | Stanford University, Available from: https://www.youtube.com/watch?v=JrOJspZ1CUw

4) https://github.com/Kchu/CURE-cluster-python
