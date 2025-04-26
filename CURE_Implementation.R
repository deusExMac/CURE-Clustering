
#############################################################################################
#
#  Implements CURE algorithm [1] as described in section 7.4.1 (pg 275) of [2] and [3].
#  CURE allows for oddly shaped clusters.
#
#
#  References:
#
#    1) Guha S, Rastogi R, Shim K. CURE: An efficient clustering algorithm for large
#       databases[J]. ACM Sigmod record, 1998, 27(2): 73-84. doi: 10.1145/276305.276312

#    2) Leskovec, J., Rajaraman, A. and Ullman, J. D.: Mining of Massive Datasets. 
#       Cambridge University Press, 2014. Available from http://infolab.stanford.edu/~ullman/mmds/book0n.pdf
#
#    3) Lecture 62 — The CURE Algorithm (Advanced) | Stanford University, 
#       Available from: https://www.youtube.com/watch?v=JrOJspZ1CUw
#
#
#  v0.7/mmt/April 2025
#
#############################################################################################

# Close any graphic device left open
graphics.off()

# Cleanup the environment
# all=TRUE in ls() is optional
rm(list = ls())

# Clear console
cat("\014")





# Package installation ----------------------------------------------------


 
# Required R packages for this module
packages <- c("debugr", "config", "LaF", "FNN", "randomcoloR")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  cat('Installing following required packages: ', paste(packages[!installed_packages], collapse=' '), '\n', sep='')
  install.packages(packages[!installed_packages])
}





# Package loading ---------------------------------------------------------


library(debugr) # debugging output


library(LaF) # reading large files in chunks
library(FNN) # for knn

# Evaluation
library(clevr) # clustering evaluation
library(clusterCrit) # Calculating clustering metrics

library(randomcoloR) # for displaying clusters in different colors


# Configuration related. Configuration is loaded from 
# a default config file into a global variable defined here.
source("configuration.R")







# Function definitions ----------------------------------------------------




#' Executes hierarchical clustering on the data.
#'
#' @details  Calculates distance matrix. Requires numerical data.
#'           Does not normalize numerical values.
#'                            
#'
#' @param d DATA.FRAME with the data to be hierarchical clustered. 
#' @param nClusters INTEGER number of clusters to generate
#' @param hcM STRING method to calculate distance between clusters. Allowed
#'        values: same values as the method argument in hclust. Defaults to 
#'        complete
#' @param lbls  VECTOR of labels for identifying clusters. Defaults to 
#'        sequential numbers starting from 1         
#'  
hierarchicalClustering <- function(d, nClusters=0, hcM="complete", lbls=NULL){
  dataHClustering <- hclust(dist(d), method=hcM)
  if (!is.null(lbls)){
      dataHClustering$labels <- lbls
  }else{
      dataHClustering$labels <- 1:nrow(d)   
  }
  
  r <- cutree(dataHClustering, k=nClusters)
  return(r)
}






#' Method of finding the  points that are as far away from one another 
#' as possible within a cluster. Implements the variation of the 1st method  
#' described in section 7.3.2 (pg 267) of Mining Massive Datasets when
#' discussing K-Means.
#' 
#' @details  Calculates distance matrix. Requires numerical data.
#'           Does not normalize numerical values.
#'                            
#'
#' @param df DATA.FRAME data points belonging to the same cluster 
#' @param nPoints INTEGER number of representative points to select. Selects
#'        points whithin the cluster that are as far apart as possible. 
#' @param lbls VECTOR of labels for identifying clusters. Defaults to 
#'        sequential numbers starting from 1         
#' @return VECTOR with the rownames of the cluster points that have been selected
#'         and are far apart within the cluster 
#'         
farthestPoints<-function(df, nPoints=getConfiguration(item='nRepresentatives', default=5)){
  
  # TODO: Check this
  if (nrow(df) <= nPoints){
      return(rownames(df)) 
  }
  
  # Calculate distance matrix. 
  # Row and column labels/names are the 
  # positions of 
  dmat <- as.matrix(stats::dist(df))
  diag(dmat) <- NA
  
  
  # Get one random point
  k <- sample.int(nrow(dmat), 1)
  r <- rownames(dmat)[k]
  
  
  while (length(r) < nPoints ){
    dmatSub <- dmat[setdiff(colnames(dmat), r), r, drop = FALSE]
    
    if (length(r) > 1){
      maxMinDistance <- max(apply(dmatSub, 1, min))
    }else{
      maxMinDistance <- max(dmatSub) 
    }
    
    mPCoord <- which(dmatSub==maxMinDistance, arr.ind=T)
    r[length(r)+1] <- rownames(dmatSub)[mPCoord[1,'row']]
    
  }
  
  return(r)
}



#' Merges all clusters whose representatives are closer than a threshold.
#'
#' @details  Calculates distance matrix of representative points. Merges 
#'           clusters of the sampled data for which representatives are closer
#'           than a threshold. Representatives are estimated 
#'           again for the new, merged cluster. Applies shrinking to new 
#'           representatives. Merge operation finishes when no merge operation 
#'           happened.
#'           The new representatives are returned.
#'                            
#'
#' @param dataPoints DATA FRAME data frame with all sampled data. 
#' @param reprPoints DATA FRAME data frame with the initial representative 
#'                   points
#' @param mergeDist DOUBLE the representative distance below of which data belonging 
#'                  to the representative's cluster have to be merged to the same
#'                  cluster.
#' @param shrink BOOLEAN whether or not to shrink the new representatives of the 
#'               created clusters after a merge. Defaults to TRUE.
#'                  
#' @return data frame with the new representative points of the new clusters      
#'  
mergeClusters<-function(dataPoints, 
                        reprPoints, 
                        mergeDist=getConfiguration(item='clusterMergeDistance', default=0.88),
                        shrink=TRUE){
  
  
  if (mergeDist <= 0){
      dwatch(crit="1==1", objs=c("mergeDist"), msg="Not merging due to negative merge distance")
      return(reprPoints)
  }
  
  # Signals if a merge happen
  merged <- FALSE
  
  while (TRUE){
    
    # Generate distance matrix between representative points
    repDM <- as.matrix(dist(reprPoints[, TARGET_COLUMNS]))
    diag(repDM) <- NA
    
    # Get representative points closer together then maximum merge distance.
    # This will return a (row, col) that is relatively referencing the cells in
    # the distance matrix.
    closeRepPoints <- which(repDM < mergeDist, arr.ind=T)
    
    # Remove duplicates. If (i, j) is below threshold, so will (j, i)
    closeRepPoints <- closeRepPoints[!duplicated(t(apply(closeRepPoints, 1, sort.int, method='quick'))), , drop=FALSE]
    
    
    if (nrow(closeRepPoints) == 0){
      return(reprPoints)
    }
    
    # Iterate over all points below threshold
    for (i in 1:(nrow(closeRepPoints))) {
      
      rpt1Index <- rownames(repDM)[ closeRepPoints[i, "row"] ]
      rpt2Index <- colnames(repDM)[ closeRepPoints[i, "col"] ]
      
      
      if (reprPoints[rpt1Index, "cluster"] == 
          reprPoints[rpt2Index, "cluster"])
         next
      
      
      
      # merge clusters here, remove representative points and add new repr
      # and start merge operation again...  
      
      # Always assigning points with greatest cluster id to cluster with smallest 
      # cluster id. I.e. if clusters 2 and 6 should be merged, points in cluster 6
      # will be assigned to cluster 2.
      minC <- min(reprPoints[rpt1Index, "cluster"], reprPoints[rpt2Index, "cluster"])
      maxC <- max(reprPoints[rpt1Index, "cluster"], reprPoints[rpt2Index, "cluster"]) 
      
      dwatch(crit="1==1", expr=c("reprPoints[rpt1Index, 'cluster']", 
                                 "reprPoints[rpt2Index, 'cluster']",
                                 "repDM[rpt1Index, rpt2Index]"),  msg="Merging clusters and distance")
      
      
      
      
      # Actual merge operation of cluster points (not representatives)
      dataPoints[dataPoints[, 'cluster']==maxC, "cluster"] <- minC
      
      # Remove representative points
      idx<-which( (reprPoints[, 'cluster'] == minC) |  
                  (reprPoints[, 'cluster'] == maxC) )
      
      reprPoints <- reprPoints[-idx, ,drop=FALSE]
      
      
      # Get representatives from new cluster.
      # NOTE: nrpts are the indexes
      nrpts <- farthestPoints(dataPoints[dataPoints$cluster==minC, , drop=FALSE ], 
                              getConfiguration(item='nRepresentatives', default=5))
      
      
      if (!shrink){
        reprPoints <- rbind(reprPoints, dataPoints[nrpts, , drop=FALSE])
      } else {
              movedPoints <- as.data.frame(t(apply(dataPoints[nrpts, TARGET_COLUMNS],
                                                       1,
                                                       movePointCloser,
                                                       referencePoint=colMeans(dataPoints[dataPoints$cluster==minC, TARGET_COLUMNS]),
                                                       pctCloser=getConfiguration(item='alpha', default=0.2)
                                          )))
      
              movedPoints$cluster <- minC
              # Append new representative points
              reprPoints <- rbind(reprPoints, movedPoints)
      }
      
      merged <- TRUE
      
      # Merge happened; need to recalculate distance matrix
      break
    } # for  
    
    
    if (!merged){
      dwatch(crit="1==1", msg="Nothing merged. Quitting merge operation")
      break
    }
    
    # Start again
    merged <- FALSE
    
  } # while true
  
  dwatch(crit="1==1", expr=c("unique(reprPoints$cluster)"), msg="Clusters remaining after merge")
  
  
  if (length(TARGET_COLUMNS) == 2){
      plot(reprPoints[, TARGET_COLUMNS[1]], reprPoints[, TARGET_COLUMNS[2]], 
           main=paste("[mergeClusters] Representative points after merge"),
           xlab=TARGET_COLUMNS[1], ylab=TARGET_COLUMNS[2], pch=16,
           col=reprPoints[, "cluster"])
    
      
      for (c in unique(reprPoints$cluster)){
           rpt <- reprPoints[ reprPoints[ ,"cluster"]==c, ]
           
           text(colMeans(rpt[, TARGET_COLUMNS])[1], colMeans(rpt[, TARGET_COLUMNS])[2], 
                labels=c( paste0('cluster\n', c)), 
                col=randomColor(1, hue="blue", luminosity="dark"), cex=0.9)
      }
  }  
   
  return(reprPoints)
  
} #mergeClusters










#' Calculates the Euclidean distance between two specific points
#'
#' @details  Does not rely on existing R function. Points can have
#'           many dimensions.
#'           TODO: use existing R function
#'                            
#'
#' @param x1 VECTOR first data point.
#' @param x2 VECTOR second data point. 
#' @param lbls  VECTOR of labels for identifying clusters. Defaults to 
#'        sequential numbers starting from 1         
#'
#' @return DOUBLE the Euclidean distance between two points
#'    
eucDist <- function(x1, x2){
  return(sqrt(sum((x1 - x2)^2)))
}




#' Makes the size/length of a vector equal to one (1). 
#'
#' @details  Divides each vector element by the vector's norm. 
#'                            
#'
#' @param x VECTOR vector with numerical data only.
#'
#' @return DOUBLE a vector with length equal to one.
#'  
unitLengthNormalization <- function(x) { return( x / sqrt(sum(x^2)) )}






#' Calculates the new coordinates of points that are closer to a
#' reference point by a specified percentage in terms of distance (e.g. 
#' 0.2 means 20% closer to reference point). Used to move representatives
#' a specified percentage closer to the center of cluster they belong to.   
#'
#' @details  Calculates the Euclidean distance between the point and the  
#'           reference point and multiplies the specified percentage to get
#'           desired distance from reference. Normalizes
#'           the point difference to the unit vector and multiplies it with
#'           the desired distance. Subtracting this from the original point,
#'           will move it closer to the reference point. If original point and
#'           reference point coincidence, the original point is returned. 
#'           This method has taken from here: 
#'           https://math.stackexchange.com/q/175906                 
#'
#' @param p VECTOR point (representative) to move closer.
#' @param referencePoint VECTOR the point moving point p closer to (cluster center). 
#' @param pctCloser DOUBLE percentage of distance to move closer to reference point. 
#'        Can be any positive number. E.g. 0.35 would mean 35% closer to reference
#'        point (center)        
#'
#' @return VECTOR the point closer to the reference point (i.e. moved representative)
#'    
movePointCloser<-function(p, referencePoint, pctCloser){
  
  # If point falls exactly on the reference point
  # don't move closer since it will cause a divide by zero
  # when calculating the unit-length vector. In this
  # case, just return the point (or reference point)
  if (all(p==referencePoint)){
      return(p)
  }
  
  pointDist<-eucDist(p, referencePoint)*pctCloser
  pt <- p-referencePoint
  return( p-(pointDist*unitLengthNormalization(pt)) )
}






#' Get a random sample (lines) from a csv file by specifying a number of
#' of random lines to read from the file. 
#'
#' @details  File must be in csv format. Reads random lines from csv file; not consecutive
#'           lines. Line order is (also) random.  
#'
#' @param file STRING name of csv file to read random number of lines from.
#' @param n INTEGER number of lines to read. If number of lines to read from
#'          file exceeds number to total number of lines in file, script terminates.
#'
#' @return DATA FRAME a data frame with the  number of lines read.
#'
readRandomSampleFromFile <- function(file, n=130) {
  lf <- laf_open(detect_dm_csv(file, sep=getConfiguration(item='csvSeparator', default=","), header=TRUE, factor_fraction=-1))
  if (getConfiguration(item='maxFileLines', default=-1) < 0)
      maxLines <- nrow(lf)
  else {
      maxLines <- getConfiguration(item='maxFileLines', default=-1L)
  }
  
  if (maxLines < n){
      message( paste('Number of total lines in file ', getConfiguration(item='maxFileLines', default=-1)), 
               paste(' set fo smaller value than number of lines sampled ', n)
             ) 
      stop('Terminating.')
  }
  
  # These random lines will be read
  rLines <- sample(1:maxLines, n)
  dwatch(crit="length(rLines) >=0", objs=c("rLines", "file"), expr=c("length(rLines)"),  msg="Random sample lines read from file")
  return(LaF::read_lines(lf, rLines))
}




#' Reads the next consecutive chunk -consisting of a specified number of lines/rows- 
#' from a csv file.
#'
#' @details  File must be in csv format. Assumes file is already open. A "chunk" is 
#'           specified in terms of the number of lines/rows to read from the csv file.
#'           Each randomly selected line is returned only once. 
#'           Reading is done this way to support reading of huge files.  
#'
#' @param lafCon FILE CONNECTION a handle to the open file. Must have been opened with 
#'               LaF::laf_open() 
#' @param n INTEGER number of lines/rows to read from file.  
#'
#' @return DATA FRAME data frame with the number of lines read from file. If no more
#'         lines/rows can be read, NULL is returned.
#'
readChunk <- function(lafCon, n){
  if (is.null(lafCon)){
    message('No valid file connection. Cannot read file')
    stop('Terminating.')
  } 
  cdata <- next_block(lafCon, nrows=n)
  if (nrow(cdata)==0){
    return(NULL)
  }
  
  return(cdata)
}




#' Plots a specific cluster. Plot will include: the cluster data points, the 
#' initial representative points, the new representative points after 
#' shrinking and the cluster center.
#'
#' @details  Used to plot the clusters resulting from the initial hierarchical 
#'           clustering step. The plot will be generated only if the data has 
#'           2-dimensions. Cluster center will be calculated from points.       
#'
#' @param ptitle STRING a string that will be used as the main plot title 
#' @param cpoints DATA FRAME data frame containing the cluster data (points)  
#' @param rcPoints DATA FRAME data frame containing the initial cluster 
#'                 representative points
#' @param nrcPoints DATA FRAME data frame containing the new/final representative
#'                  points after shrinking. 
#'                                   
#' @return NOTHING Plot will be displayed in plot panel
#'
plotCluster<-function(ptitle, cpoints, rcPoints, nrcPoints){
  
  if (ncol(cpoints) != 2){
      message( paste('\t[WARNING] Cannot plot: data does not have 2 dimensions. (actual: ', ncol(cpoints), ' dimensions).', sep=''))
      return()
  }
  
  plot(cpoints[, TARGET_COLUMNS[1]], 
       cpoints[, TARGET_COLUMNS[2]], 
       main=paste(ptitle, " (n=", as.character(nrow(cpoints)), ")", sep=''),
       xlab="F1", ylab="F2", pch=19, cex=0.6)
  
  text(cpoints[, TARGET_COLUMNS[1]], cpoints[, TARGET_COLUMNS[2]], labels=rownames(cpoints), pos=1, cex=0.5)
  
  
  # Plot cluster center
  points(colMeans(cpoints)[TARGET_COLUMNS[1]], 
         colMeans(cpoints)[TARGET_COLUMNS[2]], 
         col = "red", bg="red",  pch=4, cex=1.5, lwd=2.1)
  
  
  # Highlight representative points for cluster
  points(rcPoints[, TARGET_COLUMNS[1]], 
         rcPoints[, TARGET_COLUMNS[2]], col = "red", bg="yellow",  pch=21, cex=1.1, lwd=1.8)
  
  # Plot new representative points that are closer to cluster center
  points(nrcPoints[, TARGET_COLUMNS[1]], nrcPoints[, TARGET_COLUMNS[2]], col = "green", bg="green",  pch=4, cex=1.5, lwd=1.8)
  text(nrcPoints[, TARGET_COLUMNS[1]], nrcPoints[, TARGET_COLUMNS[2]], labels=rownames(nrcPoints), col="violetred3", pos=1, cex=0.5)
  
  # Custom, transparent color for legend box 
  legend("bottomleft", 
         legend=c("Cluster data point", "Representative", "Moved representative", "Cluster center"),
         bg='bisque', pch = c(20, 1, 4, 4), col=c('black', 'red', 'green', 'red'),  cex=0.7,pt.cex=0.7
         )
}





#' Calculates and returns INTERNAL cluster metrics. 
#'
#' @details  Internal cluster metrics are specified in the arguments. Cluster 
#'           metrics must be specified by name, according to the naming scheme as 
#'           supported by the clusterCrit package and returned by 
#'           getCriteriaNames(isInternal=TRUE). Any of the names 
#'           -as they appear in the returned list- can be used. If no metrics
#'           are specified, no metrics will be calculated.        
#'
#' @param cLabel STRING a string to label the metrics. CURRENTLY NOT USED (24/04/2025) 
#' @param data DATA FRAME data frame containing the data (points) of a single
#'             cluster  
#' @param clustData VECTOR containing the cluster where each row in data has 
#'                  been assigned. 
#' @param iML VECTOR vector containing the names of internal metrics to calculate
#'                  based on data and clustData. 
#'                                   
#' @return LIST a list with two keys: header, a vector with the names of 
#'         metrics calculated and values, a vector with the calculated values 
#'         of the respective metrics in header.
#'
internalClusteringMetrics<-function(cLabel, data, clustData, iML){
  
  if (vectorIsEmpty(iML))
      return(list('header'=c(), 'values'=c()))
  
  
  icMetrics<-intCriteria(as.matrix(data), 
                         clustData,
                         iML)
  hVec <- c()
  mVec <- c()
  for (m in iML){
      mVec <- append(mVec, icMetrics[[tolower(m)]] )
      hVec <- append(hVec, m)
  }
  
  return(list('header'=hVec, 'values'=mVec))
}



#' Calculates and returns EXTERNAL cluster metrics. External here refers to 
#' a true/desired clustering of the same data.
#'
#' @details  External cluster metrics are specified in the arguments. Cluster 
#'           metrics must be specified by name, according to the naming scheme as 
#'           returned by getCriteriaNames(isInternal=FALSE). Any of the names 
#'           -as they appear in the returned list- can be used. The only metric
#'           that is calculated by this function but not supported by 
#'           clusterCrit is the v_measure (See https://towardsdatascience.com/v-measure-an-homogeneous-and-complete-clustering-ab5b1823d0ad/).         
#'
#' @param cLabel STRING a string to label the metrics. CURRENTLY NOT USED (24/04/2025) 
#' @param trueClust VECTOR vector containing the true clusters where each data/
#'                  point belongs to. 
#' @param calcClust VECTOR containing the cluster where each row/point has 
#'                  been assigned to by the clustering algorithm. 
#' @param eML VECTOR vector containing the names of external metrics to calculate
#'                  based on trueClust and calcClust. 
#'                                   
#' @return LIST a list with two keys: header, a vector with the names of 
#'         the metrics calculated and values, a vector with the calculated values 
#'         of the respective metrics in header.
#'
externalClusteringMetrics<-function(cLabel, trueClust, calcClust, eML){
  
  if (vectorIsEmpty(eML)){
    return(list('header'=c(), 'values'=c()))
  }
  
 
  # Filter out v_measure which is not supported
  # by clusterCrit package. v_measure will be 
  # calculated using the clevr package 
  filteredML <- eML[!tolower(eML) == "v_measure"]
  
  icMetrics <- NULL
  if (length(filteredML)>0) {
     if (!vectorIsEmpty(trueClust) & 
         !vectorIsEmpty(calcClust))
         icMetrics<-extCriteria(trueClust, 
                                calcClust,
                                filteredML)
  }
  
  
  
  hVec <- c()
  mVec <- c()
  for (m in tolower(eML)){
      if (m=='v_measure')
          mVec <- append(mVec, v_measure(trueClust, calcClust, beta = 1) )
      else
          mVec <- append(mVec, icMetrics[[m]] )
    
    hVec <- append(hVec, paste(m, '(ex.)', collapse='')) 
  }
  
  return(list('header'=hVec, 'values'=mVec))
}






#' Executes the CURE (Clustering Using REpresentatives) algorithm using the
#' data in a csv file.  
#'
#' @details  Settings for the execution of CURE (such as number of representatives)
#'           merge threshold, shrinking factor (alpha value) and others) are 
#'           specified in the global variable CONFIGURATION defined and loaded in
#'           utils.R.   
#'
#' @param fileName STRING the path to the csv file containing the data to be
#'                 clustered using CURE. 
#'                                   
#' @return LIST a list with three keys: clusteredOutputFileName, containing the 
#'         file name with the output of CURE clustering if the configuration 
#'         specifies saving the output to file (is empty if no such option is 
#'         set), clusteredDataPoints a data frame with all rows of the CURE 
#'         clustered data file, representivePoints a data frames with the 
#'         actual representative points used for clustering.   
#'
CURE<-function(fileName){
  
  ##########################################################
  #
  # Pass 1/2 of CURE (Initialization of CURE)
  #
  # The steps in section 7.4.1 in [2] will be carried out.
  #
  ##########################################################
  
  dwatch(crit="exists('fileName')", objs=c("fileName"), msg="File name to read data from:")
  dwatch(crit="1==1", 
         expr=c("getConfiguration(item='randomDataSampleSize', default=130)" ), 
         msg="Number of random samples from file")
  
  # Take a small sample of the data to cluster it in main memory 
  randomSample <- readRandomSampleFromFile(fileName, 
                                           getConfiguration(item='randomDataSampleSize', 
                                                            default=130))
  
  # An empty vector means all columns
  if (length(TARGET_COLUMNS)==0){
      # Setting the global variable since it is referenced throughout 
      # the file
      TARGET_COLUMNS <<- colnames(randomSample)
  }else{
         for (c in TARGET_COLUMNS)
              if ( !(c %in% colnames(randomSample)))
                  stop(paste0('Error: Invalid column ', c, '. Terminating') )
  }
  
  # Filter requested columns
  randomSample<-randomSample[, TARGET_COLUMNS]
  
  
  # Perform hierarchical clustering
  randomSample$cluster <- hierarchicalClustering(randomSample, 
                                                 getConfiguration(item='nClusters', default=3), 
                                                 getConfiguration(item='agglomerativeHierarchicalClusteringMethod', default="complete"), lbls=NULL)
  
  
  
  
  
  ##########################################################
  #
  # Pass 2/2 of CURE: Selecting representative points and 
  # shrinking
  #
  ##########################################################
  
  representativePoints <- data.frame()
  
  for (i in 1:getConfiguration(item='nClusters', default=3)){

    # Farthest points apart in cluster will be the representative points for
    # this cluster. 
    clusterRepresentativePoints <- farthestPoints(randomSample[randomSample$cluster==i, ,drop=FALSE], getConfiguration(item='nRepresentatives', default=5))
    
    
    dwatch(crit="any(randomSample[clusterRepresentativePoints, 'cluster', drop=FALSE] != i)", objs=c("clusterRepresentativePoints"), msg="Error! Clusters different than exprected.")
    
    
    
    # Move each of the representative points a fixed fraction of the distance
    # between its location and the centroid of its cluster
    shrunkRepresentativePoints <- as.data.frame(t(apply(randomSample[clusterRepresentativePoints, TARGET_COLUMNS, drop=FALSE],
                                                        1,
                                                        movePointCloser,
                                                        referencePoint=colMeans(randomSample[randomSample$cluster==i, TARGET_COLUMNS, drop=FALSE]),
                                                        pctCloser=getConfiguration(item='alpha', default=0.2)
    )))
    
    dwatch(crit="(eucDist(shrunkRepresentativePoints[1,TARGET_COLUMNS], colMeans(randomSample[randomSample$cluster==i, TARGET_COLUMNS, drop=FALSE])) / eucDist(shrunkRepresentativePoints[1, TARGET_COLUMNS], colMeans(randomSample[randomSample$cluster==i, TARGET_COLUMNS, drop=FALSE]))) == getConfiguration(item='alpha', default=0.2)",  msg="Error! Shriking not accurate.")
    
    # Aggregate shrunk representative points. Add also cluster each representative
    # belongs to
    shrunkRepresentativePoints$cluster <- i
    representativePoints<-rbind(representativePoints, shrunkRepresentativePoints) 
    
    # Plot cluster (debugging purposes only)
    plotCluster(paste('Cluster ', i, sep=''), 
                randomSample[randomSample$cluster==i, TARGET_COLUMNS], 
                randomSample[clusterRepresentativePoints, TARGET_COLUMNS ], 
                shrunkRepresentativePoints)
    
  } # for  
  
  
  
  ##########################################################
  #
  # Merging clusters that have representative points closer
  # together than a threshold. Representatives are recalculated
  # for merged clusters.
  #
  ##########################################################
  representativePoints <- mergeClusters(dataPoints=randomSample, 
                                        reprPoints=representativePoints,
                                        mergeDist=getConfiguration(item='clusterMergeDistance', default=0.88),
                                        shrink=TRUE)
  
  
  
  dwatch(crit="1==1", expr=c("getConfiguration(item='nClusters', default=3)",
                             "length(unique(representativePoints$cluster))"), msg="Initial and final (after merging) number of clusters")
  
  ##########################################################
  #
  # Final step: point assignment to clusters
  #
  ##########################################################
  
  LAFCON <- laf_open(detect_dm_csv(fileName, sep=getConfiguration(item='csvSeparator', 
                                                                  default=","), 
                                   header = TRUE, factor_fraction=-1))
  
  # Make sure we are at the start of file
  goto(LAFCON, 1)
  
  # Prepare the output file name if configured that way
  outputFilename <- ""
  if (getConfiguration(item='saveOutputToFile', default=FALSE)) {
      if (getConfiguration(item='appendDateTime', default=FALSE)) 
          outputFilename=paste0(getCurrectDateTime("%d-%m-%Y_%H%M%OS_"), getConfiguration(item='outputFileName', default="out.csv"))
      else {
          outputFilename=getConfiguration(item='outputFileName', default="out.csv")
      }
  }
  
  # All data points will be stored here. Makes the assumption that
  # all data points fit into memory.
  allDataPoints <-data.frame()
  
  while(TRUE){
    batch <- readChunk(LAFCON, getConfiguration(item='chunkSize', default=100))
    if (is.null(batch)){
      break
    }
    
    # Execute k-NN (k=1) on the loaded chunk from file and assign 
    # point to cluster of nearest representative
    c<-get.knnx(representativePoints[, TARGET_COLUMNS], 
                batch[, TARGET_COLUMNS], k=1)
    
    # Add cluster to each point of the current batch (where it belongs).
    batch$CURECluster <- representativePoints[c$nn.index, "cluster"]
    
    # Append batch to the data frame that will be returned. The data frame
    # will contain all rows of the source file.
    # NOTE: For testing purposes only, the assumption is made that the 
    # data file is small and fits into memory. For large files, batches
    # can be stored to file after being clustered.
    allDataPoints<-rbind(allDataPoints, batch)
    
    if (getConfiguration(item='saveOutputToFile', default=FALSE)){
      
        write.table(batch[, append(TARGET_COLUMNS, "CURECluster")], 
                    outputFilename, 
                    append=TRUE, sep=getConfiguration(item='csvSeparator', default=","),
                    col.names=!file.exists(outputFilename), 
                    row.names = FALSE)
    }
    
    
  } # while (TRUE)
  
  close(LAFCON)
  
  return(list('clusteredOutputFileName'=outputFilename, 
              'clusteredDataPoints'=allDataPoints, 
              'representivePoints'=representativePoints) )
  
} # End of CURE function

