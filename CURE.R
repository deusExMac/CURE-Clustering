
#############################################################################################
#
#  Main script to execute the CURE algorithm. 
# 
#  The implementation of CURE follows the design described in section 7.4.1 (pg 275) 
#  of [2] and [3]. CURE allows for oddly shaped clusters.
#
#  Folder testdata/ contains csv files to test CURE clustering. Some of the 
#  files in testdata/ were downloaded from [4].
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
#    4) https://github.com/Kchu/CURE-cluster-python 
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
packages <- c("clevr", "clusterCrit", "ascii", "randomcoloR")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  cat('Installing following required packages: ', paste(packages[!installed_packages], collapse=' '), '\n', sep='')
  userResponse <- readline('Press enter to install. Press q to quit without installing >>>')
  if (userResponse=='q'){
      cat('Libraries not installed.')
      stop('Terminating due to user request')
  }
  
  install.packages(packages[!installed_packages])
}




# Package loading ---------------------------------------------------------


library(debugr) # debugging output



# printing formatted tables
library(ascii)
library(randomcoloR) # drawing representatives for each cluster

# Configuration related and some helper functions. 
# Contains also a global variable maintaining the configuration/settings 
# that will be used by CURE during execution.
source("configuration.R")

# Actual CURE implementation
source("CURE_Implementation.R")






# Start of execution ------------------------------------------------------


# Preprocessing -----------------------------------------------------------

#
#
# Initializing properly some important global variables.
# Review and modify these accordingly before any execution of this
# script.
#
#


# CSV file containing data to be clustered.
# This will also determine the settings in the yaml
# config file.
# Change this to cluster another file using CURE.
DATA_FILE_NAME <<- "testdata/jain.csv"

# This is important as it determines which section of the 
# configuration file will be active
Sys.setenv(R_CONFIG_ACTIVE = DATA_FILE_NAME)

# Load the default configuration file into the global variable defined in
# configuration.R. 
# Whenever a CURE setting is required, this variable is looked up. If no 
# configuration file exists, a default configuration is loaded.
CONFIGURATION <<- readConfigurationFile()

# Initialize the global variable specifying the columns of the data file 
# that will be used for clustering. An empty vector means all columns. 
# This variable is defined in configuration.R
TARGET_COLUMNS <<- getConfiguration(item="targetColumns", default=c())


# Actives dwatch() calls
if (getConfiguration(item="debugMessages", default=FALSE)){
    debugr_switchOn()
}else{
    debugr_switchOff()
}


cat("Executing with following configuration settings from [", DEFAULT_CONFIGURATION_FILE, ']\n', sep='')

# Configuration has already been loaded in utils.R and can be accessed using
# the global variable CONFIGURATION
printConfig(prefix='\t-')
cat('\n')
cat("Reading data from file: [", DATA_FILE_NAME, "]\n\n", sep='')




# CURE execution ----------------------------------------------------------


cat('[', getCurrectDateTime(), '] >>> Starting CURE on [', 
    DATA_FILE_NAME, ']\n', sep='')

# For the values returned by system.time and their interpretation, 
# see https://stackoverflow.com/questions/5688949/what-are-user-and-system-times-measuring-in-r-system-timeexp-output
cureBenchmark<-system.time(clusteredData <- CURE(DATA_FILE_NAME))

cat('[', getCurrectDateTime(), '] >>> CURE finished\n', sep='')


dataP <- clusteredData$clusteredDataPoints
representatives <- clusteredData$representivePoints

# Displaying clustering 
# plotting all points and cluster these were assigned to by CURE.
# ONLY for 2-D data

if (length(TARGET_COLUMNS) == 2){
    # Clustering will be displayed as 
    par(mfrow=c(1,2))
    plot(dataP[, TARGET_COLUMNS[1]], dataP[, TARGET_COLUMNS[2]], 
         main=paste("Clustering points using CURE", 
                    " (N=", as.character(nrow(dataP)), 
                    ")\n#init. clust=", getConfiguration(item='nClusters', default=3), 
                    " #Repr=", getConfiguration(item='nRepresentatives', default=5), 
                    " ALPHA=", getConfiguration(item='alpha', default=0.2), 
                    "\nMerge dist=", getConfiguration(item='clusterMergeDistance', default=0.88), sep=""),
        xlab=TARGET_COLUMNS[1], ylab=TARGET_COLUMNS[2], 
        pch=16, col=dataP$CURECluster)


    # Display also representative points
    for (i in unique(representatives$cluster)){
         clRep <- representatives[ representatives$cluster==i, ]
         points(clRep[, TARGET_COLUMNS[1]], clRep[, TARGET_COLUMNS[2]], col=randomColor(1, luminosity="bright"), pch=10, cex=2.0, lwd=2.8)
    }

} # if




# Clustering evaluation ---------------------------------------------------



cat('[', getCurrectDateTime(), '] >>> Executing kmeans n=', nrow(dataP[, TARGET_COLUMNS]), ' k=', getConfiguration(item='nClusters', default=3), '\n', sep='')

kmeansBenchmark<-system.time(kmcl <- kmeans(dataP[, TARGET_COLUMNS], centers=getConfiguration(item='nClusters', default=3), nstart=50, iter.max=530))
dataP$kmeanscluster <- kmcl$cluster

cat('[', getCurrectDateTime(), '] >>> kmeans done.\n', sep='')

if (length(TARGET_COLUMNS) == 2){
     plot(dataP[, TARGET_COLUMNS[1]], dataP[, TARGET_COLUMNS[2]], 
          main=paste("Clustering points using k-means", " (N=", as.character(nrow(dataP)), ")\n k=", getConfiguration(item='nClusters', default=3), sep=""),
          xlab=TARGET_COLUMNS[1], ylab=TARGET_COLUMNS[2], pch=16, col=dataP$kmeanscluster)
}

# Get the desired metrics/criteria from the configuration.
# Criteria are specified in the configuration file by name, as they 
# appear when calling getCriteriaNames(isInternal=F|T)
internalMetrics <- getConfiguration(item="internalClusterMetrics", default=NULL)
externalMetrics <- getConfiguration(item="externalClusterMetrics", default=NULL)



icmC<-internalClusteringMetrics('CURE', dataP, 
                                dataP$CURECluster, 
                                internalMetrics)
icmK<-internalClusteringMetrics('K-means', dataP, 
                                dataP$kmeanscluster, 
                                internalMetrics)


ecmC<-externalClusteringMetrics('CURE', dataP$cluster, 
                                dataP$CURECluster, 
                                externalMetrics)


ecmK<-externalClusteringMetrics('K-Means', dataP$cluster, 
                                dataP$kmeanscluster, 
                                externalMetrics)




# Prepare data frame to display metrics
CURERow <- append( c(cureBenchmark[3]), append(icmC$values, ecmC$values))
KMeansRow <- append( c(kmeansBenchmark[3]), append(icmK$values, ecmK$values) )

# Data frame from vector as row
nDF <- as.data.frame(t(CURERow))
nDF <- rbind(nDF, KMeansRow)

# Add header and row names
names(nDF) <- append( c("Time elapsed"), append(icmK$header, ecmK$header) )
row.names(nDF) <- c('CURE', 'k-means')


# Display evaluation data
print( ascii(nDF, digits=5), type="rest" )




