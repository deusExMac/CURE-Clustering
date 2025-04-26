
# 
# Install these packages first: testthat
#
# Run from console
#
# testthat::test_dir("tests")
#
# to execute tests



# Package installation ----------------------------------------------------



# Required R packages for this module
packages <- c("testthat", "devtools", "usethis", "here")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  cat('Installing following required packages: ', paste(packages[!installed_packages], collapse=', '), sep='')
  install.packages(packages[!installed_packages])
}



# Package loading ---------------------------------------------------------

library(testthat)

source("../CURE_Implementation.R", chdir = TRUE)






# Tests -------------------------------------------------------------------



test_that("Euclidean distance", {
  expect_equal(eucDist(c(1,2,3), c(4,5,6)), 5.196152, tolerance = 0.0000001)
  expect_equal(eucDist(c(1,2,3), c(1,2,3)), 0.0)
})



test_that("CURE point shrinking towards center", {
  expect_equal(movePointCloser(c(1,2,3), c(1,2,3)), c(1,2,3))
  expect_equal(movePointCloser(c(1,2,3), c(0.1,0.2,0.3), 0.3), c(0.73, 1.46, 2.19))
  expect_equal(movePointCloser(c(-1.67, 0.98, 3.01), c(0.73, 1.44, 2.37), 1.0), c(0.73, 1.44, 2.37))
  expect_equal(movePointCloser(c(-1.67, 0.98, 3.01), c(0.73, 1.44, 2.37), 1.3), c(1.450, 1.578, 2.178))
})




test_that("Testing configuration file reading and section (cfg) loading", {
  Sys.unsetenv("R_CONFIG_ACTIVE")
  
  expect_equal(readConfigurationFile('../CURE.config', cfg='testdata/syntheticDataN500x5.csv')$chunkSize, 150)
  expect_equal(readConfigurationFile('../CURE.config', cfg='testdata/syntheticDataN500x5.csv')$nClusters, 5)
  # Next should load default section of configuration file
  expect_equal(readConfigurationFile('../CURE.config', cfg='testdata/non_existent_file.csv')$outputFileName, "CURE_Clustered.csv")
  expect_equal(readConfigurationFile('../CURE.config', cfg='testdata/non_existent_file.csv')$nClusters, 7)
  expect_equal(readConfigurationFile('../CURE.config', cfg='testdata/non_existent_file.csv')$nRepresentatives, 16)
})


# Trivial testing of CURE clustering.

test_that("Testing CURE clustering of data with 2 variables. FILE testdata/jain.csv", {
  # Make sure...
  Sys.unsetenv("R_CONFIG_ACTIVE")
  
  CONFIGURATION <<-readConfigurationFile('../CURE.config', cfg='testdata/jain.csv')
  TARGET_COLUMNS <<- getConfiguration(item="targetColumns", default=c())
  expect_no_error(res <- CURE('../testdata/jain.csv')) 
})


test_that("Testing CURE clustering of data with 2 variables. FILE testdata/aggregation.csv", {
  # Make sure...
  Sys.unsetenv("R_CONFIG_ACTIVE")
  
  CONFIGURATION <<-readConfigurationFile('../CURE.config', cfg='testdata/aggregation.csv')
  TARGET_COLUMNS <<- getConfiguration(item="targetColumns", default=c())
  expect_no_error(res <- CURE('../testdata/aggregation.csv'))
})


test_that("Testing CURE clustering of data with 2 variables. FILE testdata/R15.csv", {
  # Make sure...
  Sys.unsetenv("R_CONFIG_ACTIVE")
  
  CONFIGURATION <<-readConfigurationFile('../CURE.config', cfg='testdata/R15.csv')
  TARGET_COLUMNS <<- getConfiguration(item="targetColumns", default=c())
  expect_no_error(res <- CURE('../testdata/R15.csv'))
})


test_that("Testing CURE clustering of data with more than 2 variables. FILE: testdata/syntheticDataN500x5.csv", {
  # Make sure...
  Sys.unsetenv("R_CONFIG_ACTIVE")
  
  CONFIGURATION <<- readConfigurationFile('../CURE.config', cfg='testdata/syntheticDataN500x5.csv')
  TARGET_COLUMNS <<- getConfiguration(item="targetColumns", default=c())
  expect_no_error(res <- CURE('../testdata/syntheticDataN500x5.csv'))
})


test_that("Testing CURE clustering of data with more than 2 variables. FILE: testdata/syntheticDataN500x5.csv", {
  # Make sure...
  Sys.unsetenv("R_CONFIG_ACTIVE")
  
  CONFIGURATION <<- readConfigurationFile('../CURE.config', cfg='testdata/syntheticData30000x6.csv')
  TARGET_COLUMNS <<- getConfiguration(item="targetColumns", default=c())
  expect_no_error(res <- CURE('../testdata/syntheticData30000x6.csv'))
})


