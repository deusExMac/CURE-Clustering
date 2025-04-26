#############################################################################################
#
#  Contains functions related to reading configuration files and defines some 
#  important global configuration variables including default settings. Moreover  
#  hosts some general purpose functions.
#
#  Must be included using source() by other files.
#
#  v0.7/mmt/April 2025
#
#############################################################################################

library(yaml)




# Global variables definition --------------------------------------------------------


# Actual configuration that will be used. Global variable.
CONFIGURATION <<- NULL

# Which columns to use when clustering
TARGET_COLUMNS <- c()


# Configuration should be a yaml formatted file.
# See https://www.cloudbees.com/blog/yaml-tutorial-everything-you-need-get-started
DEFAULT_CONFIGURATION_FILE = "CURE.config"



# Default configuration in case config cannot be read 
DEFAULT_CONFIGURATION = list('csvseparator'=',',
                             'targetcolumns'='',
                             'datasamplesize'=200,
                             'maxfilelines'=-1L,
                             'readchunksize'=130,
                             'nclusters'=3,
                             'nrepresentatives'=4,
                             'alpha'=0.72,
                             'mergedistance'=12.19)





# Function definitions ----------------------------------------------------



readConfigurationFile<-function(configFile=DEFAULT_CONFIGURATION_FILE, cfg="default"){
   
   CONFIGURATION <- c() # cleanup
   tryCatch({
              
              return(config::get(file=configFile,
                                 config=Sys.getenv("R_CONFIG_ACTIVE", cfg),
                                 use_parent = FALSE))
            },
            error=function(errMessage){
                           message('[WARNING] Error reading configuration file [', configFile, ']. Continuing with default settings')
                           message(errMessage)
                           readline('[Press enter to continue with default configuration or ctlr-C to quit] >')
                           return(DEFAULT_CONFIGURATION)
                  }
   )
}



getConfiguration = function (confList=CONFIGURATION, item, default=NULL) {
   value<-confList[[item]]
   if (is.null(value)) 
     return(default)
   else 
     return(value)
}



printConfig<-function(prefix='', cfg=CONFIGURATION){
     for (n in names(cfg))
          cat(prefix, n, ':', cfg[[n]], '\n')
}




# General purpose functions ---------------------------------------------------------




getCurrectDateTime<-function(fmt="%d/%m/%Y %H:%M:%OS"){
  return(format(Sys.time(), fmt))
}


vectorIsEmpty <- function(v) {
  if (length(v) == 0 | is.null(v)) {
      return(TRUE)
  } 
  
  return(FALSE)
}





