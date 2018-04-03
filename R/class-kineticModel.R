checkBase <- function(object) {
  errors <- character()
  if(nrow(object@data)>0){
    if(nrow(object@data)!=length(object@ids)){
      errors=c(errors,"There must be the same number of data entries as ids")
    }
    if(nrow(object@spikeIns)==0)
          errors = c(errors,"Must provide BOTH data AND spike ins. @spikeIns is empty.")
    if(ncol(object@data)!=nrow(object@expMetadata)){
      errors=c(errors,"There must be the same number of rows in the design table as there are columns in the data matrix")
    }
    if(length(object@normFactors)>0){
      if(length(object@normFactors)!= ncol(object@data)){
        errors=c(errors,"There must be the same number of normalization factors as there are columns in the data matrix.")
      }
      if(!is.numeric(object@normFactors)){
        errors=c(errors,"The normailzation factors must be numeric.")
      }
    }
  }
  if(!is.null(object@dispersionModel(1))){
    if(length(object@dispersionModel(1:10))!=10){
      errors=c(errors,"The dispersion model function must produce a vector of equal length to its input.")
    }
  }
  if (length(errors) == 0) TRUE else errors
}


methods::setClass(Class = "kineticModel",
                  representation = representation(ids="character",times="numeric",simData="matrix",equlibVals="numeric",
                                                  data="matrix",expMetadata="data.frame",normFactors="numeric",dispersionModel="function",spikeIns="matrix",
                                                  inferenceResults="list",inferedParams="matrix"),
                  prototype = methods::prototype(ids = NA_character_, times=NA_real_,equlibVals=NA_real_,dispersionModel=NULL),
                  validity = checkBase)
