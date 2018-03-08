checkBase <- function(object) {
  errors <- character()
  if(nrow(object@data)>0){
    if(nrow(object@data)!=length(object@ids)){
      errors=c(errors,"There must be the same number of data entries as ids")
    }
    if(ncol(object@data)!=nrow(object@expMetadata)){
      errors=c(errors,"There must be the same number of rows in the design table as there are columns in the data matrix")
    }
    if(length(object@sizeFactors)>0){
      if(length(object@sizeFactors)!= ncol(object@data)){
        errors=c(errors,"There must be the same number of size factors as there are columns in the data matrix.")
      }
      if(!is.numeric(object@sizeFactors)){
        errors=c(errors,"The size factors must be numeric.")
      }
    }
  }
  if(!is.null(object@errorModel(1))){
    if(length(object@errorModel(1:10))!=10){
      errors=c(errors,"The error model function must produce a vector of equal length to its input.")
    }
  }
  if (length(errors) == 0) TRUE else errors
}


methods::setClass(Class = "kineticModel",
                  representation = representation(ids="character",times="numeric",simData="matrix",equlibVals="numeric",
                                                  data="matrix",expMetadata="data.frame",sizeFactors="numeric",errorModel="function"),
                  prototype = methods::prototype(ids = NA_character_, times=NA_real_,equlibVals=NA_real_,errorModel=NULL),
                  validity = checkBase)
