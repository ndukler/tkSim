checkBase <- function(object) {
  errors <- character()
  if(nrow(object@data)>0){
    if(nrow(object@data)!=length(object@ids)){
      errors=c(errors,"There must be the same number of data entries as ids")
    }
    if(ncol(object@data)!=nrow(object@col.info)){
      errors=c(errors,"There must be the same number of rows in the design table as there are columns in the data matrix")
    }
    if(length(object@size.factors)>0){
      if(length(object@size.factors)!= ncol(object@data)){
        errors=c(errors,"There must be the same number of size factors as there are columns in the data matrix.")
      }
      if(!is.numeric(object@size.factors)){
        errors=c(errors,"The size factors must be numeric.")
      }
    }
  }
  if(!is.null(object@error.model(1))){
    if(length(object@error.model(1:10))!=10){
      errors=c(errors,"The error.model function must produce a vector of equal length to its input.")
    }
  }
  if (length(errors) == 0) TRUE else errors
}


methods::setClass(Class = "kineticExperiment",
                  representation = representation(ids="character",times="numeric",predicted.abundance="matrix",equilibrium.values="numeric",
                                                  data="matrix",col.info="data.frame",size.factors="numeric",error.model="function"),
                  prototype = methods::prototype(ids = NA_character_, times=NA_real_,equilibrium.values=NA_real_),
                  validity = checkBase)
