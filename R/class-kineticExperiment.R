checkBase <- function(object) {
  errors <- character()
   if(nrow(object@data)>0){
    if(nrow(object@data)!=length(object@ids)){
      errors=c(errors,"There must be the same number of data entries as ids")
    }
    if(ncol(object@data)!=nrow(object@col.info)){
      errors=c(errors,"There must be the same number of rows in the design table as there are columns in the data matrix")
    }
    if(length(object@design)==0){
      errors=c(errors,"There must be an expression to create the design matrix.")
    }
  }
  if (length(errors) == 0) TRUE else errors
}


methods::setClass(Class = "kineticExperiment",
                  representation = representation(ids="character",times="numeric",predicted.abundance="matrix",equilibrium.values="numeric",
                                                  data="matrix",col.info="data.frame",design="expression"),
                  prototype = methods::prototype(ids = NA_character_, times=NA_real_,equilibrium.values=NA_real_),
                  validity = checkBase)
