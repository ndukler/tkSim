checkSimple <- function(object) {
  errors <- character()
  if(!is.na(object@synthesis.rates) || !is.na(object@degredation.rates)){
    if(length(object@synthesis.rates)!=length(object@degredation.rates)){
      errors=c(errors,"There must be the same number of synthesis and degredation rates")
    }
    if(length(object@synthesis.rates)!=length(object@ids)){
      errors=c(errors,"There must be the same number of synthesis rates as ids")
    }
    if(length(object@synthesis.rates)!=length(object@initial.values)){
      errors=c(errors,"There must be the same number of synthesis rates as initial values")
    }
    if(any(object@degredation.rates>=1 | object@degredation.rates <= 0)){
      errors=c(errors,"All degredation rates mut be between 0 and 1")
    }
    if(any(object@synthesis.rates <= 0)){
      errors=c(errors,"All synthesis rates must be greater than 0")
    }
  }
  if(nrow(object@data)>0){
    if(nrow(object@data)!=length(object@ids)){
      errors=c(errors,"There must be the same number of data entries as ids")
    }
    if(ncol(object@data)!=nrow(object@design)){
      errors=c(errors,"There must be the same number of rows in the design table as there are columns in the data matrix")
    }
    if(length(object@expression)==0){
      errors=c(errors,"There must be an expression to create the design matrix.")
    }
  }
  if (length(errors) == 0) TRUE else errors
}


methods::setClass(Class = "simpleKineticExperiment",
                  representation = representation(ids="character",initial.values = "numeric", synthesis.rates = "numeric",
                                      degredation.rates="numeric",times="numeric",abundances="matrix",equilibrium.values="numeric",
                                      data="matrix",design="data.frame",formula="expression"),
                  prototype = methods::prototype(ids = NA_character_, initial.values = NA_real_,synthesis.rates = NA_real_,
                                        degredation.rates=NA_real_,times=NA_real_,equilibrium.values=NA_real_),
                  validity = checkSimple)
