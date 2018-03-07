checkBasic <- function(object) {
  errors <- character()
  if(!is.na(object@synthRates) && !is.na(object@degRates)){
    if(length(object@synthRates)!=length(object@degRates)){
      errors=c(errors,"There must be the same number of synthesis and degredation rates")
    }
    if(length(object@synthRates)!=length(object@ids)){
      errors=c(errors,"There must be the same number of synthesis rates as ids")
    }
    if(length(object@synthRates)!=length(object@initVals)){
      errors=c(errors,"There must be the same number of synthesis rates as initial values")
    }
    if(any(object@degRates>=1 | object@degRates <= 0)){
      errors=c(errors,"All degredation rates mut be between 0 and 1")
    }
    if(any(object@synthRates <= 0)){
      errors=c(errors,"All synthesis rates must be greater than 0")
    }
  } else if (nrow(object@data)==0){
    errors=c(errors,"Either synthesis AND degredation rates must be provided OR a dataset")
  }
  if (length(errors) == 0) TRUE else errors
}

#' Class basicKineticModel
#'
#' Class \code{basicKineticModel} defines a kinetic model dx/dt=a-b*x
#'
#' @name basicKineticModel-class
#' @rdname basicKineticModel-class
#' @include class-kineticModel.R
#' @exportClass basicKineticModel
methods::setClass(Class = "basicKineticModel",
                  representation = representation(initVals = "numeric", synthRates = "numeric",degRates="numeric"),
                  prototype = methods::prototype(initVals = NA_real_,synthRates = NA_real_,degRates=NA_real_),
                  validity = checkBasic,contains="kineticModel")
