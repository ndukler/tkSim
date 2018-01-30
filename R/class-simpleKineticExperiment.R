checkSimple <- function(object) {
  errors <- character()
  if(!is.na(object@synthesis.rates) && !is.na(object@degredation.rates)){
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
  } else if (nrow(object@data)==0){
    errors=c(errors,"Either synthesis AND degredation rates must be provided OR a dataset")
  }
  if (length(errors) == 0) TRUE else errors
}

#' Class simpleKineticExperiment
#'
#' Class \code{simpleKineticExperiment} defines a kinetic model dx/dt=a-b*x
#'
#' @name simpleKineticExperiment-class
#' @rdname simpleKineticExperiment-class
#' @include class-kineticExperiment.R
#' @exportClass simpleKineticExperiment
methods::setClass(Class = "simpleKineticExperiment",
                  representation = representation(initial.values = "numeric", synthesis.rates = "numeric",degredation.rates="numeric"),
                  prototype = methods::prototype(initial.values = NA_real_,synthesis.rates = NA_real_,degredation.rates=NA_real_),
                  validity = checkSimple,contains="kineticExperiment")
