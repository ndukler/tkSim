checkSimpleSim <- function(object) {
  errors <- character()
  if(is.na(object@synthesis.rates) || is.na(object@degredation.rates)){
    errors=c(errors,"Both synthesis and degredation rates must be specified")
  }
  else{
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
  if (length(errors) == 0) TRUE else errors
}


methods::setClass(Class = "simpleKineticSim",
                  representation = representation(ids="character",initial.values = "numeric", synthesis.rates = "numeric",
                                      degredation.rates="numeric",times="numeric",abundances="matrix",equilibrium.values="numeric"),
                  prototype = methods::prototype(ids = NA_character_, initial.values = NA_real_,synthesis.rates = NA_real_,
                                        degredation.rates=NA_real_,times=NA_real_,equilibrium.values=NA_real_),
                  validity = checkSimpleSim)
