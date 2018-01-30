methods::setGeneric("fitModel", function(object,error.model) {
  standardGeneric("fitModel")
})

#' Fits kinetic model
#'
#' Fits parameters for kinetic model given an error model (NOT IMPLEMENTED)
#' @param object A simpleKineticSim object
#' @param error.model A function that takes predicted and observed values and computes the probability of the observation given the prediction
#' @include  class-simpleKineticExperiment.R
#' @examples
#' ts=simpleKineticExperiment(syn.rate = 1:10,deg.rate = rep(0.3,10))
#' @export
methods::setMethod("fitModel", signature(object = "simpleKineticExperiment"), function(object,error.model) {
  ## Check if an error.model is needed. Then, if an error model is included, check for validity and update error.model
  if(is.null(error.model)){
    if(is.null(object@error.model(1))){
      stop("There is no pre-specified kineticExperiment error model so an error model must be provided.")
    }
  } else if(is.function(error.model)){
    if(length(error.model(1:10))==10 && is.numeric(error.model(1:10))){
      object@error.model=error.model
    } else {
      stop("error.model function must produce a *numeric* vector of the same length as the input.")
    }
  } else {
    stop("error.model must be function")
  }
  predictAbundance(object,times=unique(object@col.info$time))

  })

geneNLL <- function(synthesis.rate,degredation.rate,initial.values,times,data){
  e.mu=exp(-degredation.rate*object@times)*(initial.values-synthesis.rate/degredation.rate)+(synthesis.rate/degredation.rate)
}

