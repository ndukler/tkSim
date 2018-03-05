methods::setGeneric("fitModel", function(object,errorModel) {
  standardGeneric("fitModel")
})

#' Fits kinetic model
#'
#' Fits parameters for kinetic model given an error model (NOT IMPLEMENTED)
#' @param object A basicKineticModel object
#' @param errorModel A function that takes predicted and observed values and computes the probability of the observation given the prediction
#' @include  class-basicKineticModel.R
#' @examples
#' ts=basicKineticModel(synthRate = 1:10,degRate = rep(0.3,10))
#' @export
methods::setMethod("fitModel", signature(object = "basicKineticModel"), function(object,errorModel) {
  ## Check if an errorModel is needed. Then, if an error model is included, check for validity and update errorModel
  if(is.null(errorModel)){
    if(is.null(object@errorModel(1))){
      stop("There is no pre-specified kineticModel error model so an error model must be provided.")
    }
  } else if(is.function(errorModel)){
    if(length(errorModel(1:10))==10 && is.numeric(errorModel(1:10))){
      object@errorModel=errorModel #!!! Not permanent unless object is returned
    } else {
      stop("errorModel function must produce a *numeric* vector of the same length as the input.")
    }
  } else {
    stop("errorModel must be function")
  }
  predictAbundance(object,times=unique(object@expMetadata$time)) #!!! assumes first column will be named time

  })

geneNLL <- function(synthesis.rate,degredation.rate,initVals,times,data){   #!!! why is there a function defined here? Globally accessable?
  e.mu=exp(-degredation.rate*object@times)*(initVals-synthesis.rate/degredation.rate)+(synthesis.rate/degredation.rate)
}
