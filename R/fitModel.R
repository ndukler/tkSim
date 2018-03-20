methods::setGeneric("fitModel", function(object,dispersionModel) {
  standardGeneric("fitModel")
})

#' Fits kinetic model
#'
#' Fits parameters for kinetic model given an dispersion model (NOT IMPLEMENTED)
#' @param object A basicKineticModel object
#' @param dispersionModel A function that takes predicted and observed values and computes the probability of the observation given the prediction
#' @include  class-basicKineticModel.R
#' @examples
#' ts=basicKineticModel(synthRate = 1:10,degRate = rep(0.3,10))
#' @export
methods::setMethod("fitModel", signature(object = "basicKineticModel"), function(object,dispersionModel) {
  ## Check if an dispersionModel is needed. Then, if an dispersion model is included, check for validity and update dispersionModel
  if(is.null(dispersionModel)){
    if(is.null(object@dispersionModel(1))){
      stop("There is no pre-specified kineticModel dispersion model so an dispersion model must be provided.")
    }
  } else if(is.function(dispersionModel)){
    if(length(dispersionModel(1:10))==10 && is.numeric(dispersionModel(1:10))){
      object@dispersionModel=dispersionModel #!!! Not permanent unless object is returned
    } else {
      stop("dispersionModel function must produce a *numeric* vector of the same length as the input.")
    }
  } else {
    stop("dispersionModel must be function")
  }
  predictAbundance(object,times=unique(object@expMetadata$time)) #!!! assumes first column will be named time

  })

geneNLL <- function(synthesis.rate,degredation.rate,initVals,times,data){   #!!! why is there a function defined here? Globally accessable?
  e.mu=exp(-degredation.rate*object@times)*(initVals-synthesis.rate/degredation.rate)+(synthesis.rate/degredation.rate)
}
