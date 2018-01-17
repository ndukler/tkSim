methods::setGeneric("predictAbundance", function(object,times) {
    standardGeneric("predictAbundance")
})

#' Calculate transcript abundance
#'
#' Calculates transcript abundance
#' @param object A simpleKineticSim object
#' @param times Times to calculate abundance values at.
#' @name predictAbundance
#' @include  class-simpleKineticExperiment.R
#' @examples
#' ts=simpleKineticExperiment(syn.rate = 1:10,deg.rate = rep(0.3,10))
#' ts=predictAbundance(ts,0:30)
#' @export
methods::setMethod("predictAbundance", signature(object = "simpleKineticExperiment",times="numeric"), function(object,times) {
  object@times <- times
  ab <- t(apply(cbind(object@synthesis.rates,object@degredation.rates,object@initial.values),1,
                             function(x) exp(-x[2]*object@times)*(x[3]-x[1]/x[2])+x[1]/x[2]))
  rownames(ab)=object@ids
  object@predicted.abundance <- ab
  return(object)
})
