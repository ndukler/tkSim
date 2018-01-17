methods::setGeneric("calculateEquilibrium", function(object) {
  standardGeneric("calculateEquilibrium")
})

#' Calculate equilibrium transcript abundance
#'
#' Calculates equilibrium transcript abundance
#' @param object A simpleKineticSim object
#' @name calculateEquilibrium
#' @include  class-simpleKineticExperiment.R
#' @examples
#' ts=simpleKineticExperiment(syn.rate = 1:10,deg.rate = rep(0.3,10))
#' ts=calculateEquilibrium(ts)
#' @export
methods::setMethod("calculateEquilibrium", signature(object = "simpleKineticExperiment"), function(object) {
  object@equilibrium.values <- object@synthesis.rates/object@degredation.rates
  return(object)
})
