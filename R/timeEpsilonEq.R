methods::setGeneric("timeEpsilonEq", function(object,epsilon.percent) {
  standardGeneric("timeEpsilonEq")
})

#' Calculate time to near equilibrium
#'
#' Calculates the time at which each transcript will be within +/- epsilon% of its equilibrium value
#' @param object A simpleKineticSim object
#' @param epsilon.percent +/- Percent of equilibrium that must be reached
#' @name timeEpsilonEq
#' @include simpleKineticSimClass.R
#' @examples
#' ts=simpleKineticSim(syn.rate = 1:10,deg.rate = rep(0.3,10))
#' ts=timeEquilibrium(ts,0.01) # To 99% or 101% of equilibrium value
#' @export
methods::setMethod("timeEpsilonEq", signature(object = "simpleKineticSim"), function(object,epsilon.percent) {
  ev=calculateEquilibrium(object)@equilibrium.values
  ## Determine the stopping value depending on direction of approach
  fin=numeric(length(ev))
  fin[object@initial.values <= ev] = ev-(ev*epsilon.percent)
  fin[object@initial.values > ev] = ev+(ev*epsilon.percent)
  ## Return time at which value is reached
  return(log((fin-object@synthesis.rates/object@degredation.rates)/(object@initial.values-object@synthesis.rates/object@degredation.rates))/-object@degredation.rates)
})
