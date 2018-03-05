methods::setGeneric("timeEpsilonEq", function(object,epsilon.percent) {
  standardGeneric("timeEpsilonEq")
})

#' Calculate time to near equilibrium
#'
#' Calculates the time at which each transcript will be within +/- epsilon% of its equilibrium value
#' @param object A basicKineticModel object
#' @param epsilon.percent +/- Percent of equilibrium that must be reached
#' @name timeEpsilonEq
#' @include class-basicKineticModel.R
#' @examples
#' ts=basicKineticModel(synthRate = 1:10,degRate = rep(0.3,10))
#' ts=timeEquilibrium(ts,0.01) # To 99% or 101% of equilibrium value
#' @export
methods::setMethod("timeEpsilonEq", signature(object = "basicKineticModel"), function(object,epsilon.percent) {
  ev=calculateEquilibrium(object)@equlibVals
  ## Determine the stopping value depending on direction of approach
  fin=numeric(length(ev))
  fin[object@initVals < ev] = ev[object@initVals < ev]-(ev[object@initVals < ev]*epsilon.percent)
  fin[object@initVals > ev] = ev[object@initVals > ev]+(ev[object@initVals > ev]*epsilon.percent)
  fin[object@initVals == ev] = 0
  ## Return time at which value is reached
  return(log((fin-object@synthRates/object@degRates)/(object@initVals-object@synthRates/object@degRates))/-object@degRates)
})
