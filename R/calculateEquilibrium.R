methods::setGeneric("calculateEquilibrium", function(object) {
  standardGeneric("calculateEquilibrium")
})

#' Calculate equilibrium transcript abundance
#'
#' Calculates equilibrium transcript abundance
#' @param object A basicKineticModel object
#' @name calculateEquilibrium
#' @include  class-basicKineticModel.R
#' @examples
#' ts=basicKineticModel(synthRate = 1:10,degRate = rep(0.3,10))
#' ts=calculateEquilibrium(ts)
#' @export
methods::setMethod("calculateEquilibrium", signature(object = "basicKineticModel"), function(object) {
  object@equlibVals <- object@synthRates/object@degRates
  return(object)
})
