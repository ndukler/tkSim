methods::setGeneric("calculateAbundance", function(object,times) {
    standardGeneric("calculateAbundance")
})

#' Calculate transcript abundance
#'
#' Calculates transcript abundance
#' @param object A simpleKineticSim object
#' @param times Times to calculate abundance values at.
#' @name calculateAbundance
#' @include simpleKineticSimClass.R
#' @examples
#' ts=simpleKineticSim(syn.rate = 1:10,deg.rate = rep(0.3,10))
#' ts=calculateAbundance(ts,0:30)
#' @export
methods::setMethod("calculateAbundance", signature(object = "simpleKineticSim",times="numeric"), function(object,times) {
  object@times <- times
  ab <- t(apply(cbind(object@synthesis.rates,object@degredation.rates,object@initial.values),1,
                             function(x) exp(-x[2]*object@times)*(x[3]-x[1]/x[2])+x[1]/x[2]))
  rownames(ab)=object@ids
  object@abundances <- ab
  return(object)
})
