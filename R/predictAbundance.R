methods::setGeneric("predictAbundance", function(object,times) {
    standardGeneric("predictAbundance")
})

#' Calculate transcript abundance
#'
#' Calculates transcript abundance
#' @param object A basicKineticModel object
#' @param times Times to calculate abundance values at.
#' @name predictAbundance
#' @include  class-basicKineticModel.R
#' @examples
#' ts=basicKineticModel(synthRate = 1:10,degRate = rep(0.3,10))
#' ts=predictAbundance(ts,0:30)
#' @export
methods::setMethod("predictAbundance", signature(object = "basicKineticModel",times="numeric"), function(object,times) {
  object@times <- times
  if(length(times) > 1){
    ab <- t(apply(cbind(object@synthRates,object@degRates,object@initVals),1,
                             function(x){exp(-x[2]*object@times)*(x[3]-x[1]/x[2])+(x[1]/x[2])}))
  } else if(length(times)==1){
    ab <- matrix(apply(cbind(object@synthRates,object@degRates,object@initVals),1,
                  function(x){exp(-x[2]*object@times)*(x[3]-x[1]/x[2])+(x[1]/x[2])}),ncol = 1)
  }
  rownames(ab)=object@ids
  object@predictedAbundance <- ab
  return(object)
})
