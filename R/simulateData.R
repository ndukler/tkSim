setGeneric("simulateData", function(object) standardGeneric("simulateData"))

#' Simulate Theoretical Data
#'
#' Simulates theoretical data based on given time points, synthesis rate, and degredation rate.
#' @param object A basicKineticModel object
#' @name simulateData
#' @include  class-basicKineticModel.R
#' @examples
#' bkm=basicKineticModel(synthRate = 1:10,degRate = rep(0.3,10))
#' bkm=simulateData(bkm)
#' @export
#'
setMethod('simulateData',signature(object='basicKineticModel'),function(object)
{
  #calculate fractional label abundance for each pair of syn and deg rates over all time points
  # e^(-bt)*(C - a/b) + a/b,   a=synrate b=degrate C=integration const
  if(length(times) > 1)
  {
    datout <- t( apply( cbind(object@synthRates, object@degRates, object@initVals), 1,
                  function(x){ exp(-x[2] * object@times) * (x[3] - x[1] / x[2]) + (x[1] / x[2]) }))

  } else if(length(times)==1)
  {
    datout <- matrix( apply( cbind(object@synthRates, object@degRates, object@initVals), 1,
                       function(x){ exp(-x[2] * object@times) * (x[3] - x[1] / x[2]) + (x[1] / x[2]) }),ncol = 1)
  }

  rownames(datout)=object@ids
  object@predictedAbundance <- dataout

  #calculate equlibirum state
  object@equlibVals <- object@synthRates/object@degRates
  return(object)
})
