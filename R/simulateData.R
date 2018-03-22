setGeneric("simulateData", function(object) standardGeneric("simulateData"))

#' Simulate Theoretical Data
#'
#' Simulates theoretical data based on given time points, synthesis rate, and degredation rate.
#' @param object A \linkS4class{basicKineticModel} object
#' @name simulateData
#' @include  class-basicKineticModel.R
#' @return  Returns a \linkS4class{basicKineticModel} object containing the simulated data with any parameters
#' overwritten by those supplied by the user.
#' @examples
#' bkm=basicKineticModel(synthRate = 1:10,degRate = rep(0.3,10), times=0:30)
#' bkm=simulateData(bkm)
#' @export
#'
setMethod('simulateData',signature(object='basicKineticModel'),function(object)
{
  #quick check for object validity
  validObject(object)
  #calculate fractional label abundance for each pair of syn and deg rates over all time points
  ##OLD: X(t) = e^(-bt)*(X(0) - a/b) + a/b,   a=synrate b=degrate X(0)=x at time 0
  ##If X(0)=0, X(t) = a(1-e^(-bt))/b
  if(is.na(object@times[1]))
  {
    stop("Error: Please set times in the basic kinetic model before running simulate data.")
  } else if(length(object@times) > 1)
  {
    dataout <- t( apply( cbind(object@synthRates, object@degRates, object@initVals), 1,
                  function(x){ exp(-x[2] * object@times) * (x[3] - x[1] / x[2]) + (x[1] / x[2]) }))

  } else if(length(times)==1)
  {
    dataout <- matrix( apply( cbind(object@synthRates, object@degRates, object@initVals), 1,
                       function(x){ exp(-x[2] * object@times) * (x[3] - x[1] / x[2]) + (x[1] / x[2]) }) ,ncol = 1)
  }

  rownames(dataout) = object@ids
  object@simData <- dataout

  #calculate equlibirum state
  object@equlibVals <- object@synthRates/object@degRates
  return(object)
})
