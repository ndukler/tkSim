checkBasic <- function(object) {
  errors <- character()
  if(!is.na(object@synthRates) && !is.na(object@degRates)){
    if(length(object@synthRates)!=length(object@degRates)){
      errors=c(errors,"There must be the same number of synthesis and degredation rates")
    }
    if(length(object@synthRates)!=length(object@ids)){
      errors=c(errors,"There must be the same number of synthesis rates as ids")
    }
    if(length(object@synthRates)!=length(object@initVals)){
      errors=c(errors,"There must be the same number of synthesis rates as initial values")
    }
    if(any(object@degRates>=1 | object@degRates <= 0)){
      errors=c(errors,"All degredation rates mut be between 0 and 1")
    }
    if(any(object@synthRates <= 0)){
      errors=c(errors,"All synthesis rates must be greater than 0")
    }
  } else if (nrow(object@data)==0){
    errors=c(errors,"Either synthesis AND degredation rates must be provided OR a dataset")
  }
  if (length(errors) == 0) TRUE else errors
}

#' Class basicKineticModel
#'
#' Class \code{basicKineticModel} defines a kinetic model dx/dt=a-b*x
#'
#' @name basicKineticModel-class
#' @rdname basicKineticModel-class
#' @include class-kineticModel.R
#' @exportClass basicKineticModel
methods::setClass(Class = "basicKineticModel",
                  representation = representation(initVals = "numeric", synthRates = "numeric",degRates="numeric",posteriors="list"),
                  prototype = methods::prototype(initVals = NA_real_,synthRates = NA_real_,degRates=NA_real_,posteriors=list()),
                  validity = checkBasic,contains="kineticModel")


#' Basic Kinetic Model Constructor
#'
#' Contructs a simple kinetic model object from synthesis and degredation rates
#' @param synthRate Synthesis rates. Must be greater than 0.
#' @param degRate Degredation rate. Must be between 0 and 1.
#' @param initAbund Initial abundance at time 0.
#' @param ids A character string id.
#' @param data Experimental or simulated matrix of count data
#' @param expMetadata A data.frame where each row corresponds to the column in data with the same index and the columns are different properties
#' @name basicKineticModel
#' @export
basicKineticModel <- function(times=NA_real_,synthRate=NA_real_,degRate=NA_real_,initAbund=NA_real_,ids=NA,data=NULL,expMetadata=data.frame(),dispersionModel=function(x){}){
  if(length(initAbund) == 1 && is.na(initAbund))
    initAbund=rep(0,length(synthRate))
  if(length(ids) == 1 && is.na(ids))
    ids=as.character(1:max(length(synthRate),nrow(data)))
  if(is.null(data))
    data=matrix(nrow=0,ncol=0)
  new("basicKineticModel",times=times,synthRates=synthRate,degRates=degRate,initVals=initAbund,ids=ids,expMetadata=expMetadata,
      data=data,dispersionModel=dispersionModel)
}
