setGeneric("inferParameters", function(object) standardGeneric("inferParameters"))

#' Infer Parameters from Read Data
#'
#' Infers \code{alpha} (synthesis rate) and \code{beta} (degredation rate) from sequencing read data (stored in the \code{data} slot) using the
#' basic kinetic model: \code{dx/dt = alpha - beta * x}.  If no read data exists in the provided \linkS4class{basicKineticModel} then simulated read data
#' will be generated using simulation data stored in \code{@@simData} before parameter inference. If no simulated data exists, it will be generated
#' before simulating read data and inferring parameters.
#' @param object A \linkS4class{basicKineticModel} object
#' @name inferParameters
#' @include  class-basicKineticModel.R
#' @examples
#' bkm=basicKineticModel(synthRate = 1:10,degRate = rep(0.3,10), times=0:30)
#' bkm=simulateData(bkm) #optional
#' bkm=simulateReads(bkm) #optional
#' bkm=inferParameters(bkm)
#' @export
#'

setMethod("inferParameters", signature(object="basicKineticModel"), function(object)
{
  ## Check if an errorModel is needed. Then, if an error model is included, check for validity and update errorModel
  errorModel = object@errorModel
  if(is.null(errorModel)){
    if(is.null(errorModel(1))){
      stop("There is no pre-specified error model for this object so an error model must be provided.")
    }
  } else if(is.function(errorModel)){
    if(length(errorModel(1:10))==10 && is.numeric(errorModel(1:10))){
      object@errorModel=errorModel #!!! Not permanent unless object is returned
    } else {
      stop("errorModel function must produce a *numeric* vector of the same length as the input.")
    }
  } else {
    stop("errorModel must be function")
  }
  #predictAbundance(object,times=unique(object@expMetadata$time)) #!!! assumes first column will be named time
  if(nrow(object@data)==0) #data not present
  {
    cat("\nNo experimental data detected, will now attempt to use simulated data.\n")
    if(nrow(object@simData)==0) #simulated data not present
    {
      cat("\nNo simulated data detected. Will now attempt to generate simulated data\n")
      object = simulateData(object)
      cat("\nData simulation successful. Now attempting to simulate reads from data.\n")
      object = simulateReads(object,expectedLibSize = 3,replicates = 1,errorModel=function(x){rep(2,length(x))}) #cludge
      cat("\nRead simulation sucessful. Now inferring parameters from simulated data.\n")
    } else {

      cat("\nNo read data detected. Will now attempt to simulate read data.\n")
      object = simulateReads(object,expectedLibSize = 3,replicates = 1,errorModel=function(x){rep(2,length(x))}) #cludge
      cat("\nRead simulation sucessful. Now inferring parameters from simulated data.\n")
    }
  }
  if(object@times[1]==0)
    stop("Cannot have a time point of zero")
  # negbinomNLL = function(obs,times,initVals,FINISH)
  # {
  #
  # }

  ##temp test for one gene
  nllFactory = function(object,geneIdx)
  {
    obs = object@data[geneIdx,] #-1 to remove NaN at 0 from read simulation function
    time=object@times
    initVal = object@initVals[geneIdx]

    getAbund = function(alpha,beta,time,initVal)
    {
      return(exp(-time * beta) * (initVal - alpha / beta) + alpha / beta)
    }

    dispersion = function(mu)
    {
      return(rep(2,length(mu))) #cluge
    }

    return(function(params)
    {
      expMu = getAbund(params[1],params[2],time,initVal)
      # print(expMu)
      logProb = dnbinom(obs,mu=expMu,size=dispersion(expMu), log = T)
      # print(logProb)
      return(-sum(logProb,na.rm=T)) #currently generating NaN for t=0
    })
  }

  # geneNLL = function(params=c(alpha,beta),time,initVals,obs)
  # {
  #   #define model of theoretical data
  #   getAbund = function(alpha=params[1],beta=params[2],time,initVal)
  #   {
  #     return(exp(-time * beta) * (initVal - alpha / beta) + alpha / beta)
  #   }
  #
  #   #define size function (dispersion function)
  #   dispersion = function(mu)
  #   {
  #     ##TODO
  #   }
  #
  #   expMu = getAbund(alpha,beta,time,initVal)
  #   logProb = dnbinom(obs,mu=expMu,size=dispersion(expMu), log = T)
  #   return(-sum(logProb))
  # }

  # optimize(geneNLL)
  testNLL = nllFactory(object,10)
  optim(par=c(1,0.2),fn=testNLL,method="L-BFGS-B",lower = c(10^-5,10^-5), upper=c(Inf,1))
})

test=inferParameters(bkm)

# geneNLL <- function(synthesis.rate,degredation.rate,initVals,times,data){   #!!! why is there a function defined here? Globally accessable?
#   e.mu=exp(-degredation.rate*object@times)*(initVals-synthesis.rate/degredation.rate)+(synthesis.rate/degredation.rate)
# }
#

