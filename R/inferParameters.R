setGeneric("inferParameters", function(object) standardGeneric("inferParameters"))

#' Infer Parameters from Read Data
#'
#' Infers \code{alpha} (synthesis rate) and \code{beta} (degredation rate) from sequencing read data (stored in the \code{data} slot) using the
#' basic kinetic model: \code{dx/dt = alpha - beta * x}.  If no read data exists in the provided \linkS4class{basicKineticModel} then simulated read data
#' will be generated using simulation data stored in \code{@@simData} before parameter inference. If no simulated data exists, it will be generated
#' before simulating read data and inferring parameters.
#' @param object A \linkS4class{basicKineticModel} object
#' @name inferParameters
#' @include  class-basicKineticModel.R getAbund.R
#' @examples
#' bkm=basicKineticModel(synthRate = 1:10,degRate = rep(0.3,10), times=0:30)
#' bkm=simulateData(bkm) #optional
#' bkm=simulateReads(bkm) #optional
#' bkm=inferParameters(bkm)
#' @export
#'

setMethod("inferParameters", signature(object="basicKineticModel"), function(object)
{
  ## Check if an dispersionModel is needed. Then, if an dispersion model is included, check for validity and update dispersionModel
  dispersionModel = object@dispersionModel
  if(is.null(dispersionModel)){
    if(is.null(dispersionModel(1))){
      stop("There is no pre-specified dispersion model for this object so an dispersion model must be provided.")
    }
  } else if(is.function(dispersionModel)){
    if(length(dispersionModel(1:10))==10 && is.numeric(dispersionModel(1:10))){
      object@dispersionModel=dispersionModel #!!! Not permanent unless object is returned
    } else {
      stop("dispersionModel function must produce a *numeric* vector of the same length as the input.")
    }
  } else {
    stop("dispersionModel must be function")
  }
  
  if(nrow(object@data)==0) #data not present
  {
    cat("\nNo experimental data detected, will now attempt to use simulated data.\n")
    if(nrow(object@simData)==0) #simulated data not present
    {
      cat("\nNo simulated data detected. Will now attempt to generate simulated data\n")
      object = simulateData(object)
      cat("\nData simulation successful. Now attempting to simulate reads from data.\n")
      object = simulateReads(object,expectedLibSize=10^6,replicates=3,spikeInSizes = 200,dispersionModel=function(x){rep(10^6,length(x))}) #cludge
      cat("\nRead simulation sucessful. Now inferring parameters from simulated data.\n")
    } else {

      cat("\nNo read data detected. Will now attempt to simulate read data.\n")
      object = simulateReads(object,expectedLibSize=10^6,replicates=3,spikeInSizes = 200,dispersionModel=function(x){rep(10^6,length(x))}) #cludge
      cat("\nRead simulation sucessful. Now inferring parameters from simulated data.\n")
    }
  }

  ##temp test for one gene
  nllFactory = function(geneIdx,object)
  {
    obs = object@data[geneIdx,]
    time=object@expMetadata$time
    initVal = object@initVals[geneIdx]
    normFactors = object@sizeFactors

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
      expMu = getAbund(params[1],params[2],time,initVal)*normFactors
      logProb = dnbinom(obs,mu=expMu,size=dispersion(expMu), log = T)
      return(-sum(logProb,na.rm=T)) #currently generating NaN for t=0
    })
  }

  #optimize to find Max Likelyhood of params
  nLL = lapply(X=1:nrow(object@data), FUN=nllFactory,object=object)
  paramRes = lapply(X=nLL,FUN=function(x){
              optim(par=c(1,0.2), fn=x, method="L-BFGS-B", lower=c(10^-5,10^-5), upper=c(Inf,1),hessian=T)
            })
  paramSummary = t(vapply(X=paramRes,FUN.VALUE=numeric(7),FUN=function(x){
                    #calculate 95% CI using hessian
                    fisher = solve(x$hessian) #returns inverse matrix  !!!FIX, should be -hessian
                    sigma = sqrt(diag(fisher))
                    upper = x$par+1.96*sigma
                    lower = x$par-1.96*sigma
                    c(upper[1],x$par[1],lower[1],upper[2],x$par[2],lower[2],x$convergence)
                  }))
  names(paramRes) = object@ids
  colnames(paramSummary) = c("a+CI","alpha","a-CI","b+CI","beta","b-CI","errorCode")
  rownames(paramSummary) = object@ids

  object@inferenceResults = paramRes
  object@inferedParams = paramSummary
  invisible(object)
})
