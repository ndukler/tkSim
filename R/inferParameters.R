setGeneric("inferParameters", function(object,...) standardGeneric("inferParameters"))

#' Infer Parameters from Read Data
#' @description
#' Infers \code{alpha} (synthesis rate) and \code{beta} (degredation rate) from sequencing read data (stored in the \code{data} slot) using the
#' basic kinetic model: \code{dx/dt = alpha - beta * x}.
#' @description
#' If no read data exists in the provided \linkS4class{basicKineticModel} then simulated read data
#' will be generated using simulation data stored in \code{@@simData} before parameter inference. If no simulated data exists, it will be generated
#' before simulating read data and inferring parameters.
#'
#' @param object A \linkS4class{basicKineticModel} object
#' @param dispersionModel  A disperson model to use for inference. If not specified, the dispersion model stored in \code{object} will be used instead.
#' Must be specified as a function of the gene being analized if \code{dispByGene = TRUE} or as a function of the mean of the distribution if \code{dispByGene = FALSE}.
#' See the return values of \code{\link{estimateDispersions}} for examples of these two kinds of functions.
#' @param dispByGene Boolean controlling the expected nature of the \code{dispersionModel}. See \code{dispersionModel} description for more details.
#'
#' @name inferParameters
#' @include  class-basicKineticModel.R getAbund.R nllFactory.R
#' @examples
#' ##setup
#' bkm=basicKineticModel(times=0:30, synthRate = 1:10,degRate = rep(0.3,10))
#' bkm=simulateData(bkm) #optional
#' bkm=simulateReads(bkm,expectedLibSize=10^6,replicates=3,spikeInSizes = 200,dispersionModel=function(x){rep(10^3,length(x))}, dispByGene=F)
#'
#' ##infer params using same dispersion as simulated data
#' bkm=inferParameters(bkm,byGene=F)
#'
#' ##infer params using per-gene dispersion estimates from read data (dispersion estimates for each gene based on that gene's data alone)
#' bkm@dispersionModel = estimateDispersions(bkm,byGene=T)
#' bkm=inferParameters(bkm)
#'
#' ##infer params using mean-based dispersion estimates from read data (dispersion estimates based on entire data set)
#' bkm@dispersionModel = estimateDispersions(bkm,byGene=F)
#' bkm=inferParameters(bkm,byGene=F)
#'
#' @export
#'

setMethod("inferParameters", signature(object="basicKineticModel"), function(object,dispersionModel=NULL,dispByGene=T)
{
  validObject(object)

  ## Check if an dispersionModel is needed. Then, if an dispersion model is included, check for validity and update dispersionModel
  if(is.null(dispersionModel)){
    if(is.null(object@dispersionModel(1))){
      stop("There is no pre-specified dispersion model for this object so a dispersion model must be provided.")
    }
  } else if(is.function(dispersionModel)){
    if(length(dispersionModel(1:10))==10 && is.numeric(dispersionModel(1:10))){
      cat("\nUsing Provided Dispersion Model.  Replacing current dispersion model for returned object with user provided dispersion model.\n")
      object@dispersionModel=dispersionModel
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
      cat("Using:\n ExpectedLibSize:\t10^6\nReplicates:\t3\nSpikeInSizes:\t200\nDispersionModel:\tUser Provided Model (or if none then object's dispersionModel)")
      object = simulateReads(object,expectedLibSize=10^6,replicates=3,spikeInSizes = 200)
      cat("\nRead simulation sucessful. Now inferring parameters from simulated data.\n")
    } else {

      cat("\nNo read data detected. Will now attempt to simulate read data.\n")
      cat("Using:\n ExpectedLibSize:\t10^6\nReplicates:\t3\nSpikeInSizes:\t200\nDispersionModel:\tUser Provided Model (or if none then object's dispersionModel)")
      object = simulateReads(object,expectedLibSize=10^6,replicates=3,spikeInSizes = 200)
      cat("\nRead simulation sucessful. Now inferring parameters from simulated data.\n")
    }
  }


  #optimize to find Max Likelyhood of params
  nLL = lapply(X=1:nrow(object@data), FUN=nllFactory,object=object,dispByGene=dispByGene) #see nllFactory.R
  paramRes = lapply(X=nLL,FUN=function(x){
              optim(par=c(1,0.2), fn=x, method="L-BFGS-B", lower=c(10^-5,10^-5), upper=c(Inf,1),hessian=T)#,control=list(ndeps=c(10^-6,10^-6)))
            })
  # value = 1
  paramSummary = t(vapply(X=paramRes,FUN.VALUE=numeric(7),FUN=function(x){
                    #calculate 95% CI using hessian
                    fisher = solve(x$hessian) #returns inverse matrix
                    # tryCatch(sqrt(diag(fisher)),warning=function(w){cat("\nParams:\t",x$par,"\nFisher:\t",fisher,
                    #                                                      "\nGeneIdx:\t",value,"\nDat:\t",as.character(rate.comb[value,]),"\n")})
                    # value <<- value+1
                    sigma = sqrt(diag(fisher))
                    upper = x$par+1.96*sigma
                    lower = x$par-1.96*sigma
                    c(upper[1],x$par[1],lower[1],upper[2],x$par[2],lower[2],x$convergence)
                  }))
  #error messages for bad param inference
  badCI = which(rowSums(is.na(paramSummary))>0)
  if(length(badCI))
    cat("\nError: Bad parameter estimates for ",paste(object@ids[badCI],collapse=", "),".\nPlease check data quality for these samples.\n",sep="")

  names(paramRes) = object@ids
  colnames(paramSummary) = c("a+CI","alpha","a-CI","b+CI","beta","b-CI","errorCode")
  rownames(paramSummary) = object@ids

  object@inferenceResults = paramRes
  object@inferedParams = paramSummary
  invisible(object)
})
