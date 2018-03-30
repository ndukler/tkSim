setGeneric("calculatePosteriors", function(object,...) standardGeneric("calculatePosteriors"))

#' Calculate Posterior Probabilities for Infered Parameters
#'
#' Uses numeric methods to estimate posteriors for the infered parameters \code{alpha} (synthesis rate) and \code{beta} (degredation rate).
#' Currently uses a flat prior for both alpha and beta as a default. This function is written such that non-flat priors may be used in the future.
#' @param object A \linkS4class{basicKineticModel} object
#' @param alphaRange Scale factors used to calculate the upper and lower bounds of the parameter range explored for \code{alpha}.  These scale factors will
#' be applied to the infered value of \code{alpha}.  Must be defined as \code{c(lower,upper)}.
#' @param betaRange Scale factors used to calculate the upper and lower bounds of the parameter range explored for \code{beta}.  These scale factors will
#' be applied to the infered value of \code{beta}.  Must be defined as \code{c(lower,upper)}.
#' @param paramSpaceSize The total size of parameter space to numerically integrate over. Half of the parameter space will be given to each parameter.
#' @param logProbAlpha  A function that returns the log probability for a given value of \code{alpha}
#' @param lobProbBeta A function that returns the log probability for a given value of \code{beta}
#'
#' @name calculatePosteriors
#' @include  class-basicKineticModel.R llFactory.R logSumExp.R
#' @examples
#' EXAMPLE HERE
#' @export
setMethod("calculatePosteriors",signature(object="basicKineticModel"), function(object,alphaRange=numeric(2),betaRange=numeric(2),paramSpaceSize=10^4,dispByGene=T,logProbAlpha=NULL,logProbBeta=NULL)
{
  if(alphaRange[1]==0)
  {
    warning("No range specified for alpha, using default range of .25x to 2x alpha.")
    alphaRange=c(.25,2)
  }else if(alphaRange[1]>alphaRange[2])
    stop("Alpha range specified incorrectly. "+alphaRange[1]+" > "+alphaRange[2]+". Must be specified as (lower, upper)")
  else if(length(alphaRange)>2)
    stop("Too many arguments supplied for alpha range.  Must be a vector of length 2 in the form (lower, upper).")

  if(betaRange[1]==0)
  {
    warning("No range specified for beta, using default range of .25x to 2x beta.")
    betaRange=c(.25,2)
  }else if(betaRange[1]>betaRange[2])
    stop("beta range specified incorrectly. "+betaRange[1]+" > "+betaRange[2]+". Must be specified as (lower, upper)")
  else if(length(betaRange)>2)
    stop("Too many arguments supplied for beta range.  Must be a vector of length 2 in the form (lower, upper).")

  #    | betaRange[1]==0)
  # {
  #   stop("Must enter ranges for both alpha and beta. They should be in the form c(lower,upper) where lower and upper are multipliers that are applied to the inferred parameter to calculate the bounds")
  # }

  if(is.null(logProbAlpha))
  {
    logProbAlpha = function(x){rep(log(1/sqrt(paramSpaceSize)),length(x))}
  }

  if(is.null(logProbBeta))
  {
    logProbBeta = function(x){rep(log(1/sqrt(paramSpaceSize)),length(x))}
  }

  #generate likelyhood esitmators for each gene
  logLH = lapply(X=1:nrow(object@inferedParams), FUN=llFactory, object=object,dispByGene=dispByGene)
  posteriors = lapply(X=1:nrow(object@inferedParams),object=object,logLH=logLH,FUN=function(x,object,logLH)
  {
    alpha = object@inferedParams[x,"alpha"]
    beta = object@inferedParams[x,"beta"]
    aMax = alphaRange[2]*alpha
    aMin = alphaRange[1]*alpha
    bMax = betaRange[2]*beta
    bMin = betaRange[1]*beta

    paramRange = expand.grid(seq(aMin,aMax,length.out = sqrt(paramSpaceSize)), seq(bMin,bMax,length.out = sqrt(paramSpaceSize)))
    numerator = apply(paramRange,1,function(y) logLH[[x]](y)) + logProbAlpha(paramRange[,1]) + logProbBeta(paramRange[,2])
    marginal = logSumExp(numerator)
    posterior = exp(numerator-marginal)
    res=cbind(paramRange,posterior=posterior)
    colnames(res) = c("alpha","beta","posterior")
    return(res)
  })
  object@posteriors = posteriors
  invisible(object)
})


