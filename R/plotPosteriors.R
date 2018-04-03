#' Plot Posteriors for Infered Parameters
#' @description
#' Plots the results of the calculations from \code{\link{calculatePosteriors}} for a single gene. Allows for the selection of custom plotting bounds as
#' well as the recalculation of the posterior distribution within those bounds.
#' @description
#' Recalculation of the posterior is useful for obtaining better resolution for
#' plots that are significantly 'zoomed in' compared to the original posterior calculation. This also allows you to avoid having to call \code{\link{calculatePosteriors}}
#' if you are only interested in looking at one gene at a time.
#' @description
#' Note: \code{\link{calculatePosteriors}} must be run before calling this method unless \code{recalculate = T} in which case it \strong{does not} need to be run.
#' @name plotPosteriors
#' @usage plotPosteriors(object, geneIdx, alphaRange=NULL,
#'   betaRange=NULL, relative=T, recalculate=F, ...)
#'
#' @param object A \linkS4class{basicKineticModel} object
#' @param geneIdx The gene for which to plot the posterior distribution
#' @param alphaRange Determines the range of alpha values to plot over. Must be specified as \code{c(lower,upper)}. See \code{relative} for more information on how to specify the range. If not specified,
#' will default to the original range used when \code{\link{calculatePosteriors}} was called.
#' @param betaRange Determines the range of beta values to plot over. Must be specified as \code{c(lower,upper)}. See \code{relative} for more information on how to specify the range. If not specified,
#' will default to the original range used when \code{\link{calculatePosteriors}} was called.
#' @param relative Boolean used to determine whether parameter ranges are specified as relative or absolute. If \code{TRUE}, the ranges are treated as multiples of the
#' estimated parameter to plot over (eg. \code{alphaRange = c(.25,2)} will plot from .25*\code{alpha} to 2*\code{alpha}).  If \code{FALSE} then the ranges are treated
#' as absolute values of the specified parameter.
#' @param recalculate Boolean used to determine if the posterior distribution should be recalculated for the selected gene and parameter boundaries.
#' @param ... Additional arguments needed for recalculating the posterior distribution. Only used if \code{recalculate = TRUE}. See \code{\link{calculatePosteriors}}
#' for details on which parameters need to be supplied.
#'
#' @return A \code{ggplot2 plot} object containing the posterior plot
#'
#' @examples
#' ##setup
#' bkm=basicKineticModel(times=0:30,synthRate=1:10,degRate = rep(0.5,10))
#' bkm=simulateData(bkm)
#' bkm=simulateReads(bkm,expectedLibSize=10^6,replicates=3,spikeInSizes = 200,dispersionModel=function(x){rep(10^3,length(x))}) #CLUGE
#' bkm=inferParameters(bkm)
#' bkm=calculatePosteriors(bkm,alphaRange=c(.25,2))
#'
#' ##plot entire distribution
#' plotPosteriors(bkm,geneIdx=3)
#'
#' ##plot absolute subset of distribution
#' plotPosteriors(object=bkm,geneIdx=3,alphaRange=c(3000,4500),betaRange=c(.45,.55),relative=F)
#'
#' ##plot relative subset of distribution
#' plotPosteriors(object=bkm,geneIdx=3,alphaRange=c(.9,1.1))
#'
#' ##plot subset and recalculate posterior distribution (to obtain higher resolution)
#' #note: calculatePosteriors does NOT need to be run prior to plotting if recalculate=T
#' plotPosteriors(object=bkm,geneIdx=3,alphaRange=c(3000,4500),betaRange=c(.45,.55),relative=F,recalculate=T,paramSpaceSize=10^5)
#'
#'
#' @include llFactory.R logSumExp.R
#' @export

#error check param ranges
#make sure param ranges do not exceed posterior calcs
setGeneric("plotPosteriors", function(object,geneIdx,...) standardGeneric("plotPosteriors"))

setMethod("plotPosteriors", signature(object="basicKineticModel",geneIdx="numeric"), function(object,geneIdx,alphaRange=NULL,betaRange=NULL,relative=T,recalculate=F,...)
{
  if(is.null(object@posteriors[1])&!recalculate)
      stop("No posteriors detected.  Please run calculatePosteriors first.")
  # if(is.null(geneIdx))
  #   stop("Error: No gene selected.  Please specify the gene you want to plot using geneIdx.")
  if(length(geneIdx)>1)# else if(length(geneIdx)>1)
  {
    warning("Multiple genes specified, only plotting the first gene.")
    geneIdx=geneIdx[1]
  }

  if(!is.null(alphaRange) && (alphaRange[1]>alphaRange[2]))
    stop(paste0("Alpha range specified incorrectly. Lower: ",alphaRange[1]," > Upper: ",alphaRange[2],". Upper must be greater than Lower"))
  else if(length(alphaRange)>2)
    stop("Too many arguments supplied for alpha range.  Must be a vector of length 2 in the form (lower, upper).")

  if(!is.null(betaRange) && betaRange[1]>betaRange[2])
    stop(paste0("Beta range specified incorrectly. Lower: ",betaRange[1]," > Upper: ",betaRange[2],". Upper must be greater than Lower"))
  else if(length(betaRange)>2)
    stop("Too many arguments supplied for beta range.  Must be a vector of length 2 in the form (lower, upper).")

  posteriorData = object@posteriors[[geneIdx]]

  if(recalculate)
  {
    params = list(...)

    optionalParamNames = c("logProbAlpha","dispByGene","logProbBeta","paramSpaceSize")
    unusedParams = setdiff(names(params),optionalParamNames)
    if(length(unusedParams))
      warning('unused arguments: ', paste(lapply(unusedParams,function(x) paste0("(",x," = ",params[x],")")),collapse = ', '))

    if(!("dispByGene"%in%names(params)))
      error("If recalculating posterior, dispByGene must be provided.")

    #check for supplied params, provide defaults if necessary
    if(!("logProbAlpha"%in%names(params)))
    {
      #global basis
      aMax = alphaRange[2]*max(object@inferedParams[,'alpha'])
      aMin = alphaRange[1]*min(object@inferedParams[,'alpha'])
      logProbAlpha = function(x){rep(log(1/sqrt(paramSpaceSize)),length(x))} #log(1/(aMax-aMin)) == -log(aMax-aMin)
    }else
    {
      logProbAlpha = params$logProbAlpha
    }

    if(!('logProbBeta'%in%names(params)))
        logProbBeta = function(x){rep(log(1),length(x))}
    else
        logProbBeta = params$logProbBeta
    if(!("paramSpaceSize"%in%names(params)))
    {
      warning("No parameter space size specified, using the default of 10^4.")
      paramSpaceSize = 10^4
    }else
    {
      paramSpaceSize = params$paramSpaceSize
    }

    #define ranges for alpha and beta
    if(!is.null(alphaRange))
    {
      if(relative)
      {
        aMax = alphaRange[2]*object@inferedParams[geneIdx,"alpha"]
        aMin = alphaRange[1]*object@inferedParams[geneIdx,"alpha"]
      }else
      {
        aMax = alphaRange[2]
        aMin = alphaRange[1]
      }

    }else
    {
      aMax = max(posteriorData$alpha)
      aMin = min(posteriorData$alpha)
    }

    if(!is.null(betaRange))
    {
      if(relative)
      {
        bMax = betaRange[2]*object@inferedParams[geneIdx,"beta"]
        bMin = betaRange[1]*object@inferedParams[geneIdx,"beta"]
      }else
      {
        bMax = betaRange[2]
        bMin = betaRange[1]
      }
    } else
    {
      bMax = max(posteriorData$beta)
      bMin = min(posteriorData$beta)
    }

    #recalculate posteriors for this gene using newly defined bounds
    logLH = llFactory(geneIdx,object,params$dispByGene)
    paramRange = expand.grid(seq(aMin,aMax,length.out = sqrt(paramSpaceSize)), seq(bMin,bMax,length.out = sqrt(paramSpaceSize)))
    numerator = apply(paramRange,1,function(y) logLH(y)) + logProbAlpha(paramRange[,1]) + logProbBeta(paramRange[,2])
    marginal = logSumExp(numerator)
    posterior = exp(numerator-marginal)
    res=cbind(paramRange,posterior=posterior)
    colnames(res) = c("alpha","beta","posterior")
    posteriorData = res

  }



  #determine plot range
  if(relative)
  {
    if(!is.null(alphaRange))
    {
      xmin = object@inferedParams[geneIdx,"alpha"]*alphaRange[1]
      xmax = object@inferedParams[geneIdx,"alpha"]*alphaRange[2]
    }
    if(!is.null(betaRange))
    {
      ymin = object@inferedParams[geneIdx,"beta"]*betaRange[1]
      ymax = object@inferedParams[geneIdx,"beta"]*betaRange[2]
    }

    if(is.null(alphaRange) & is.null(betaRange))
    {
      idx = 1:nrow(posteriorData)
    }
    else if(is.null(betaRange))
    {
      idx = which(posteriorData$alpha<xmax & posteriorData$alpha>xmin)
    }else if(is.null(alphaRange))
    {
      idx = which(posteriorData$beta<ymax & posteriorData$beta>ymin)
    }else
    {
      idx = which((posteriorData$alpha<xmax & posteriorData$alpha>xmin)&(posteriorData$beta<ymax & posteriorData$beta>ymin))
    }
  } else
  {
    if(is.null(alphaRange) & is.null(betaRange))
    {
      idx = 1:nrow(posteriorData)
    }
    else if(is.null(betaRange))
    {
      idx = which(posteriorData$alpha<alphaRange[2] & posteriorData$alpha>alphaRange[1])
    }else if(is.null(alphaRange))
    {
      idx = which(posteriorData$beta<betaRange[2] & posteriorData$beta>betaRange[1])
    }else
    {
      idx = which((posteriorData$alpha<alphaRange[2] & posteriorData$alpha>alphaRange[1])&(posteriorData$beta<betaRange[2] & posteriorData$beta>betaRange[1]))
    }
  }

  plot=ggplot2::ggplot(posteriorData[idx,],ggplot2::aes(x=alpha,y=beta,fill=posterior))+ggplot2::geom_raster(interpolate = F)+viridis::scale_fill_viridis()+cowplot::theme_cowplot()
  return(plot)
})
