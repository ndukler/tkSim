#' Calculate Posterior Probabilities for Infered Parameters
#'
#' Uses numeric methods to estimate posteriors for the infered parameters \code{alpha} (synthesis rate) and \code{beta} (degredation rate).
#' Currently uses a flat prior for both alpha and beta as a default. This function is written such that non-flat priors may be used in the future.
#' @param object A \linkS4class{basicKineticModel} object
#' @param alphaRange Scale factors used to calculate the upper and lower bounds of the parameter range explored for \code{alpha}.  These scale factors will
#' be applied to the infered value of \code{alpha}.  Must be defined as \code{c(lower,upper)}.
#' @param paramSpaceSize The total size of parameter space to numerically integrate over. Half of the parameter space will be given to each parameter.
#' @param logProbAlpha  A function that returns the log probability for a given value of \code{alpha}
#' @param lobProbBeta A function that returns the log probability for a given value of \code{beta}
#'
#' @name calculatePosteriors
#' @include  class-basicKineticModel.R
#' @examples
#' EXAMPLE HERE
#' @export
calculatePosteriors = function(object,alphaRange=numeric(2),paramSpaceSize=10^4,logProbAlpha=NULL,logProbBeta=NULL)
{
  logSumExp = function(x)
  {
    a = max(x)
    return(a+log(sum(exp(x-a))))
  }

  logLFactory = function(geneIdx,object)
  {
    obs = object@data[geneIdx,] #-1 to remove NaN at 0 from read simulation function
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
      return(sum(logProb,na.rm=T))
    })
  }

  if(is.null(logProbAlpha))
  {
    # #defined as 1/(max(a)-min(a)) on a per-gene basis
    # logProbAlpha = function(geneIdx,object,alphaRange)
    # {
    #   aMax = alphaRange[2]*object@inferedParams[geneIdx,"alpha"]
    #   aMin = alphaRange[1]*object@inferedParams[geneIdx,"alpha"]
    #   return(function(x){log(1/(aMax-aMin))})
    # }
    # logProbAlpha = lapply(X=1:nrow(object@data),FUN=logProbAlpha,object=object,alphaRange=alphaRange)

    #global basis
    aMax = alphaRange[2]*max(object@inferedParams[,'alpha'])
    aMin = alphaRange[1]*min(object@inferedParams[,'alpha'])
    logProbAlpha = function(x){rep(log(1/sqrt(paramSpaceSize)),length(x))} #log(1/(aMax-aMin)) == -log(aMax-aMin)
  }

  if(is.null(logProbBeta))
  {
    logProbBeta = function(x){rep(log(1),length(x))}
  }

  #generate likelyhood esitmators for each gene
  logLH = lapply(X=1:nrow(object@inferedParams), FUN=logLFactory,object=object)
  posteriors = lapply(X=1:nrow(object@inferedParams),object=object,logLH=logLH,FUN=function(x,object,logLH)
    {
      alpha = object@inferedParams[x,"alpha"]
      beta = object@inferedParams[x,"beta"]
      aMax = alphaRange[2]*alpha
      aMin = alphaRange[1]*beta

      # numerator = logLH[[x]](c(alpha,beta))+logProbAlpha(alpha)+logProbBeta(beta)

      paramRange = expand.grid(seq(aMin,aMax,length.out = sqrt(paramSpaceSize)), seq(10^-5,1,length.out = sqrt(paramSpaceSize)))
      numerator = apply(paramRange,1,function(y) logLH[[x]](y)) + logProbAlpha(paramRange[,1]) + logProbBeta(paramRange[,2])
      marginal = logSumExp(numerator)
      posterior = exp(numerator-marginal)
      foo=cbind(paramRange,posterior=posterior)
      library(ggplot2)
      library(viridis)
      ggplot(foo,aes(x=Var1,y=Var2,fill=posterior))+geom_raster(interpolate = T)+scale_fill_viridis()+cowplot::theme_cowplot()
      ##FINISH
      # object=bkm
      # paramSpaceSize=10^4
      # logProbAlpha=NULL
      # logProbBeta=NULL
      # alphaRange = c(.25,2)
      # >data.table::as.data.table(foo)[,sum(posterior),by=Var2]
      # > plot(data.table::as.data.table(foo)[,sum(posterior),by=Var2]$V1)
      # > plot(data.table::as.data.table(foo)[,sum(posterior),by=Var1]$V1)
      # > plot(data.table::as.data.table(foo)[,sum(posterior),by=Var2]$V1)
    })

}

