#produces a negative log lilkelyhood function for a given (set of) gene(s)
#should be used for only one gene at a time for purposes of param optimization
#' @include getAbund.R
nllFactory = function(geneIdx,object,dispByGene)
{
  obs = object@data[geneIdx,]
  time=object@expMetadata$time
  initVal = object@initVals[geneIdx]
  normFactors = object@normFactors
  dispersion = object@dispersionModel

  if(dispByGene)
  {
    return(function(params)
    {
      expMu = getAbund(params[1],params[2],time,initVal)*normFactors
      tryCatch(dnbinom(obs,mu=expMu,size=dispersion(geneIdx), log = T),warning=function(w)
        cat("\nWarning: Poor data quality causing optimizer to step outside of bounds\nAttempted: alpha =",params[1],"\tbeta =",params[2],"\n"))
      logProb = dnbinom(obs,mu=expMu,size=dispersion(geneIdx), log = T)
      return(-sum(logProb,na.rm=T)) #currently generating NaN for t=0
    })
  }
  else
  {
    return(function(params)
    {
      expMu = getAbund(params[1],params[2],time,initVal)*normFactors
      tryCatch(dnbinom(obs,mu=expMu,size=dispersion(expMu), log = T),warning=function(w)
        cat("\nWarning: Poor data quality causing optimizer to step outside of bounds\nAttempted: alpha=",params[1],"\tbeta=",params[2],"\n"))
      logProb = dnbinom(obs,mu=expMu,size=dispersion(expMu), log = T)
      return(-sum(logProb,na.rm=T)) #currently generating NaN for t=0
    })
  }

}
