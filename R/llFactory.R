#produces a log lilkelyhood function for a given (set of) gene(s)
#should be used for only one gene at a time for purposes of posterior estimation
#' @include getAbund.R
llFactory = function(geneIdx,object,dispByGene)
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
      logProb = dnbinom(obs,mu=expMu,size=dispersion(geneIdx), log = T)
      return(sum(logProb,na.rm=T))
    })
  }else
  {
    return(function(params)
    {
      expMu = getAbund(params[1],params[2],time,initVal)*normFactors
      logProb = dnbinom(obs,mu=expMu,size=dispersion(expMu), log = T)
      return(sum(logProb,na.rm=T))
    })
  }


}
