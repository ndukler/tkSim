#' Build simulation object
#'
#' Contructs a simple kinetic simulation object from synthesis and degredation rates
#' @param syn.rate Synthesis rates. Must be greater than 0.
#' @param deg.rate Degredation rate. Must be between 0 and 1.
#' @param init.abundance Initial abundance at time 0.
#' @param ids A character string id.
#' @param data Experimental or simulated matrix of count data
#' @param col.info A data.frame where each row corresponds to the column in data with the same index and the columns are different properties
#' @name simpleKineticExperiment
#' @include class-simpleKineticExperiment.R
#' @export
simpleKineticExperiment <- function(syn.rate=NA_real_,deg.rate=NA_real_,init.abundance=NA_real_,ids=NA,data=NULL,col.info=data.frame(),error.model=function(x){}){
  if(length(init.abundance) == 1 && is.na(init.abundance))
    init.abundance=rep(0,length(syn.rate))
  if(length(ids) == 1 && is.na(ids))
    ids=as.character(1:max(length(syn.rate),nrow(data)))
  if(is.null(data))
    data=matrix(nrow=0,ncol=0)
  new("simpleKineticExperiment",synthesis.rates=syn.rate,degredation.rates=deg.rate,initial.values=init.abundance,ids=ids,col.info=col.info,
      data=data,error.model=error.model)
}
