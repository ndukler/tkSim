#' Build simulation object
#'
#' Contructs a simple kinetic simulation object from synthesis and degredation rates
#' @param syn.rate Synthesis rates. Must be greater than 0.
#' @param deg.rate Degredation rate. Must be between 0 and 1.
#' @param init.abundance Initial abundance at time 0.
#' @param ids A character string id.
#' @name simpleKineticExperiment
#' @include class-simpleKineticExperiment.R
#' @export
simpleKineticExperiment <- function(syn.rate,deg.rate,init.abundance=NA,ids=NA){
  if(length(init.abundance) == 1 && is.na(init.abundance))
    init.abundance=rep(0,length(syn.rate))
  if(length(ids) == 1 && is.na(ids))
    ids=as.character(1:length(syn.rate))
  new("simpleKineticSim",synthesis.rates=syn.rate,degredation.rates=deg.rate,initial.values=init.abundance,ids=ids)
}
