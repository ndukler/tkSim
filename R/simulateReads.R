methods::setGeneric("simulateReads",signature=c('object','expected.lib.size','replicates','times','error.model'),
                    def = function(object,expected.lib.size=10^6,replicates=2,times=numeric(),error.model=NULL) {
  standardGeneric("simulateReads")
})

#' Simulate Reads
#'
#' Simulate Reads from kineticExperiment object
#' @param object A kineticSim object
#' @param expected.lib.size The expected library sequencing depth. May vary from this based on sampling.
#' @param replicates Replicates per condition
#' @param times A vector of times
#' @param error.model A function that takes in a vector of values and returns a vector of dispersions for a negative bionomial model
#' @name simulateReads
#' @include  class-kineticExperiment.R
#' @examples
#' ts=simpleKineticExperiment(syn.rate = 1:10,deg.rate = rep(0.3,10))
#' ts=calculateEquilibrium(ts)
#' @export
methods::setMethod("simulateReads", signature(object = "kineticExperiment"),function(object,expected.lib.size=10^6,replicates=2,times=numeric(),error.model=NULL) {
  ## Check if an error.model is needed. Then, if an error model is included, check for validity and update error.model
  if(is.null(error.model)){
    if(is.null(object@error.model(1))){
      stop("There is no pre-specified kineticExperiment error model so an error model must be provided.")
    }
  } else if(is.function(error.model)){
    if(length(error.model(1:10))==10 && is.numeric(error.model(1:10))){
      object@error.model=error.model
    } else {
      stop("error.model function must produce a *numeric* vector of the same length as the input.")
    }
  } else {
    stop("error.model must be function")
  }
  ## Check all the other arguments
  if(!is.numeric(expected.lib.size) && length(expected.lib.size)==1)
    stop("lib.size must be a numeric of length one")
  if(!is.numeric(replicates) || replicates < 1)
    stop("replicates must be an integer greater than 1")
  if(!is.numeric(times) || length(times)==0){
    if(length(object@times) < 1){
      stop("times must be a numeric vector with at least one element if it is not already specified in the kineticExperiment object.")
    }
  } else {
    object@times=times
  }
  ## Create a matrix of the appropriate size
  sim.dat=matrix(nrow=length(object@synthesis.rates),ncol=length(object@times)*replicates)
  ## Get expected number of reads for each transcript
  object=predictAbundance(ts,object@times)
  ## Rescale for library size
  z=prop.table(object@predicted.abundance,2)*expected.lib.size
  ind=1
  for(s in 1:ncol(z)){
    for(tr in 1:nrow(z)){
      sim.dat[tr,ind:(ind+replicates-1)]=rnbinom(n=replicates,size = object@error.model(z[tr,s]),mu = z[tr,s])
    }
    ind=ind+replicates
  }
  rownames(sim.dat)=object@ids
  colnames(sim.dat)=paste0("time.",as.character(rep(object@times,each=replicates)))
  object@data=sim.dat
  object@col.info=data.frame(time=rep(object@times,each=replicates))
  object@size.factors=colSums(object@data)/max(colSums(object@data)) ## used trimmed mean as default
  return(object)
})

