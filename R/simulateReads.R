methods::setGeneric("simulateReads",signature=c('object','expectedLibSize','replicates','times','errorModel'),
                    def = function(object,expectedLibSize=10^6,replicates=2,times=numeric(),errorModel=NULL) {
  standardGeneric("simulateReads")
})

#' Simulate Reads
#'
#' Simulate Reads from kineticModel object
#' @param object A kineticSim object
#' @param expectedLibSize The expected library sequencing depth. May vary from this based on sampling.
#' @param replicates Replicates per condition
#' @param times A vector of times
#' @param errorModel A function that takes in a vector of values and returns a vector of dispersions for a negative bionomial model
#' @name simulateReads
#' @include  class-kineticModel.R
#' @examples
#' ts=basicKineticModel(synthRate = 1:10,degRate = rep(0.3,10))
#' ts=calculateEquilibrium(ts)
#' @export
methods::setMethod("simulateReads", signature(object = "kineticModel"),function(object,expectedLibSize=10^6,replicates=2,times=numeric(),errorModel=NULL) {
  ## Check if an errorModel is needed. Then, if an error model is included, check for validity and update errorModel
  if(is.null(errorModel)){
    if(is.null(object@errorModel(1))){
      stop("There is no pre-specified kineticModel error model so an error model must be provided.")
    }
  } else if(is.function(errorModel)){
    if(length(errorModel(1:10))==10 && is.numeric(errorModel(1:10))){
      object@errorModel=errorModel
    } else {
      stop("errorModel function must produce a *numeric* vector of the same length as the input.")
    }
  } else {
    stop("errorModel must be function")
  }
  ## Check all the other arguments
  if(!is.numeric(expectedLibSize) && length(expectedLibSize)==1)
    stop("lib.size must be a numeric of length one")
  if(!is.numeric(replicates) || replicates < 1)
    stop("replicates must be an integer greater than 1")
  if(!is.numeric(times) || length(times)==0){
    if(length(object@times) < 1){
      stop("times must be a numeric vector with at least one element if it is not already specified in the kineticModel object.")
    }
  } else {
    object@times=times
  }
  ## Create a matrix of the appropriate size
  simReads=matrix(nrow=length(object@synthRates),ncol=length(object@times)*replicates)
  ## Get expected number of reads for each transcript
  object=predictAbundance(ts,object@times)
  ## Rescale for library size
  z=prop.table(object@simData,2)*expectedLibSize
  ind=1
  for(s in 1:ncol(z)){
    for(tr in 1:nrow(z)){
      simReads[tr,ind:(ind+replicates-1)]=rnbinom(n=replicates,size = object@errorModel(z[tr,s]),mu = z[tr,s])
    }
    ind=ind+replicates
  }
  rownames(simReads)=object@ids
  colnames(simReads)=paste0("time.",as.character(rep(object@times,each=replicates)))
  object@data=simReads
  object@expMetadata=data.frame(time=rep(object@times,each=replicates))
  object@sizeFactors=colSums(object@data)/max(colSums(object@data)) ## used trimmed mean as default
  return(object)
})
