setGeneric("simulateReads",signature=c('object'),
                    def = function(object,...) {
  standardGeneric("simulateReads")})

#' Simulate Reads
#'
#' Simulate Reads from kineticModel object
#' @param object A kineticModel object
#' @param expectedLibSize The expected library sequencing depth. May vary from this based on sampling.
#' @param replicates Replicates per condition
#' @param times A vector of times
#' @param errorModel A function that takes in a vector of values and returns a vector of dispersions for a negative bionomial model
#' @name simulateReads
#' @include  class-kineticModel.R
#' @examples
#' bkm = basicKineticModel(synthRate = 1:10,degRate = rep(0.3,10))
#' bkm = simulateData(bkm) #optional
#' bkm = simulateReads(bkm,expectedLibSize=10^6,replicates=1,errorModel=function(x){rep(2,length(x))})
#' @export
setMethod("simulateReads", signature(object = "kineticModel"),function(object,expectedLibSize=10^6,replicates=2,numSpikeIns=4,spikeInSizes=numeric(),errorModel=NULL){
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

  if(length(spikeInSizes)==0)
    stop("Must specify the size of spike-ins.  May be a single number used for all spike-in transcripts or an array of abundances for each unique spike-in transcript.")

  ## Create a matrix of the appropriate size
  simReads=matrix(0,nrow=length(object@synthRates),ncol=length(object@times)*replicates)
  spikeIns=matrix(0,ncol=length(object@times)*replicates,nrow=numSpikeIns)
  ## Get expected number of reads for each transcript
  object=predictAbundance(object,object@times)
  temp=rbind(object@simData,matrix(spikeInSizes,nrow=numSpikeIns,ncol=ncol(object@simData)))
  ## Rescale for library size
  z=prop.table(temp,margin = 2)*expectedLibSize #prop.table(margin=2) => column percentages
  # print(z/expectedLibSize)
  indx=1
  for(t in 1:ncol(z)){
    for(tr in 1:nrow(object@simData)){
      simReads[tr, indx:(indx+replicates-1)] = rnbinom(n=replicates,size = object@errorModel(z[tr,t]),mu = z[tr,t])
    }
    spikeIndx=1
    for(tr in (nrow(object@simData)+1):nrow(z)){
      spikeIns[spikeIndx, indx:(indx+replicates-1)] = rnbinom(n=replicates,size = object@errorModel(z[tr,t]),mu = z[tr,t])
      spikeIndx = spikeIndx+1
    }
    indx = indx+replicates
  }
  rownames(simReads) = object@ids
  colnames(simReads) = paste0("time_",as.character(rep(object@times,each=replicates)))
  rownames(simReads) = object@ids
  object@data = simReads
  object@expMetadata = data.frame(time=rep(object@times,each=replicates))
  #OLD object@sizeFactors=(colSums(object@data)+colSums(spikeIns))/colSums(object@data) ## used trimmed mean as default,  FINISH
  ##normalize by first sample and then average per 'replicate'
  object@spikeIns = spikeIns
  object@sizeFactors = colMeans(spikeIns/spikeIns[,1])

  return(object)
})
