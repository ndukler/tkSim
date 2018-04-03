setGeneric("simulateReads",signature=c('object'),def = function(object,...) {standardGeneric("simulateReads")})

#' Simulate Reads
#'
#' Simulate reads from a \linkS4class{kineticModel} object. Will simulate data if the input \linkS4class{kineticModel} does not already contain simulated data.
#' @param object A \linkS4class{kineticModel} object
#' @param expectedLibSize The expected total number of reads per sequencing run/batch. May be supplied as a single number which will be applied uniformly,
#' or a vector of library sizes where each element corresponds to the library depth at a given time point.
#' @param replicates Number of replicates per condition
#' @param numSpikeIns The number of unique spike in transcripts used.
#' @param spikeInSizes The expected number of reads for each type of spike in used. May be a single number used for all spike-in transcripts
#'  or an array of abundances for each unique spike-in transcript.
#' @param dispersionModel A function for handling dispersions for a negative bionomial model. If \code{dispByGene = TRUE} then the function should take a gene
#' index as an input and return a sigle dispersion estimate as an output.  If \code{dispByGene = FALSE} then the function should take a mean as an input and return
#' a dispersion estimate based on the value of that mean as an output. See the return values from \code{\link{estimateDispersions}} for examples.
#' @param dispByGene Boolean controlling the expected nature of the \code{dispersionModel}. See \code{dispersionModel} description for more details.
#'
#' @name simulateReads
#' @include  class-kineticModel.R
#' @examples
#' ##setup
#' bkm = basicKineticModel(times=0:30,synthRate = 1:10,degRate = rep(0.3,10))
#' bkm = simulateData(bkm) #optional
#'
#' ##mean based dispersion estimation, equal spikeIns
#' bkm = simulateReads(bkm,expectedLibSize=10^6,replicates=3,numSpikeIns=4,spikeInSizes=200,dispersionModel=function(x){rep(10^4,length(x))},dispByGene=F)
#'
#' ##unequal spike ins
#' bkm = simulateReads(bkm,expectedLibSize=10^6,replicates=3,numSpikeIns=4,spikeInSizes=c(100,200,100,300),dispersionModel=function(x){rep(10^4,length(x))},dispByGene=F)
#'
#' ##unequal library sizes
#' els = seq(10^6,10^8,length.out=31)
#' bkm = simulateReads(bkm,expectedLibSize=els,replicates=3,numSpikeIns=4,spikeInSizes=200,dispersionModel=function(x){rep(10^4,length(x))},dispByGene=F)
#'
#' @export
setMethod("simulateReads", signature(object = "kineticModel"),function(object,expectedLibSize=10^6,replicates=2,numSpikeIns=4,spikeInSizes=NULL,dispersionModel=NULL,dispByGene=F){
  validObject(object)
  ## Check if a dispersionModel is needed. Then, if an dispersion model is included, check for validity and update dispersionModel
  if(is.null(dispersionModel)){
    if(is.null(object@dispersionModel(1))){
      stop("There is no pre-specified kineticModel dispersion model so an dispersion model must be provided.")
    }
  } else if(is.function(dispersionModel)){
    if(length(dispersionModel(1:10))==10 && is.numeric(dispersionModel(1:10))){
      cat("\nUsing Provided Dispersion Model.  Replacing current dispersion model for returned object with user provided dispersion model.\n")
      object@dispersionModel=dispersionModel
    } else {
      stop("dispersionModel function must produce a *numeric* vector of the same length as the input.")
    }
  } else {
    stop("dispersionModel must be function")
  }
  ## Check all the other arguments
  if(!is.numeric(expectedLibSize))
    stop("expectedLibsize must be numeric")
  else if(length(expectedLibSize)!=1 && length(expectedLibSize)!=length(object@times))
    stop("Error: expectedLibSize must either be a single value or a vector of equal length to object@times")
  if(!is.numeric(replicates) || replicates < 1)
    stop("replicates must be an integer greater than 1")
  if(!is.numeric(numSpikeIns) || numSpikeIns < 1)
    stop("numSpikeIns must be an integer greater than 1")
  if(!is.numeric(spikeInSizes) || numSpikeIns < 1 || length(spikeInSizes)<numSpikeIns && length(spikeInSizes)!=1)
    stop("Must specify the size of spike-ins.  May be a single number (greater than 1) used for all spike-in transcripts or an array of abundances for each unique spike-in transcript.")
  if(length(object@times) < 1)
    stop("Times not specified. Time must be specified in the kineticModel object as a numeric vector with at least one element.")
  if(length(object@synthRates)==0)
      stop("Error: Synthesis or degredation rates not defined in kineticModel.")


  ## Create a matrix of the appropriate size
  simReads=matrix(0,nrow=length(object@synthRates),ncol=length(object@times)*replicates)
  spikeIns=matrix(0,ncol=length(object@times)*replicates,nrow=numSpikeIns)
  ## Get expected number of reads for each transcript
  if(length(object@simData) ==0)
  {
    cat("\nNo simulated data detected.  Simulating data before simulating reads.\n")
    object=simulateData(object)
  }

  temp=rbind(object@simData,matrix(spikeInSizes,nrow=numSpikeIns,ncol=ncol(object@simData)))
  ## Rescale for library size
  if(length(expectedLibSize)==1)
    z=prop.table(temp,margin = 2)*expectedLibSize #prop.table(margin=2) => column percentages
  else
    z=prop.table(temp,margin = 2)%*%diag(expectedLibSize) #multiply each col of Z by the coresponding element in the els vector

  if(dispByGene)
  {
    indx=1
    for(t in 1:ncol(z)){
      for(tr in 1:nrow(object@simData)){
        simReads[tr, indx:(indx+replicates-1)] = rnbinom(n=replicates,size = object@dispersionModel(tr),mu = z[tr,t])
      }
      spikeIndx=1
      for(tr in (nrow(object@simData)+1):nrow(z)){
        spikeIns[spikeIndx, indx:(indx+replicates-1)] = rnbinom(n=replicates,size = object@dispersionModel(tr),mu = z[tr,t]) #HANDLE SI differently for DISP???
        spikeIndx = spikeIndx+1
      }
      indx = indx+replicates
    }
  }

  if(!dispByGene)
  {
    indx=1
    for(t in 1:ncol(z)){
      for(tr in 1:nrow(object@simData)){
        simReads[tr, indx:(indx+replicates-1)] = rnbinom(n=replicates,size = object@dispersionModel(z[tr,t]),mu = z[tr,t])
      }
      spikeIndx=1
      for(tr in (nrow(object@simData)+1):nrow(z)){
        spikeIns[spikeIndx, indx:(indx+replicates-1)] = rnbinom(n=replicates,size = object@dispersionModel(z[tr,t]),mu = z[tr,t])
        spikeIndx = spikeIndx+1
      }
      indx = indx+replicates
    }
  }

  rownames(simReads) = object@ids
  colnames(simReads) = paste0("time_",as.character(rep(object@times,each=replicates)))
  rownames(simReads) = object@ids
  object@data = simReads
  object@expMetadata = data.frame(time=rep(object@times,each=replicates))

  ##normalize by first sample and then average per 'replicate'
  object@spikeIns = spikeIns
  object@normFactors = colMeans(spikeIns/spikeIns[,1])

  return(object)
})
