setGeneric("estimateDispersions",function(object,...) standardGeneric("estimateDispersions"))

#'Estimate Dispersions Using DESeq2
#'
#'Estimates the dispersions for each gene using DESeq2. If \code{byGene = TRUE} then a vector containing a single dispersion
#'estimate for each gene will be returned. If \code{FALSE} then a general function will be returned that gives a dispersion estimate
#'for a given mean using the composite information from all genes, time points, and replicates.
#'
#'@param object A \linkS4class{basicKineticModel} object
#'@param byGene Boolean controlling fuction return value. See \code{description} or \code{value} for more information
#'
#'@return If \code{byGene = TRUE} then a function that takes a gene index and returns a single dispersion
#'estimate for that gene. If \code{FALSE} then a general function will be returned that gives a dispersion estimate
#'for a given read value (assumed to be the mean of the distribution) using the composite information from all genes, time points, and replicates.
#'
#'@name estimateDispersions
#'@export

setMethod("estimateDispersions", signature(object="basicKineticModel"), function(object,byGene=T)
{
  validObject(object)
  if(is.na(object@data[1]))
    stop("Error: No data detected.  Please run simulateReads first or set @data slot in the provided basicKineticModel object")

  data = object@data
  colnames(data) = 1:ncol(data)
  expMetadata = object@expMetadata
  expMetadata$time = factor(expMetadata$time)
  dds = DESeq2::DESeqDataSetFromMatrix(countData = data, colData = expMetadata,  design= ~ time)
  DESeq2::normalizationFactors(dds) = matrix(object@normFactors,nrow=nrow(data),ncol=ncol(data),byrow = T)
  dds = DESeq2::estimateDispersions(dds)
  disp = DESeq2::dispersions(dds)
  bob<<-dds@dispersionFunction

  if(byGene)
  {
    fn = function(disp)
      {
        return(function(x) 1/disp[x])
      }
    fn2 = fn(disp)
    return(fn2)
  }
  else
    return(function(x) 1/dds@dispersionFunction(x))
})
