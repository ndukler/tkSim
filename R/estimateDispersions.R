setGeneric("estimateDispersions",function(object,byGene) standardGeneric("estimateDispersions"))

#'Estimate Dispersions Using DESeq2
#'
#'Estimates the dispersions for each gene using DESeq2. If \code{byGene = TRUE} then a vector containing a single dispersion
#'estimate for each gene will be returned. If \code{FALSE} then a general function will be returned that gives a dispersion estimate
#'for a given mean using the composite information from all genes, time points, and replicates.
#'
#'@param object A \linkS4class{basicKineticModel} object
#'@param byGene Boolean controlling fuction return value. See \code{description} or \code{value} for more information
#'
#'@return If \code{byGene = TRUE} then a vector containing a single dispersion
#'estimate for each gene will be returned. If \code{FALSE} then a general function will be returned that gives a dispersion estimate
#'for a given mean using the composite information from all genes, time points, and replicates.
#'
#'@name estimateDispersions
#'@export

setMethod("estimateDispersions", signature(object="basicKineticModel"), function(object,byGene=T)
{
  if(is.na(object@data[1]))
    stop("Error: No data detected.  Please run simulateReads first or set @data slot in the provided basicKineticModel object")

  data = object@data
  colnames(data) = 1:ncol(data)
  dds = DESeq2::DESeqDataSetFromMatrix(countData = data, colData = factor(object@expMetadata),  design= ~ time)
  DESeq2::normalizationFactors(dds) = matrix(object@normFactors,nrow=nrow(data),ncol=ncol(data),byrow = T)
  dds = DESeq2::estimateDispersions(dds)
  disp = DESeq2::dispersions(dds)

  if(byGene)
    return(disp)
  else
    return(dds@dispersionFunction)
})
