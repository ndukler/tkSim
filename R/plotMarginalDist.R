#' Plot the marginal posterior distribution of an estimated parameter
#'
#' Plots the marginal distribution of a selected parameter for a given gene
#' @name plotMarginalDist
#' @param object A \linkS4class{basicKineticModel} object
#' @param geneIdx The row index in \code{@@data} of the gene from which you want to plot the marginal posterior distribution from
#' @param param A string of the name of the parameter you want to plot.
#' @export
plotMarginalDist = function(object,geneIdx,param)
{
  # > plot(data.table::as.data.table(bkm@posteriors[[1]])[,sum(posterior),by=alpha]$V1)
  posteriorData = object@posteriors[[geneIdx]]
  vals = unique(posteriorData[,param])
  margin=vapply(vals, FUN.VALUE = numeric(1), posteriorData=posteriorData, param=param, FUN = function(x,posteriorData,param)
  {
    sum(posteriorData[which(posteriorData[,param]%in%x),"posterior"])
  })
  plot(vals,margin)
}

