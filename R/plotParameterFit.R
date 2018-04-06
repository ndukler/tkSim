setGeneric("plotParameterFit", function(object,...) standardGeneric("plotParameterFit"))

#' Plot Data Against the Fit Parameters
#'
#' Plot the curve defined by the estimated parameters against the raw data for a selected set of genes
#' @name plotParameterFit
#'
#' @param object A \linkS4class{basicKineticModel} object
#' @param geneIdx The row index from \code{@@data} of the gene to plot. May be defined as an integer for one gene or a vector for multiple genes.
#' @param legend Boolean controlling if the plot should contain a legend.
#' @param scaleData If true, will scale the data by a given genes's corresponding normalization factors. If false, will scale the fit function by the average
#' normalization factor across replicates for each time point.
#' @param plotTimes A vector of times specifying the times over which to plot the data and parameter curve. All time points in \code{object@@data} that are also in
#' \code{plotTimes} will be plotted. Time points \strong{only} in \code{plotTimes} will be ignored.  If \code{NULL} the times specified by \code{object@@times}
#' will be used.
#' @include  class-kineticModel.R getAbund.R
#' @export
setMethod("plotParameterFit",signature(object="basicKineticModel"), function(object,geneIdx=NULL,legend=F,scaleData=T,plotTimes=NULL){
  validObject(object)

  if(is.null(geneIdx))
    geneIdx = nrow(object@data)
  times = object@expMetadata$time

  fitData = mapply(getAbund,alpha=object@inferedParams[geneIdx,"alpha"],beta=object@inferedParams[geneIdx,"beta"],initVal=object@initVals[geneIdx],MoreArgs = list(time=unique(times)))

  colnames(fitData) = object@ids[geneIdx]
  #compute average normalization factors per time point

  if(!scaleData)
  {
    fitNormFactors = sapply(unique(times),function(x,times,normFactors){mean(normFactors[which(times==x)])},times=times,normFactors=object@normFactors)
    fitData = fitData*fitNormFactors
  }
  rownames(fitData) = unique(times)
  fitData = reshape2::melt(fitData)
  colnames(fitData) = c("time","gene","value")

  if(scaleData)
  {
    plotData = reshape2::melt(object@data[geneIdx,]/object@normFactors)
  }else{
    plotData = reshape2::melt(object@data[geneIdx,])
  }

  #plot only the time points specified in plotTimes
  if(!is.null(plotTimes))
  {
    plotData = plotData[which(as.numeric(gsub("time_","",plotData$Var2))%in%plotTimes),]
    fitData = fitData[which(fitData$time%in%plotTimes),]
  }

  #plot data and fit
  plot = ggplot2::ggplot()
  plot = plot+ggplot2::geom_line(data=fitData,ggplot2::aes(x=time,y=value,group=factor(gene)))
  if(ncol(plotData)==1)
  {
    plot=plot+ggplot2::geom_point(data=plotData,ggplot2::aes(x=as.numeric(gsub("time_","",colnames(object@data))),y=value),colour="blue")
  }else{
    plot=plot+ggplot2::geom_point(data=plotData,ggplot2::aes(x=as.numeric(gsub("time_","",Var2)),y=value,color=as.factor(Var1),group=as.factor(Var1)))
  }
  plot=plot+cowplot::theme_cowplot()+
    ggplot2::xlab("Time")+
    ggplot2::ylab("Labeled Transcripts")
  if(legend)
  {
    plot=plot+ggplot2::labs(color="Transcript Ids")
  }else{
    plot=plot+ggplot2::guides(color=FALSE)
  }
  return(plot)
})

