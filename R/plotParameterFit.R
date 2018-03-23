setGeneric("plotParameterFit", function(object,...) standardGeneric("plotParameterFit"))

#' @name plotParameterFit
#' @param scaleData If true, will scale the data by it's corresponding normalization factors. If false, will scale the fit function by the average
#' normalization factor across replicates for each time point.
#' @include  class-kineticModel.R getAbund.R
#' @export
setMethod("plotParameterFit",signature(object="basicKineticModel"), function(object,geneIdx=NULL,legend=F,scaleData=T){
  if(is.null(geneIdx))
    geneIdx = nrow(object@data)

  fitData = mapply(getAbund,alpha=object@inferedParams[geneIdx,"alpha"],beta=object@inferedParams[geneIdx,"beta"],initVal=object@initVals[geneIdx],MoreArgs = list(time=object@times))

  colnames(fitData) = object@ids[geneIdx]
  #compute average normalization factors per time point
  times = object@expMetadata$time
  normFactor = sapply(unique(times),function(x,times,normFactors){mean(normFactors[which(times==x)])},times=times,normFactors=object@normFactors)
  # fitData = fitData*normFactor
  rownames(fitData) = object@times
  # fitData = fitData*object@normFactors
  fitData = reshape2::melt(fitData)
  colnames(fitData) = c("time","gene","value")

  plotData = reshape2::melt(object@data[geneIdx,]/object@normFactors)
  print(head(plotData))
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

