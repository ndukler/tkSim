setGeneric("plotParameterFit", function(object,...) standardGeneric("plotParameterFit"))

#' @name plotParameterFit
#' @include  class-kineticModel.R
#' @export
setMethod("plotParameterFit",signature(object="basicKineticModel"), function(object,geneIdx=NULL,legend=F){
  if(is.null(geneIdx))
    geneIdx = nrow(object@data)

  getAbund = function(alpha,beta,time,initVal)
  {
    return(exp(-time * beta) * (initVal - alpha / beta) + alpha / beta)
  }

  fitData = mapply(getAbund,alpha=object@inferedParams[geneIdx,1],beta=object@inferedParams[geneIdx,2],initVal=object@initVals[geneIdx],MoreArgs = list(time=object@times))

  colnames(fitData) = object@ids[geneIdx]
  #compute average normalization factors per time point
  times = object@expMetadata$time
  normFactor = sapply(unique(times),function(x,times,sizeFactors){mean(sizeFactors[which(times==x)])},times=times,sizeFactors=object@sizeFactors)
  fitData = fitData*normFactor
  rownames(fitData) = object@times
  fitData = reshape2::melt(fitData)
  colnames(fitData) = c("time","gene","value")

  plotData = reshape2::melt(object@data[geneIdx,])
  plot = ggplot2::ggplot()
  plot = plot+ggplot2::geom_line(data=fitData,aes(x=time,y=value,group=factor(gene)))+
    ggplot2::geom_point(data=plotData,aes(x=as.numeric(gsub("time_","",Var2)),y=value,color=as.factor(Var1),group=as.factor(Var1)))+
    cowplot::theme_cowplot()+
    ggplot2::xlab("Time")+
    ggplot2::ylab("Labeled Transcripts")
  if(legend)
  {
    plot=plot+ggplot2::labs(color="Transcript Ids")
  }else{
    plot=plot+ggplot2::guides(color=FALSE)
  }
  print(plot)
})

