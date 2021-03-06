methods::setGeneric("plotAbundances", function(object) {
  standardGeneric("plotAbundances")
})

#' Plot transcript abundance
#'
#' Plots transcript abundance. Must run calculateAbundance first.
#' @param object A basicKineticModel object
#' @name plotAbundances
#' @include  class-basicKineticModel.R
#' @examples
#' ts=basicKineticModel(synthRate = 1:10,degRate = rep(0.3,10))
#' ts=calculateAbundance(ts,0:30)
#' plotAbundances(ts)
#' @export
methods::setMethod("plotAbundances", signature(object = "kineticModel"), function(object) {
  if(ncol(object@simData)>0){
    ab.m=data.table::as.data.table(reshape2::melt(object@simData,value.name="Abundance"))
    ab.m[,Var2:=object@times[Var2]]
    data.table::setnames(ab.m,c("Var1","Var2"),c("Id","Time"))
    g <- ggplot2::ggplot(ab.m,ggplot2::aes(x=Time,y=Abundance,group=Id))+
      ggplot2::geom_line()+
      ggplot2::theme_classic()
    return(g)
  } else {
    stop("Must run simulateData first.")
  }
})
