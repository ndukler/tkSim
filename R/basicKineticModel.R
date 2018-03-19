#' Basic Kinetic Model Constructor
#'
#' Contructs a simple kinetic model object from synthesis and degredation rates
#' @param synthRate Synthesis rates. Must be greater than 0.
#' @param degRate Degredation rate. Must be between 0 and 1.
#' @param initAbund Initial abundance at time 0.
#' @param ids A character string id.
#' @param data Experimental or simulated matrix of count data
#' @param expMetadata A data.frame where each row corresponds to the column in data with the same index and the columns are different properties
#' @name basicKineticModel
#' @include class-basicKineticModel.R
#' @export
basicKineticModel <- function(times=NA_real_,synthRate=NA_real_,degRate=NA_real_,initAbund=NA_real_,ids=NA,data=NULL,expMetadata=data.frame(),dispersionModel=function(x){} ){
  if(length(initAbund) == 1 && is.na(initAbund))
    initAbund=rep(0,length(synthRate))
  if(length(ids) == 1 && is.na(ids))
    ids=as.character(1:max(length(synthRate),nrow(data)))
  if(is.null(data))
    data=matrix(nrow=0,ncol=0)
  new("basicKineticModel",times=times,synthRates=synthRate,degRates=degRate,initVals=initAbund,ids=ids,expMetadata=expMetadata,
      data=data,dispersionModel=dispersionModel)
}
