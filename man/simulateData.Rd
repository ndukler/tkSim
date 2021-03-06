% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulateData.R
\docType{methods}
\name{simulateData}
\alias{simulateData}
\title{Simulate Theoretical Data}
\usage{
\S4method{simulateData}{basicKineticModel}(object)
}
\arguments{
\item{object}{A \linkS4class{basicKineticModel} object}
}
\value{
Returns a \linkS4class{basicKineticModel} object containing the simulated data with any parameters
overwritten by those supplied by the user.
}
\description{
Simulates theoretical data based on given time points, synthesis rate, and degredation rate.
}
\examples{
bkm=basicKineticModel(synthRate = 1:10,degRate = rep(0.3,10), times=0:30)
bkm=simulateData(bkm)
}
