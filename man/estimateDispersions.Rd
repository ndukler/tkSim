% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateDispersions.R
\docType{methods}
\name{estimateDispersions}
\alias{estimateDispersions}
\title{Estimate Dispersions Using DESeq2}
\usage{
\S4method{estimateDispersions}{basicKineticModel}(object, byGene = T)
}
\arguments{
\item{object}{A \linkS4class{basicKineticModel} object}

\item{byGene}{Boolean controlling fuction return value. See \code{description} or \code{value} for more information}
}
\value{
If \code{byGene = TRUE} then a function that takes a gene index and returns a single dispersion
estimate for that gene. If \code{FALSE} then a general function will be returned that gives a dispersion estimate
for a given read value (assumed to be the mean of the distribution) using the composite information from all genes, time points, and replicates.
}
\description{
Estimates the dispersions for each gene using DESeq2. Overwrites the existing dipersion model in the supplied \code{\linkS4class{basicKinetic Model}}
and then returns the object. If \code{byGene = TRUE} then a function returning  a single dispersion
estimate for each gene will be set to \code{@dispersionModel}. If \code{FALSE} then a general function that gives a dispersion estimate
for a given mean using the composite information from all genes, time points, and replicates will be set to \code{@dispersionModel}.
}
