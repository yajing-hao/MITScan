% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MITScan_ZIP.R
\name{scorep.zip}
\alias{scorep.zip}
\title{Compute score test statistics for a given OTU across host genes using zero-inflated Poisson models}
\usage{
scorep.zip(Y, X, Z, A, Ap)
}
\arguments{
\item{Y}{A vector of abundances for a given OTU.}

\item{X}{A matrix of covariates related to the count component.}

\item{Z}{A matrix of covariates related to the zero component.}

\item{A}{A matrix of host gene expression from RNA-Seq.
Rows and columns represent host gene names and samples, respectively.}

\item{Ap}{A matrix of permuted host gene expression.}
}
\value{
A vector of score test statistics for a given OTU across all host genes
and the median of permuted score test statistics.
}
\description{
scorep.zip() computes score test statistics for a given OTU
across all host genes using zero-inflated Poisson models,
and return a vector where the i-th element of the vector represents
the score test statistic for the i-th host gene and the given OTU.
}
