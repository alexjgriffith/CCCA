% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcaUtils.r
\name{qstem}
\alias{qstem}
\title{Quick Stem Plot}
\usage{
qstem(y, xlim = c(-32, 32), ...)
}
\arguments{
\item{y}{all locations before determining pileup}

\item{xlim}{This range that is to be ploted}

\item{...}{remaining argument to stem and plot}
}
\description{
Abbreviates stem plotting for the analysis of composite motifs
}
\examples{
x<-runif(1000)*100
qstem(x)
}

