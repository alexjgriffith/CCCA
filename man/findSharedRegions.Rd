% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motifComparision.r
\name{findSharedRegions}
\alias{findSharedRegions}
\title{Find Shared Regions}
\usage{
findSharedRegions(lm, lc, n1, n2, reg = NULL)
}
\description{
Finds the intersection of the motifs and optionaly a subset region
this function is meant to be called after the locations have been
found (i.e. using grep("motif",Sequence,ignore.case=TRUE)
}

