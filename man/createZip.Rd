% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{createZip}
\alias{createZip}
\title{Create Zip}
\usage{
createZip(x)
}
\arguments{
\item{x}{even length input vector}
}
\value{
a list of pairs
}
\description{
takes a vector and retuns a list of pairs
}
\examples{
createZip(c(1,2,3,4,5,6)) == list(c(1,2),c(3,4),c(5,6))
}

