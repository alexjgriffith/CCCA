% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcaUtils.r
\name{stringToSwap}
\alias{stringToSwap}
\title{string to swap fun}
\usage{
stringToSwap(x)
}
\arguments{
\item{x}{input string}
}
\value{
a function mapping the even values of x to the odd
}
\description{
takes a string and produces a swap fun from it. The string is
broken allong spaces.
}
\examples{
stringToSwap("a b c d")("b")=="a"
}

