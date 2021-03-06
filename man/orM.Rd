% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{orM}
\alias{orM}
\title{orM}
\usage{
orM(mat)
}
\arguments{
\item{mat}{input matrix of logicals}
}
\value{
a vector of logicals
}
\description{
applies the or (|) function accross a matrix recursivly
}
\examples{
a<-as.matrix(cbind(c(TRUE,TRUE,TRUE),
                    c(FALSE,TRUE,FALSE)))
result<-orM(a)
result==c(TRUE,TRUE,TRUE)
}
\seealso{
~\code{\link{ors} \link{"|"}}
}

