% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/afs.r
\name{shiftFromZero}
\alias{shiftFromZero}
\title{Shift from zero}
\source{
\url{https://github.com/alexjgriffith/mulcal}
}
\usage{
shiftFromZero(summit)
}
\arguments{
\item{value}{integer or double which cannot be less than 0}

\item{width}{integer or double}
}
\value{
a range value+width value-width where value-width >0
}
\description{
Checks if the result of subtracting width from value is zero
}
\author{
Alexander Griffith griffita@gmail.com
}

