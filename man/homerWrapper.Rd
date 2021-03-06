% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/homerWrapper.r
\name{homerWrapper}
\alias{homerWrapper}
\title{Homer Wrapper}
\usage{
homerWrapper(sequences, foreground, background, homerLocation,
  motifsFile = FALSE, opts = "-S 25 -len 6")
}
\arguments{
\item{sequences}{Biostring dnastring object}

\item{foreground}{The foreground subset of sequences, list}

\item{background}{The background subset of sequences, list}

\item{homerLocation}{The locatoin of the homer executable on your system}

\item{motifsFile}{The file to save the motifs in,
if none then the motifs are note saved}

\item{opts}{A string of additional parameters to be passed to homer}
}
\value{
a pwm of the motifs found.
}
\description{
This is only a wrapper for the DENOVO functionality of homer2, homer must be installed on the system before this can be used.
}

