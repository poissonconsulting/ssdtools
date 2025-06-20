% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gompertz.R
\name{qgompertz}
\alias{qgompertz}
\title{Quantile Function for Gompertz Distribution
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
qgompertz(p, llocation = 0, lshape = 0, lower.tail = TRUE, log.p = FALSE)
}
\arguments{
\item{p}{vector of probabilities.}

\item{llocation}{location parameter on the log scale.}

\item{lshape}{shape parameter on the log scale.}

\item{lower.tail}{logical; if TRUE (default), probabilities are \code{P[X <= x]}, otherwise, \code{P[X > x]}.}

\item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
}
\description{
Deprecated for \code{ssd_qgompertz()}.
}
\keyword{internal}
