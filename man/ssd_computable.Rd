% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computable.R
\name{ssd_computable}
\alias{ssd_computable}
\alias{ssd_computable.tmbfit}
\alias{ssd_computable.fitdists}
\title{Is Computable Standard Errors}
\usage{
ssd_computable(x, ...)

\method{ssd_computable}{tmbfit}(x, ...)

\method{ssd_computable}{fitdists}(x, ...)
}
\arguments{
\item{x}{The object.}

\item{...}{Unused.}
}
\value{
A flag for each distribution indicating if all parameters have numerically computable standard errors.

A flag indicating if all parameters have numerically computable standard errors.

A logical vector for each distribution indicating if all parameters have numerically computable standard errors.
}
\description{
Generic function to test if all parameters have numerically computable standard errors.
}
\section{Methods (by class)}{
\itemize{
\item \code{ssd_computable(tmbfit)}: Is Computable Standard for tmbfit Object

\item \code{ssd_computable(fitdists)}: Is At Boundary for fitdists Object

}}
\examples{
fits <- ssd_fit_dists(ssddata::ccme_boron, 
  dists = c("lnorm", "lnorm_lnorm", "burrIII3")
)
ssd_computable(fits$lnorm)
ssd_computable(fits$lnorm_lnorm)
ssd_computable(fits$burrIII3)

fits <- ssd_fit_dists(ssddata::ccme_boron, 
  dists = c("lnorm", "lnorm_lnorm", "burrIII3")
)
ssd_computable(fits)

}
