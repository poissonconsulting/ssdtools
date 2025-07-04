% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcanz.R
\name{ssd_hp_bcanz}
\alias{ssd_hp_bcanz}
\title{BCANZ Hazard Proportion}
\usage{
ssd_hp_bcanz(
  x,
  conc = 1,
  ...,
  nboot = 10000,
  min_pboot = 0.95,
  proportion = FALSE
)
}
\arguments{
\item{x}{The object.}

\item{conc}{A numeric vector of concentrations to calculate the hazard proportions for.}

\item{...}{Unused.}

\item{nboot}{A count of the number of bootstrap samples to use to estimate the confidence limits. A value of 10,000 is recommended for official guidelines.}

\item{min_pboot}{A number between 0 and 1 of the minimum
proportion of bootstrap samples that must successfully fit (return a likelihood)
to report the confidence intervals.}

\item{proportion}{A numeric vector of proportion values to estimate hazard concentrations for.}
}
\value{
A tibble of corresponding hazard concentrations.
}
\description{
Gets  proportion of species affected at specified concentration(s)
using settings adopted by BC, Canada, Australia and New Zealand for official guidelines.
This function can take several minutes to run with recommended 10,000 iterations.
}
\examples{
fits <- ssd_fit_bcanz(ssddata::ccme_boron)
ssd_hp_bcanz(fits, nboot = 100)
}
\seealso{
\code{\link[=ssd_hp]{ssd_hp()}}.

Other BCANZ: 
\code{\link{ssd_fit_bcanz}()},
\code{\link{ssd_hc_bcanz}()}
}
\concept{BCANZ}
