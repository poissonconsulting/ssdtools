% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcanz.R
\name{ssd_hc_bcanz}
\alias{ssd_hc_bcanz}
\title{BCANZ Hazard Concentrations}
\usage{
ssd_hc_bcanz(x, ..., nboot = 10000, min_pboot = 0.95)
}
\arguments{
\item{x}{The object.}

\item{...}{Unused.}

\item{nboot}{A count of the number of bootstrap samples to use to estimate the confidence limits. A value of 10,000 is recommended for official guidelines.}

\item{min_pboot}{A number between 0 and 1 of the minimum
proportion of bootstrap samples that must successfully fit (return a likelihood)
to report the confidence intervals.}
}
\value{
A tibble of corresponding hazard concentrations.
}
\description{
Gets hazard concentrations with confidence intervals that protect
1, 5, 10 and 20\% of species using settings adopted by
BC, Canada, Australia and New Zealand for official guidelines.
This function can take several minutes to run with recommended 10,000 iterations.
}
\examples{
fits <- ssd_fit_bcanz(ssddata::ccme_boron)
ssd_hc_bcanz(fits, nboot = 100)
}
\seealso{
\code{\link[=ssd_hc]{ssd_hc()}}.

Other BCANZ: 
\code{\link{ssd_fit_bcanz}()},
\code{\link{ssd_hp_bcanz}()}
}
\concept{BCANZ}
