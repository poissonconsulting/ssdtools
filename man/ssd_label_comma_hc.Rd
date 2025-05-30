% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scales.R
\name{ssd_label_comma_hc}
\alias{ssd_label_comma_hc}
\title{Label numbers with significant digits and comma.
If \code{hc_value} is present in breaks, put on new line and make bold.}
\usage{
ssd_label_comma_hc(hc_value, digits = 3, ..., big.mark = ",")
}
\arguments{
\item{hc_value}{A number of the hazard concentration value to offset.}

\item{digits}{A whole number specifying the number of significant figures.}

\item{...}{Unused.}

\item{big.mark}{A string specifying used between every 3 digits to separate thousands on the x-axis.}
}
\value{
A "labelling" function that takes a vector x and
returns a character vector of \code{length(x)} giving a label for each input value.
}
\description{
Label numbers with significant digits and comma.
If \code{hc_value} is present in breaks, put on new line and make bold.
}
\examples{
ggplot2::ggplot(data = ssddata::anon_e, ggplot2::aes(x = Conc / 10)) +
  geom_ssdpoint() +
  ggplot2::scale_x_log10(labels = ssd_label_comma_hc(1.26))
}
\seealso{
\code{\link[scales:label_number]{scales::label_comma()}}
}
