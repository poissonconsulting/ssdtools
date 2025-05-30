% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot.R
\name{scale_colour_ssd}
\alias{scale_colour_ssd}
\alias{scale_color_ssd}
\alias{scale_fill_ssd}
\title{Discrete color-blind scale for SSD Plots}
\usage{
scale_colour_ssd(...)

scale_color_ssd(...)

scale_fill_ssd(...)
}
\arguments{
\item{...}{Arguments passed to \code{\link[ggplot2:discrete_scale]{ggplot2::discrete_scale()}}.}
}
\description{
The functions were designed for coloring different groups in a plot of SSD data.
}
\section{Functions}{
\itemize{
\item \code{scale_color_ssd()}: Discrete color-blind scale for SSD Plots

\item \code{scale_fill_ssd()}: Discrete color-blind scale for SSD Plots

}}
\examples{
# Use the color-blind palette for a SSD plot
ssd_plot(ssddata::ccme_boron, boron_pred, shape = "Group", color = "Group") +
  scale_colour_ssd()
# Use the color-blind palette for a histogram of concentrations
ggplot2::ggplot(ssddata::ccme_boron, ggplot2::aes(x = Species, y = Conc, fill = Group)) +
  ggplot2::geom_col() +
  scale_fill_ssd() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
}
\seealso{
Other ggplot: 
\code{\link{geom_hcintersect}()},
\code{\link{geom_ssdpoint}()},
\code{\link{geom_ssdsegment}()},
\code{\link{geom_xribbon}()},
\code{\link{ssd_pal}()}
}
\concept{ggplot}
