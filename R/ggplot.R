# Copyright 2015-2023 Province of British Columbia
# Copyright 2021 Environment and Climate Change Canada
# Copyright 2023-2024 Australian Government Department of Climate Change,
# Energy, the Environment and Water
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#       https://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

#' Color-blind Palette for SSD Plots
#'
#' @return A character vector of a color blind palette with 8 colors.
#' @family ggplot
#' @export
#' @examples
#' ssd_pal()
ssd_pal <- function() {
  values <- c(
    "#999999", "#E69F00", "#56B4E9", "#009E73",
    "#0072B2", "#D55E00", "#CC79A7", "#F0E442"
  )
  f <- manual_pal(values)
  attr(f, "max_n") <- length(values)
  f
}

#' Discrete color-blind scale for SSD Plots
#'
#' The functions were designed for coloring different groups in a plot of SSD data.
#'
#' @param ... Arguments passed to [ggplot2::discrete_scale()].
#' @family ggplot
#' @export
#' @examples
#' # Use the color-blind palette for a SSD plot
#' ssd_plot(ssddata::ccme_boron, boron_pred, shape = "Group", color = "Group") +
#'   scale_colour_ssd()
scale_colour_ssd <- function(...) {
  discrete_scale("colour", palette = ssd_pal(), ...)
}

#' @describeIn scale_colour_ssd Discrete color-blind scale for SSD Plots
#' @export
scale_color_ssd <- function(...) {
  discrete_scale("colour", palette = ssd_pal(), ...)
}

#' @describeIn scale_colour_ssd Discrete color-blind scale for SSD Plots
#' @export
#' @examples
#' # Use the color-blind palette for a histogram of concentrations
#' ggplot2::ggplot(ssddata::ccme_boron, ggplot2::aes(x = Species, y = Conc, fill = Group)) +
#'   ggplot2::geom_col() +
#'   scale_fill_ssd() +
#'   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
scale_fill_ssd <- function(...) {
  discrete_scale("fill", palette = ssd_pal(), ...)
}

#' Species Sensitivity Data Points
#'
#' Uses the empirical cumulative distribution to create scatterplot of points `x`.
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#' @seealso [`ssd_plot_cdf()`]
#' @family ggplot
#' @export
#' @examples
#' ggplot2::ggplot(ssddata::ccme_boron, ggplot2::aes(x = Conc)) +
#'   geom_ssdpoint()
geom_ssdpoint <- function(mapping = NULL,
                          data = NULL,
                          stat = "ssdpoint",
                          position = "identity",
                          ...,
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSsdpoint,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

#' Species Sensitivity Censored Segments
#'
#' Uses the empirical cumulative distribution to draw lines between points `x` and `xend`.
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_segment
#' @seealso [`ssd_plot_cdf()`]
#' @family ggplot
#' @export
#' @examples
#' ggplot2::ggplot(ssddata::ccme_boron, ggplot2::aes(x = Conc, xend = Conc * 2)) +
#'   geom_ssdsegment()
geom_ssdsegment <- function(mapping = NULL,
                            data = NULL,
                            stat = "ssdsegment",
                            position = "identity",
                            ...,
                            arrow = NULL,
                            arrow.fill = NULL,
                            lineend = "butt",
                            linejoin = "round",
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSsdsegment,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      arrow = arrow, arrow.fill = arrow.fill,
      lineend = lineend, linejoin = linejoin, na.rm = na.rm, ...
    )
  )
}

#' Species Sensitivity Hazard Concentration Intersection
#'
#' Plots the intersection between each `xintercept` and `yintercept` value.
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_path
#' @inheritParams params
#' @seealso [`ssd_plot_cdf()`]
#' @family ggplot
#' @export
#' @examples
#' ggplot2::ggplot(ssddata::ccme_boron, ggplot2::aes(x = Conc)) +
#'   geom_ssdpoint() +
#'   geom_hcintersect(xintercept = 1.5, yintercept = 0.05)
geom_hcintersect <- function(mapping = NULL,
                             data = NULL,
                             ...,
                             xintercept,
                             yintercept,
                             na.rm = FALSE,
                             show.legend = NA) {
  if (!missing(xintercept)) {
    data <- data.frame(xintercept = xintercept)
    mapping <- aes(xintercept = xintercept)
    show.legend <- FALSE
  }

  if (!missing(yintercept)) {
    if (!missing(xintercept)) {
      data$yintercept <- yintercept
      mapping$yintercept <- yintercept
    } else {
      data <- data.frame(yintercept = yintercept)
      mapping <- aes(yintercept = yintercept)
    }
    show.legend <- FALSE
  }

  layer(
    data = data, mapping = mapping, stat = StatIdentity, geom = GeomHcintersect,
    position = PositionIdentity, show.legend = show.legend, inherit.aes = FALSE,
    params = list(na.rm = na.rm, ...)
  )
}

#' Ribbon on X-Axis
#'
#' Plots the `x` interval defined by `xmin` and `xmax`.
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#' @seealso [`ssd_plot_cdf()`]
#' @family ggplot
#' @export
#' @examples
#' gp <- ggplot2::ggplot(boron_pred) +
#'   geom_xribbon(ggplot2::aes(xmin = lcl, xmax = ucl, y = proportion))
geom_xribbon <- function(mapping = NULL,
                         data = NULL,
                         stat = "identity",
                         position = "identity",
                         ...,
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomXribbon,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

#' Species Sensitivity Data Points
#' `r lifecycle::badge('deprecated')`
#'
#' Deprecated for `geom_ssdpoint()`.
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#' @keywords internal
#' @export
geom_ssd <- function(mapping = NULL,
                     data = NULL,
                     stat = "ssdpoint",
                     position = "identity",
                     ...,
                     na.rm = FALSE,
                     show.legend = NA,
                     inherit.aes = TRUE) {
  lifecycle::deprecate_stop("0.3.5", "geom_ssd()", "geom_ssdpoint()")
}

#' Plot Species Sensitivity Data
#' `r lifecycle::badge('deprecated')`
#'
#' Uses the empirical cumulative density/distribution to visualize species sensitivity data.
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#' @seealso [`geom_ssdpoint()`]
#' @family ggplot2
#' @keywords internal
#' @export
stat_ssd <- function(mapping = NULL,
                     data = NULL,
                     geom = "point",
                     position = "identity",
                     ...,
                     na.rm = FALSE,
                     show.legend = NA,
                     inherit.aes = TRUE) {
  lifecycle::deprecate_stop("0.3.5", "stat_ssd()")
}
