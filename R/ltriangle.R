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

#' @describeIn ssd_p Cumulative Distribution Function for Log-Triangular Distribution
#' @export
#' @examples
#'
#' ssd_pltriangle(1)
ssd_pltriangle <- function(
  q,
  locationlog = 0,
  scalelog = 3,
  lower.tail = TRUE,
  log.p = FALSE
) {
  pdist(
    "triangle",
    q = q,
    location = locationlog,
    scale = scalelog,
    lower.tail = lower.tail,
    log.p = log.p,
    .lgt = TRUE
  )
}

#' @describeIn ssd_q Quantile Function for Log-Triangular Distribution
#' @export
#' @examples
#'
#' ssd_qltriangle(0.5)
ssd_qltriangle <- function(
  p,
  locationlog = 0,
  scalelog = 3,
  lower.tail = TRUE,
  log.p = FALSE
) {
  qdist(
    "triangle",
    p = p,
    location = locationlog,
    scale = scalelog,
    lower.tail = lower.tail,
    log.p = log.p,
    .lgt = TRUE
  )
}

#' @describeIn ssd_r Random Generation for Log-Triangular Distribution
#' @export
#' @examples
#'
#' withr::with_seed(50, {
#'   x <- ssd_rltriangle(10000)
#' })
#' hist(x, breaks = 1000)
ssd_rltriangle <- function(n, locationlog = 0, scalelog = 3, chk = TRUE) {
  rdist(
    "triangle",
    n = n,
    location = locationlog,
    scale = scalelog,
    .lgt = TRUE,
    chk = chk
  )
}

#' @describeIn ssd_e Default Parameter Values for Log-Triangular Distribution
#' @export
#' @examples
#'
#' ssd_eltriangle()
ssd_eltriangle <- function() {
  list(locationlog = 0, scalelog = 3)
}

sltriangle <- function(data, pars = NULL) {
  if (!is.null(pars)) {
    return(pars)
  }

  x <- mean_weighted_values(data)
  logx <- log(x)
  location <- mean(logx, na.rm = TRUE)
  # half-width must exceed the largest deviation so the support covers the data
  halfwidth <- max(abs(logx - location), na.rm = TRUE)

  list(
    locationlog = location,
    log_scalelog = log(halfwidth * 1.1)
  )
}

# Symmetric triangular distribution with mode `location` and half-width
# `scale` (support is `location` +/- `scale`).
ptriangle_ssd <- function(q, location, scale) {
  if (scale <= 0) {
    return(rep(NaN, length(q)))
  }
  z <- (q - location) / scale
  ifelse(
    z <= -1,
    0,
    ifelse(
      z <= 0,
      (z + 1)^2 / 2,
      ifelse(z < 1, 1 - (1 - z)^2 / 2, 1)
    )
  )
}

qtriangle_ssd <- function(p, location, scale) {
  if (scale <= 0) {
    return(rep(NaN, length(p)))
  }
  z <- ifelse(p <= 0.5, sqrt(2 * p) - 1, 1 - sqrt(2 * (1 - p)))
  location + scale * z
}

rtriangle_ssd <- function(n, location, scale) {
  if (scale <= 0) {
    return(rep(NaN, n))
  }
  location + scale * (stats::runif(n) + stats::runif(n) - 1)
}

pltriangle_ssd <- function(q, locationlog, scalelog) {
  ptriangle_ssd(log(q), location = locationlog, scale = scalelog)
}

qltriangle_ssd <- function(p, locationlog, scalelog) {
  exp(qtriangle_ssd(p, location = locationlog, scale = scalelog))
}
