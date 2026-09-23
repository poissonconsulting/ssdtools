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

#' @describeIn ssd_p Cumulative Distribution Function for Log-Asymmetric-Triangular Distribution
#' @export
#' @examples
#'
#' ssd_pltriangle3(1)
ssd_pltriangle3 <- function(
  q,
  locationlog = 0,
  scalelog = 3,
  skew = 0,
  lower.tail = TRUE,
  log.p = FALSE
) {
  pdist(
    "triangle3",
    q = q,
    location = locationlog,
    scale = scalelog,
    skew = skew,
    lower.tail = lower.tail,
    log.p = log.p,
    .lgt = TRUE
  )
}

#' @describeIn ssd_q Quantile Function for Log-Asymmetric-Triangular Distribution
#' @export
#' @examples
#'
#' ssd_qltriangle3(0.5)
ssd_qltriangle3 <- function(
  p,
  locationlog = 0,
  scalelog = 3,
  skew = 0,
  lower.tail = TRUE,
  log.p = FALSE
) {
  qdist(
    "triangle3",
    p = p,
    location = locationlog,
    scale = scalelog,
    skew = skew,
    lower.tail = lower.tail,
    log.p = log.p,
    .lgt = TRUE
  )
}

#' @describeIn ssd_r Random Generation for Log-Asymmetric-Triangular Distribution
#' @export
#' @examples
#'
#' withr::with_seed(50, {
#'   x <- ssd_rltriangle3(10000)
#' })
#' hist(x, breaks = 1000)
ssd_rltriangle3 <- function(
  n,
  locationlog = 0,
  scalelog = 3,
  skew = 0,
  chk = TRUE
) {
  rdist(
    "triangle3",
    n = n,
    location = locationlog,
    scale = scalelog,
    skew = skew,
    .lgt = TRUE,
    chk = chk
  )
}

#' @describeIn ssd_e Default Parameter Values for Log-Asymmetric-Triangular Distribution
#' @export
#' @examples
#'
#' ssd_eltriangle3()
ssd_eltriangle3 <- function() {
  list(locationlog = 0, scalelog = 3, skew = 0)
}

sltriangle3 <- function(data, pars = NULL) {
  if (!is.null(pars)) {
    return(pars)
  }

  x <- mean_weighted_values(data)
  logx <- log(x)
  logx <- logx[is.finite(logx)]
  # fall back to a small positive value when the data have no spread on the
  # log scale so the starting value stays finite
  halfwidth <- diff(range(logx)) / 2 * 1.1
  if (!length(logx) || !is.finite(halfwidth) || halfwidth <= 0) {
    return(c(sltriangle(data), list(skew = 0)))
  }
  # place the support just beyond the data and set the mode from the mean of
  # a triangular distribution, (lower + mode + upper) / 3, keeping it strictly
  # inside the support. The log-likelihood has a kink at every observation, so
  # a start on the correct side of the mode matters more than for the smooth
  # distributions.
  lower <- mean(range(logx)) - halfwidth
  upper <- mean(range(logx)) + halfwidth
  location <- 3 * mean(logx) - lower - upper
  location <- min(
    max(location, lower + 0.05 * halfwidth),
    upper - 0.05 * halfwidth
  )

  list(
    locationlog = location,
    log_scalelog = log(halfwidth),
    skew = 1 - (location - lower) / halfwidth
  )
}

bltriangle3 <- function(...) {
  list(
    lower = list(locationlog = -Inf, log_scalelog = -Inf, skew = -1),
    upper = list(locationlog = Inf, log_scalelog = Inf, skew = 1)
  )
}

# Asymmetric triangular distribution with mode `location`, half-width `scale`
# and skewness `skew` in [-1, 1]. The support is
#   [location - (1 - skew) * scale, location + (1 + skew) * scale]
# so it has total width 2 * scale, `skew = 0` is the symmetric triangular and
# `skew = 1` (`skew = -1`) puts the mode at the lower (upper) limit.
triangle3_limits <- function(location, scale, skew) {
  list(
    lower = location - (1 - skew) * scale,
    upper = location + (1 + skew) * scale
  )
}

ptriangle3_ssd <- function(q, location, scale, skew) {
  if (scale <= 0 || skew < -1 || skew > 1) {
    return(rep(NaN, length(q)))
  }
  limits <- triangle3_limits(location, scale, skew)
  a <- limits$lower
  b <- limits$upper
  w <- 2 * scale
  # the innermost branches are only evaluated for a < q <= location and
  # location < q < b, so their denominators are positive
  ifelse(
    q <= a,
    0,
    ifelse(
      q >= b,
      1,
      ifelse(
        q <= location,
        (q - a)^2 / (w * (location - a)),
        1 - (b - q)^2 / (w * (b - location))
      )
    )
  )
}

qtriangle3_ssd <- function(p, location, scale, skew) {
  if (scale <= 0 || skew < -1 || skew > 1) {
    return(rep(NaN, length(p)))
  }
  limits <- triangle3_limits(location, scale, skew)
  a <- limits$lower
  b <- limits$upper
  w <- 2 * scale
  # cumulative probability at the mode
  pmode <- (1 - skew) / 2
  q <- ifelse(
    p <= pmode,
    a + sqrt(p * w * (location - a)),
    b - sqrt((1 - p) * w * (b - location))
  )
  # return the limits exactly rather than to floating point tolerance
  q[p == 0] <- a
  q[p == 1] <- b
  q
}

rtriangle3_ssd <- function(n, location, scale, skew) {
  if (scale <= 0 || skew < -1 || skew > 1) {
    return(rep(NaN, n))
  }
  qtriangle3_ssd(stats::runif(n), location, scale, skew)
}

pltriangle3_ssd <- function(q, locationlog, scalelog, skew) {
  ptriangle3_ssd(log(q), location = locationlog, scale = scalelog, skew = skew)
}

qltriangle3_ssd <- function(p, locationlog, scalelog, skew) {
  exp(qtriangle3_ssd(p, location = locationlog, scale = scalelog, skew = skew))
}
