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

#' @describeIn ssd_p Cumulative Distribution Function for Weibull Distribution
#' @export
#' @examples
#'
#' ssd_pweibull(1)
ssd_pweibull <- function(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  pdist("weibull",
    q = q, shape = shape, scale = scale,
    lower.tail = lower.tail, log.p = log.p
  )
}

#' @describeIn ssd_q Cumulative Distribution Function for Weibull Distribution
#' @export
#' @examples
#'
#' ssd_qweibull(0.5)
ssd_qweibull <- function(p, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  qdist("weibull",
    p = p, shape = shape, scale = scale,
    lower.tail = lower.tail, log.p = log.p
  )
}

#' @describeIn ssd_r Random Generation for Weibull Distribution
#' @export
#' @examples
#'
#' withr::with_seed(50, {
#'   x <- ssd_rweibull(10000)
#' })
#' hist(x, breaks = 1000)
ssd_rweibull <- function(n, shape = 1, scale = 1, chk = TRUE) {
  rdist("weibull", n = n, shape = shape, scale = scale, chk = chk)
}

#' @describeIn ssd_e Default Parameter Values for Log-Normal Distribution
#' @export
#' @examples
#'
#' ssd_eweibull()
ssd_eweibull <- function() {
  list(shape = 1, scale = 1)
}

sweibull <- function(data, pars = NULL) {
  if (!is.null(pars)) {
    return(pars)
  }

  x <- mean_weighted_values(data)
  n <- length(x)
  p <- (1:n - 0.3) / (n + 0.4)
  m <- stats::lm(log(-log(1 - p)) ~ log(sort(x)))

  shape <- m$coefficients[2]
  log_scale <- -m$coefficients[1] / shape
  list(
    log_scale = log_scale,
    log_shape = log(shape)
  )
}

pweibull_ssd <- function(q, shape, scale) {
  if (shape <= 0 || scale <= 0) {
    return(NaN)
  }
  stats::pweibull(q, shape, scale)
}

qweibull_ssd <- function(p, shape, scale) {
  if (shape <= 0 || scale <= 0) {
    return(NaN)
  }
  stats::qweibull(p, shape, scale)
}

rweibull_ssd <- function(n, shape, scale) {
  if (shape <= 0 || scale <= 0) {
    return(rep(NaN, n))
  }
  stats::rweibull(n, shape, scale)
}
