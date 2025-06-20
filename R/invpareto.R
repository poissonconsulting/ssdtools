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

#' @describeIn ssd_p Cumulative Distribution Function for Inverse Pareto Distribution
#' @export
#' @examples
#'
#' ssd_pinvpareto(1)
ssd_pinvpareto <- function(q, shape = 3, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  pdist("invpareto",
    q = q, shape = shape, scale = scale,
    lower.tail = lower.tail, log.p = log.p
  )
}

#' @describeIn ssd_q Quantile Function for Inverse Pareto Distribution
#' @export
#' @examples
#'
#' ssd_qinvpareto(0.5)
ssd_qinvpareto <- function(p, shape = 3, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  qdist("invpareto",
    p = p, shape = shape, scale = scale,
    lower.tail = lower.tail, log.p = log.p
  )
}

#' @describeIn ssd_r Random Generation for Inverse Pareto Distribution
#' @export
#' @examples
#'
#' withr::with_seed(50, {
#' x <- ssd_rinvpareto(10000)
#' })
#' hist(x, breaks = 1000)
ssd_rinvpareto <- function(n, shape = 3, scale = 1, chk = TRUE) {
  rdist("invpareto", n = n, shape = shape, scale = scale, chk = chk)
}

#' @describeIn ssd_e Default Parameter Values for Inverse Pareto Distribution
#' @export
#' @examples
#'
#' ssd_einvpareto()
ssd_einvpareto <- function() {
  list(shape = 3, scale = 1)
}

sinvpareto <- function(data, pars = NULL) {
  scale <- max(data$right)
  shape <- 1 / mean(log(scale / data$right))

  n <- length(data$right)
  scale <- scale * (shape * n) / (shape * n - 1)
  shape <- 1 / mean(log(scale / data$right))

  spars <- list(log_scale = log(scale), log_shape = log(shape))
  if (!is.null(pars)) { # use new bias corrected order statistic
    pars$log_scale <- spars$log_scale
    return(pars)
  }
  spars
}

minvpareto <- function() {
  list(log_scale = factor(NA))
}

pinvpareto_ssd <- function(q, shape, scale) {
  if (shape <= 0 || scale <= 0) {
    return(NaN)
  }
  pow((q / scale), shape)
}

qinvpareto_ssd <- function(p, shape, scale) {
  if (shape <= 0 || scale <= 0) {
    return(NaN)
  }
  pow(p, 1 / shape) * scale
}
