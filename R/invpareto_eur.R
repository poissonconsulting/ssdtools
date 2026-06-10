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

#' @describeIn ssd_p Cumulative Distribution Function for European Inverse Pareto Distribution
#' @export
#' @examples
#'
#' ssd_pinvpareto_eur(1)
ssd_pinvpareto_eur <- function(
  q,
  shape = 3,
  scale = 1,
  lower.tail = TRUE,
  log.p = FALSE
) {
  pdist(
    "invpareto_eur",
    q = q,
    shape = shape,
    scale = scale,
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @describeIn ssd_q Quantile Function for European Inverse Pareto Distribution
#' @export
#' @examples
#'
#' ssd_qinvpareto_eur(0.5)
ssd_qinvpareto_eur <- function(
  p,
  shape = 3,
  scale = 1,
  lower.tail = TRUE,
  log.p = FALSE
) {
  qdist(
    "invpareto_eur",
    p = p,
    shape = shape,
    scale = scale,
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @describeIn ssd_r Random Generation for European Inverse Pareto Distribution
#' @export
#' @examples
#'
#' withr::with_seed(50, {
#'   x <- ssd_rinvpareto_eur(10000)
#' })
#' hist(x, breaks = 1000)
ssd_rinvpareto_eur <- function(n, shape = 3, scale = 1, chk = TRUE) {
  rdist("invpareto_eur", n = n, shape = shape, scale = scale, chk = chk)
}

#' @describeIn ssd_e Default Parameter Values for European Inverse Pareto Distribution
#' @export
#' @examples
#'
#' ssd_einvpareto_eur()
ssd_einvpareto_eur <- function() {
  list(shape = 3, scale = 1)
}

sinvpareto_eur <- function(data, pars = NULL) {
  if (!is.null(pars)) {
    return(pars)
  }
  # The European inverse Pareto has no closed-form mle, so use the median as
  # an initial scale; with shape = 1 this matches the empirical median exactly
  # (F(scale) = 0.5) and provides a stable starting point for optimisation.
  scale <- stats::median(data$right, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    scale <- 1
  }
  list(log_shape = log(1), log_scale = log(scale))
}

pinvpareto_eur_ssd <- function(q, shape, scale) {
  if (shape <= 0 || scale <= 0) {
    return(NaN)
  }
  pow(q / (q + scale), shape)
}

qinvpareto_eur_ssd <- function(p, shape, scale) {
  if (shape <= 0 || scale <= 0) {
    return(NaN)
  }
  u <- pow(p, 1 / shape)
  scale * u / (1 - u)
}
