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

#' @describeIn ssd_p Cumulative Distribution Function for Log-Normal Distribution
#' @export
#' @examples
#'
#' ssd_plnorm(1)
ssd_plnorm <- function(q, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE) {
  pdist("lnorm",
    q = q, meanlog = meanlog, sdlog = sdlog,
    lower.tail = lower.tail, log.p = log.p
  )
}

#' @describeIn ssd_q Cumulative Distribution Function for Log-Normal Distribution
#' @export
#' @examples
#'
#' ssd_qlnorm(0.5)
ssd_qlnorm <- function(p, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE) {
  qdist("lnorm",
    p = p, meanlog = meanlog, sdlog = sdlog,
    lower.tail = lower.tail, log.p = log.p
  )
}

#' @describeIn ssd_r Random Generation for Log-Normal Distribution
#' @export
#' @examples
#'
#' withr::with_seed(50, {
#'   x <- ssd_rlnorm(10000)
#' })
#' hist(x, breaks = 1000)
ssd_rlnorm <- function(n, meanlog = 0, sdlog = 1, chk = TRUE) {
  rdist("lnorm", n = n, meanlog = meanlog, sdlog = sdlog, chk = chk)
}

#' @describeIn ssd_e Default Parameter Values for Log-Normal Distribution
#' @export
#' @examples
#'
#' ssd_elnorm()
ssd_elnorm <- function() {
  list(meanlog = 0, sdlog = 1)
}

slnorm <- function(data, pars = NULL) {
  if (!is.null(pars)) {
    return(pars)
  }

  x <- mean_weighted_values(data)

  list(
    meanlog = mean(log(x), na.rm = TRUE) * (1 + 1e-2),
    log_sdlog = log(sd(log(x), na.rm = TRUE)) * (1 - 1e-2)
  )
}

plnorm_ssd <- function(q, meanlog, sdlog) {
  if (sdlog <= 0) {
    return(NaN)
  }
  stats::plnorm(q, meanlog, sdlog)
}

qlnorm_ssd <- function(p, meanlog, sdlog) {
  if (sdlog <= 0) {
    return(NaN)
  }
  stats::qlnorm(p, meanlog, sdlog)
}

rlnorm_ssd <- function(n, meanlog, sdlog) {
  if (sdlog <= 0) {
    return(rep(NaN, n))
  }
  stats::rlnorm(n, meanlog, sdlog)
}
