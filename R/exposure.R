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

#' Proportion Exposure
#'
#' Calculates average proportion exposed based on log-normal distribution of concentrations.
#'
#' @inheritParams params
#' @param meanlog The mean of the exposure concentrations on the log scale.
#' @param sdlog The standard deviation of the exposure concentrations on the log scale.
#' @param nboot The number of samples to use to calculate the exposure.
#' @return The proportion exposed.
#' @export
#' @examples
#' \dontrun{
#' fits <- ssd_fit_dists(ssddata::ccme_boron, dists = "lnorm")
#' withr::with_seed(50, {
#'   ssd_exposure(fits)
#'   ssd_exposure(fits, meanlog = 1)
#'   ssd_exposure(fits, meanlog = 1, sdlog = 1)
#' })
#' }
ssd_exposure <- function(x, meanlog = 0, sdlog = 1, ..., nboot = 1000) {
  chk_number(meanlog)
  chk_number(sdlog)
  chk_unused(...)
  chk_whole_number(nboot)
  chk_gt(nboot)
  conc <- ssd_rlnorm(nboot, meanlog = meanlog, sdlog = sdlog)
  mean(ssd_hp(x, conc, proportion = TRUE)$est)
}
