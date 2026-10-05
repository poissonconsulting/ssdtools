# Copyright 2015-2023 Province of British Columbia
# Copyright 2021 Environment and Climate Change Canada
# Copyright 2023-2025 Australian Government Department of Climate Change,
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

#' Hazard Proportion
#'
#' Calculates proportion of species affected at specified concentration(s)
#' with quantile based bootstrap confidence intervals for
#' individual or model-averaged distributions
#' using parametric or non-parametric bootstrapping.
#' For more information see the inverse function [`ssd_hc()`].
#'
#' @inheritParams params
#' @param proportion A flag specifying whether to return hazard proportions
#' (`proportion = TRUE`) or hazard percentages (`proportion = FALSE`).
#' To not break existing code the default value is `FALSE` but
#' will be switching the default to `TRUE` in a future version.
#' The user is recommended to manually set to `TRUE` now to avoid
#' unexpected changes in future versions.
#' @return A tibble of corresponding hazard proportions.
#' @seealso [`ssd_hc()`]
#' @export
#' @examples
#' fits <- ssd_fit_dists(ssddata::ccme_boron)
#' ssd_hp(fits, conc = 1)
ssd_hp <- function(x, ...) {
  UseMethod("ssd_hp")
}

#' @describeIn ssd_hp Hazard Proportions for fitdists Object
#' @export
ssd_hp.fitdists <- function(
  x,
  conc = 1,
  ...,
  average = TRUE,
  ci = FALSE,
  level = 0.95,
  nboot = 1000,
  min_pboot = 0.8,
  multi_est = deprecated(),
  est_method = "multi",
  ci_method = "weighted_samples",
  parametric = TRUE,
  delta = 9.21,
  proportion = FALSE,
  samples = FALSE,
  save_to = NULL,
  control = NULL
) {
  chk_vector(conc)
  chk_numeric(conc)
  chk_unused(...)

  est_method <- .hcp_est_method(est_method, multi_est, "ssd_hp")
  ci_method <- .hcp_ci_method(ci_method, "ssd_hp")
  proportion <- .hp_proportion(
    proportion,
    missing(proportion),
    "ssd_hp",
    id = "ssd_hp"
  )

  hcp <- hcp(
    x = x,
    value = conc,
    ci = ci,
    level = level,
    nboot = nboot,
    average = average,
    est_method = est_method,
    delta = delta,
    min_pboot = min_pboot,
    parametric = parametric,
    ci_method = ci_method,
    control = control,
    save_to = save_to,
    samples = samples,
    hc = FALSE
  )
  hcp <- dplyr::rename(hcp, conc = "value")
  if (!proportion) {
    hcp <- hcp_unscale(hcp, 100)
  }
  hcp
}


#' @describeIn ssd_hp Hazard Proportions for fitburrlioz Object
#' @export
#' @examples
#'
#' fit <- ssd_fit_burrlioz(ssddata::ccme_boron)
#' ssd_hp(fit)
ssd_hp.fitburrlioz <- function(
  x,
  conc = 1,
  ...,
  ci = FALSE,
  level = 0.95,
  nboot = 1000,
  min_pboot = 0.8,
  parametric = FALSE,
  proportion = FALSE,
  samples = FALSE,
  save_to = NULL
) {
  chk_length(x, upper = 1L)
  chk_named(x)
  chk_subset(names(x), c("burrIII3", "invpareto", "llogis", "lgumbel"))
  chk_vector(conc)
  chk_numeric(conc)
  chk_flag(ci)
  chk_unused(...)

  proportion <- .hp_proportion(
    proportion,
    missing(proportion),
    "ssd_hp",
    id = "ssd_hp"
  )

  fun <- if (names(x) == "burrIII3") fit_burrlioz else fit_tmb

  hcp <- hcp(
    x = x,
    value = conc,
    ci = ci,
    level = level,
    nboot = nboot,
    average = FALSE,
    est_method = "multi",
    delta = Inf,
    min_pboot = min_pboot,
    parametric = parametric,
    ci_method = "multi_free",
    control = NULL,
    save_to = save_to,
    samples = samples,
    hc = FALSE,
    fun = fun
  )

  hcp <- dplyr::rename(hcp, conc = "value")
  if (!proportion) {
    hcp <- hcp_unscale(hcp, 100)
  }
  hcp
}
