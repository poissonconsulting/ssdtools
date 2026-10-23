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

#' Hazard Concentrations for Species Sensitivity Distributions
#'
#' Calculates concentration(s) with bootstrap confidence intervals
#' that protect specified proportion(s) of species for
#' individual or model-averaged distributions
#' using parametric or non-parametric bootstrapping.
#'
#' Model-averaged estimates and/or confidence intervals (including standard error)
#' can be calculated  by treating the distributions as
#' constituting a single mixture distribution
#' versus 'taking the mean'.
#' When calculating the model averaged estimates treating the
#' distributions as constituting a single mixture distribution
#' ensures that `ssd_hc()` is the inverse of `ssd_hp()`.
#'
#' Distributions with an absolute AIC difference greater
#' than a delta of by default 7 have considerably less support (wt < 0.01)
#' and are excluded
#' prior to calculation of the hazard concentrations to reduce the run time.
#'
#' @references
#'
#' Burnham, K.P., and Anderson, D.R. 2002. Model Selection and Multimodel Inference. Springer New York, New York, NY. doi:10.1007/b97636.
#'
#' @inheritParams params
#' @return A tibble of corresponding hazard concentrations.
#' @seealso [`predict.fitdists()`] and [`ssd_hp()`].
#' @export
ssd_hc <- function(x, ...) {
  UseMethod("ssd_hc")
}

.ssd_hc_dist <- function(x, dist, proportion) {
  fun <- ns_fun(paste0("ssd_q", dist))
  args <- list(p = proportion)
  args <- c(as.list(x), args)
  est <- do.call(fun, args)
  tibble(
    dist = dist,
    proportion = proportion,
    est = est,
    se = NA_real_,
    lcl = NA_real_,
    ucl = NA_real_,
    wt = 1,
    nboot = 0L,
    pboot = NA_real_
  )
}

# Resolve the deprecated `percent` argument into a `proportion` and
# validate it. Shared by the ssd_hc() and predict() methods.
# Shared argument handling for ssd_hc() and ssd_hp(). `fun` names the calling
# function in the deprecation messages. lifecycle::deprecate_soft() only warns
# when the deprecated function was called from the global environment or from a
# test, and it works that out from its calling frames, so each helper passes
# the frames of the exported function (`env`) and of its caller (`user_env`)
# rather than its own.
.hcp_est_method <- function(est_method, multi_est, fun) {
  env <- parent.frame()
  user_env <- parent.frame(2)
  if (lifecycle::is_present(multi_est)) {
    lifecycle::deprecate_soft(
      "2.3.1",
      paste0(fun, "(multi_est)"),
      paste0(fun, "(est_method)"),
      env = env,
      user_env = user_env
    )
    chk_flag(multi_est)
    est_method <- if (multi_est) "multi" else "arithmetic"
  }
  chk_string(est_method)
  chk_subset(est_method, ssd_est_methods())
  est_method
}

.hcp_ci_method <- function(ci_method, fun) {
  env <- parent.frame()
  user_env <- parent.frame(2)
  chk_string(ci_method)
  if (ci_method == "weighted_arithmetic") {
    lifecycle::deprecate_soft(
      "2.3.1",
      I(paste0(fun, "(ci_method = 'weighted_arithmetic')")),
      I(paste0(fun, "(ci_method = 'MACL')")),
      env = env,
      user_env = user_env
    )
    ci_method <- "MACL"
  }
  chk_subset(ci_method, ssd_ci_methods())
  ci_method
}

# `deprecated` is whether to warn that percentages are still the default and
# `fun` names the function whose argument the user should set.
.hp_proportion <- function(proportion, deprecated, fun, id = NULL) {
  env <- parent.frame()
  user_env <- parent.frame(2)
  if (deprecated) {
    lifecycle::deprecate_soft(
      "2.3.1",
      I("ssd_hp(proportion = FALSE)"),
      I("ssd_hp(proportion = TRUE)"),
      paste0(
        "Please set the `proportion` argument to `",
        fun,
        "()` to be TRUE which will cause it to return hazard proportions instead of percentages then update your downstream code accordingly."
      ),
      id = id,
      env = env,
      user_env = user_env
    )
  }
  chk_flag(proportion)
  proportion
}

.hc_proportion <- function(percent, proportion) {
  if (lifecycle::is_present(percent)) {
    lifecycle::deprecate_soft(
      "2.0.0",
      "ssd_hc(percent)",
      "ssd_hc(proportion)",
      id = "hc"
    )
    chk_vector(percent)
    chk_numeric(percent)
    chk_range(percent, c(0, 100))
    proportion <- percent / 100
  }
  chk_vector(proportion)
  chk_numeric(proportion)
  chk_range(proportion)
  proportion
}

#' @describeIn ssd_hc Hazard Concentrations for Distributional Estimates
#' @export
#' @examples
#'
#' ssd_hc(ssd_match_moments())
ssd_hc.list <- function(
  x,
  percent,
  proportion = 0.05,
  ...
) {
  chk_list(x)
  chk_named(x)
  chk_unique(names(x))
  chk_unused(...)

  proportion <- .hc_proportion(percent, proportion)

  if (!length(x)) {
    hc <- no_hcp()
    hc <- dplyr::rename(hc, proportion = "value")
    return(hc)
  }
  hc <- mapply(
    .ssd_hc_dist,
    x,
    names(x),
    MoreArgs = list(proportion = proportion),
    SIMPLIFY = FALSE
  )
  dplyr::bind_rows(hc)
}

#' @describeIn ssd_hc Hazard Concentrations for fitdists Object
#' @export
#' @examples
#'
#' fits <- ssd_fit_dists(ssddata::ccme_boron)
#' ssd_hc(fits)
ssd_hc.fitdists <- function(
  x,
  percent = deprecated(),
  proportion = 0.05,
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
  samples = FALSE,
  save_to = NULL,
  control = NULL
) {
  chk_unused(...)

  proportion <- .hc_proportion(percent, proportion)
  est_method <- .hcp_est_method(est_method, multi_est, "ssd_hc")
  ci_method <- .hcp_ci_method(ci_method, "ssd_hc")

  hcp <- hcp(
    x = x,
    value = proportion,
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
    samples = samples,
    save_to = save_to,
    hc = TRUE
  )

  hcp <- dplyr::rename(hcp, proportion = "value")
  hcp
}

#' @describeIn ssd_hc Hazard Concentrations for fitburrlioz Object
#' @export
#' @examples
#'
#' fit <- ssd_fit_burrlioz(ssddata::ccme_boron)
#' ssd_hc(fit)
ssd_hc.fitburrlioz <- function(
  x,
  percent,
  proportion = 0.05,
  ...,
  ci = FALSE,
  level = 0.95,
  nboot = 1000,
  min_pboot = 0.8,
  parametric = FALSE,
  samples = FALSE,
  save_to = NULL
) {
  chk_length(x, upper = 1L)
  chk_named(x)
  chk_subset(names(x), c("burrIII3", "invpareto", "llogis", "lgumbel"))
  chk_unused(...)

  proportion <- .hc_proportion(percent, proportion)

  fun <- if (names(x) == "burrIII3") fit_burrlioz else fit_tmb

  hcp <- hcp(
    x = x,
    value = proportion,
    ci = ci,
    level = level,
    nboot = nboot,
    average = FALSE,
    est_method = "multi",
    delta = Inf,
    min_pboot = min_pboot,
    parametric = parametric,
    ci_method = "multi_free",
    save_to = save_to,
    samples = samples,
    control = NULL,
    hc = TRUE,
    fun = fun
  )

  hcp <- dplyr::rename(hcp, proportion = "value")
  hcp
}
