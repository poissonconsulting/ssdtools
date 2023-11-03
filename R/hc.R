#    Copyright 2021 Province of British Columbia
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
#' Gets concentration(s) that protect specified percentage(s) of species.
#'
#' If `ci = TRUE` uses parameteric bootstrapping to get confidence intervals on the
#' hazard concentrations(s).
#'
#' @inheritParams params
#' @return A tibble of corresponding hazard concentrations.
#' @seealso [`predict.fitdists()`] and [`ssd_hp()`].
#' @export
#' @examples
#' fits <- ssd_fit_dists(ssddata::ccme_boron)
#' ssd_hc(fits)
#' ssd_hc(estimates(fits))
#' ssd_hc(ssd_match_moments())
ssd_hc <- function(x, ...) {
  UseMethod("ssd_hc")
}

no_ssd_hc <- function() {
  tibble(
    dist = character(0),
    percent = numeric(0),
    est = numeric(0),
    se = numeric(0),
    lcl = numeric(0),
    ucl = numeric(0),
    wt = numeric(0),
    nboot = integer(0),
    pboot = numeric(0)
  )
}

.ssd_hc_dist <- function(x, dist, proportion) {
  fun <- paste0("ssd_q", dist)
  args <- list(p = proportion)
  args <- c(as.list(x), args)
  est <- do.call(fun, args)
  tibble(
    dist = dist,
    percent = proportion * 100, est = est,
    se = NA_real_, lcl = NA_real_, ucl = NA_real_,
    wt = 1,
    nboot = 0L, pboot = NA_real_
  )
}

.ssd_hc_tmbfit <- function(x, proportion, ci, level, nboot, min_pboot,
                           data, rescale, weighted, censoring, min_pmix,
                           range_shape1, range_shape2, parametric,
                           keep_samples, control) {
  args <- estimates(x)
  args$p <- proportion
  dist <- .dist_tmbfit(x)
  what <- paste0("ssd_q", dist)

  est <- do.call(what, args)
  if (!ci) {
    na <- rep(NA_real_, length(proportion))
    return(tibble(
      dist = rep(dist, length(proportion)),
      percent = proportion * 100,
      est = est * rescale,
      se = na,
      lcl = na,
      ucl = na,
      wt = rep(1, length(proportion)),
      nboot = rep(0L, length(proportion)),
      pboot = na,
      samples = na
    ))
  }
  censoring <- censoring / rescale
  fun <- safely(fit_tmb)

  estimates <- boot_estimates(x,
    fun = fun, nboot = nboot, data = data, weighted = weighted,
    censoring = censoring, min_pmix = min_pmix,
    range_shape1 = range_shape1,
    range_shape2 = range_shape2,
    parametric = parametric,
    control = control
  )

  samples <- sample_estimates(estimates, what, x = proportion)
  samples <- lapply(samples, FUN = function(x) x * rescale)

  cis <- cis_estimates(estimates, what, level = level, x = proportion)
  hc <- tibble(
    dist = dist,
    percent = proportion * 100, est = est * rescale,
    se = cis$se * rescale, lcl = cis$lcl * rescale, ucl = cis$ucl * rescale,
    wt = rep(1, length(proportion)),
    nboot = nboot, pboot = length(estimates) / nboot,
    samples = samples
  )
  replace_min_pboot_na(hc, min_pboot)
}

.ssd_hc_fitdists <- function(
    x,
    percent,
    ci,
    level,
    nboot,
    average,
    min_pboot,
    parametric,
    averaging_method,
    keep_samples,
    control) {
  if (!length(x) || !length(percent)) {
    return(no_ssd_hc())
  }

  if (is.null(control)) {
    control <- .control_fitdists(x)
  }

  data <- .data_fitdists(x)
  rescale <- .rescale_fitdists(x)
  censoring <- .censoring_fitdists(x)
  min_pmix <- .min_pmix_fitdists(x)
  range_shape1 <- .range_shape1_fitdists(x)
  range_shape2 <- .range_shape2_fitdists(x)
  weighted <- .weighted_fitdists(x)
  unequal <- .unequal_fitdists(x)
  wt_est_nest <- wt_est_nest(x)
  weight <- wt_est_nest$weight
  return_vars <- c(
    "dist", "percent", "est", "se", "lcl", "ucl",
    "wt", "method", "nboot", "pboot"
  )
  if (keep_samples) {
    return_vars <- c(return_vars, "samples")
  }

  if (parametric && ci && identical(censoring, c(NA_real_, NA_real_))) {
    wrn("Parametric CIs cannot be calculated for inconsistently censored data.")
    ci <- FALSE
  }

  if (keep_samples && !ci) {
    wrn("Bootstrap samples are only generated when ci = TRUE. Setting ci = TRUE.")
    ci <- TRUE
  }

  if (parametric && ci && unequal) {
    wrn("Parametric CIs cannot be calculated for unequally weighted data.")
    ci <- FALSE
  }

  if (averaging_method == "weighted_sample" && !ci) {
    wrn("The weighted sample averaging method can only be used when bootstrapping CIs is used. Setting ci = TRUE")
    ci <- TRUE
  }

  if (!ci) {
    nboot <- 0L
  }

  seeds <- seed_streams(length(x))

  if (!average) {
    hc <- future_map(x, .ssd_hc_tmbfit,
      proportion = percent / 100, ci = ci, level = level, nboot = nboot,
      min_pboot = min_pboot,
      data = data, rescale = rescale, weighted = weighted,
      censoring = censoring,
      min_pmix = min_pmix, range_shape1 = range_shape1,
      range_shape2 = range_shape2,
      parametric = parametric, keep_samples = keep_samples, control = control,
      .options = furrr::furrr_options(seed = seeds)
    )
    hc <- mapply(
      function(x, y) {
        x$wt <- y
        x
      },
      x = hc, y = weight,
      USE.NAMES = FALSE, SIMPLIFY = FALSE
    )
    hc <- bind_rows(hc)
    hc$method <- if (parametric) "parametric" else "non-parametric"

    hc <- hc[c(return_vars)]
    return(hc)
  }

  if (averaging_method == "weighted_sample") {
    nboot_vals <- round(nboot * weight)
    hc <- furrr::future_map2(
      .x = x, .y = nboot_vals,
      ~ .ssd_hc_tmbfit(
        x = .x,
        proportion = percent / 100, ci = ci,
        level = level,
        nboot = .y,
        min_pboot = min_pboot,
        data = data, rescale = rescale, weighted = weighted,
        censoring = censoring,
        min_pmix = min_pmix, range_shape1 = range_shape1,
        range_shape2 = range_shape2,
        parametric = parametric, keep_samples = keep_samples,
        control = control
      ),
      .options = furrr::furrr_options(seed = seeds)
    )

    hc <- lapply(hc, FUN = function(x) x |>
      tidyr::unnest_longer(samples)) |>
      bind_rows()
    pboot_chk <- hc |>
      dplyr::select(dist, pboot) |>
      unique() |>
      dplyr::filter(pboot < min_pboot)
    dists_fail <- paste(pboot_chk$dist, collapse = "; ")

    if (nrow(pboot_chk) > 0) {
      warning(paste(dists_fail, " fail the minimum bootstrap convergence criteria of ",
        min_pboot, ".",
        sep = ""
      ))
    }

    new_pboot <- nrow(hc) / length(percent) / nboot
    method <- if (parametric) "parametric" else "non-parametric"

    if (ci) {
      sample_vals <- hc |> dplyr::select(samples)
      hc <- hc |>
        dplyr::select(percent, samples) |>
        dplyr::group_by(percent) |>
        dplyr::summarise(
          est = median(samples),
          lcl = quantile(samples, probs = probs(level)[1], na.rm = TRUE),
          ucl = quantile(samples, probs = probs(level)[2], na.rm = TRUE),
          se = sd(samples)
        ) |>
        dplyr::mutate(
          percent = as.integer(percent),
          dist = "average",
          method = method,
          nboot = nboot,
          pboot = new_pboot,
          wt = 1,
          samples = list(sample_vals)
        ) |>
        dplyr::select(tidyselect::all_of(return_vars))
      return(hc)
    }
  }

  if (averaging_method == "uniroot") {
    if (!ci) {
      seeds <- seed_streams(length(percent))
      hc <- future_map(
        percent / 100, .ssd_hc_root,
        wt_est_nest = wt_est_nest, ci = ci, level = level, nboot = nboot,
        min_pboot = min_pboot,
        data = data, rescale = rescale, weighted = weighted, censoring = censoring,
        min_pmix = min_pmix, range_shape1 = range_shape1, range_shape2 = range_shape2,
        parametric = parametric,
        control = control,
        .options = furrr::furrr_options(seed = seeds)
      ) |>
        dplyr::bind_rows() |>
        dplyr::mutate(
          method = ifelse(parametric == TRUE, "parametric", "non-parametric"),
          dist = "average",
          percent = as.integer(percent),
          wt = rep(1, length(percent)),
          nboot = nboot
        ) |>
        dplyr::select(tidyselect::all_of(return_vars))
      return(hc)
    }

    if (ci) {
      seeds <- seed_streams(length(x))
      estimates <- future_map(
        x, .ssd_estimates,
        nboot = nboot,
        data = data, rescale = rescale, weighted = weighted,
        censoring = censoring,
        min_pmix = min_pmix, range_shape1 = range_shape1,
        range_shape2 = range_shape2,
        parametric = parametric, control = control,
        .options = furrr::furrr_options(seed = seeds)
      )

      wt_est_nest_boot <- lapply(1:nboot, FUN = function(n) {
        boot_estimates <- lapply(estimates, FUN = function(y) {
          y[[n]] |>
            dplyr::bind_rows(.id = "term") |>
            tidyr::pivot_longer(cols = everything(), names_to = "term", values_to = "est")
        }) |> tibble()
        colnames(boot_estimates) <- "data"
        boot_estimates |>
          dplyr::mutate(
            dist = names(estimates),
            weight = weight
          ) |>
          dplyr::select(dist, weight, data)
      })

      hcs_boot <- lapply(wt_est_nest_boot, FUN = function(x) {
        seeds <- seed_streams(length(percent))
        hcs <- lapply(percent, FUN = function(p) {
          .ssd_hc_root(
            proportion = p / 100,
            wt_est_nest = x, ci = ci, level = level, nboot = nboot,
            min_pboot = min_pboot,
            data = data, rescale = rescale, weighted = weighted, censoring = censoring,
            min_pmix = min_pmix, range_shape1 = range_shape1, range_shape2 = range_shape2,
            parametric = parametric, control = control
          )
        }) |>
          dplyr::bind_rows() |>
          dplyr::mutate(
            method = ifelse(parametric == TRUE, "parametric", "non-parametric"),
            dist = "average",
            percent = as.integer(percent),
            wt = rep(1, length(percent)),
            nboot = nboot
          ) |>
          dplyr::select(tidyselect::all_of(setdiff(return_vars, "samples")))
        return(hcs)
      })

      hc <- hcs_boot |>
        dplyr::bind_rows() |>
        dplyr::select(percent, est, pboot, method) |>
        dplyr::group_by(percent, pboot, method) |>
        dplyr::summarise(
          lcl = quantile(est, probs = probs(level)[1], na.rm = TRUE),
          ucl = quantile(est, probs = probs(level)[2], na.rm = TRUE),
          se = sd(est),
          est = median(est)
        ) |>
        dplyr::mutate(
          dist = "average",
          method = method,
          nboot = nboot,
          wt = 1
        )
      if (keep_samples) {
        sample_list <- hcs_boot |>
          dplyr::bind_rows() |>
          dplyr::select(percent, est) |>
          tidyr::nest(samples = est, .by = percent) |>
          dplyr::left_join(hc) |>
          dplyr::select(tidyr::any_of(return_vars))

        return(sample_list)
      }

      return(hc |> dplyr::select(tidyr::any_of(return_vars)))
    }
  }

  if (averaging_method == "arithmetic" || averaging_method == "geometric") {
    seeds <- seed_streams(length(x))
    hc_boot <- future_map(x, .ssd_hc_tmbfit,
      proportion = percent / 100, ci = ci, level = level, nboot = nboot,
      min_pboot = min_pboot,
      data = data, rescale = rescale, weighted = weighted, censoring = censoring,
      min_pmix = min_pmix, range_shape1 = range_shape1, range_shape2 = range_shape2,
      parametric = parametric, keep_samples = keep_samples, control = control,
      .options = furrr::furrr_options(seed = seeds)
    )

    hc <- lapply(hc_boot, function(x) x[c("percent", "est", "se", "lcl", "ucl", "pboot")])
    hc <- lapply(hc, as.matrix)
    hc <- Reduce(function(x, y) {
      abind(x, y, along = 3)
    }, hc)
    suppressMessages(min <- apply(hc, c(1, 2), min))

    if (averaging_method == "arithmetic") {
      suppressMessages(hc <- apply(hc, c(1, 2), weighted.mean, w = weight))
    } else {
      suppressMessages(hc <- apply(hc, c(1, 2), weighted.geomean, w = weight))
    }

    min <- as.data.frame(min)
    hc <- as.data.frame(hc) |>
      tibble() |>
      dplyr::mutate(
        method = ifelse(parametric, "parametric", "non-parametric"),
        dist = "average",
        percent = as.integer(percent),
        wt = rep(1, length(percent)),
        nboot = nboot, pboot = min$pboot
      )

    if (keep_samples) {
      sample_dat <- hc_boot |>
        dplyr::bind_rows() |>
        dplyr::select(dist, percent, samples) |>
        tidyr::nest(.by = percent) |>
        dplyr::rename(samples = data) |>
        dplyr::left_join(hc) |>
        dplyr::select(tidyr::any_of(return_vars))

      return(sample_dat)
    }

    return(hc |> dplyr::select(tidyr::any_of(return_vars)))
  }
}

#' @describeIn ssd_hc Hazard Concentrations for Distributional Estimates
#' @export
ssd_hc.list <- function(x, percent = 5, ...) {
  chk_list(x)
  chk_named(x)
  chk_unique(names(x))
  chk_unused(...)

  if (!length(x)) {
    return(no_ssd_hc())
  }
  hc <- mapply(.ssd_hc_dist, x, names(x),
    MoreArgs = list(proportion = percent / 100),
    SIMPLIFY = FALSE
  )
  bind_rows(hc)
}

#' @describeIn ssd_hc Hazard Concentrations for fitdists Object
#' @export
ssd_hc.fitdists <- function(
    x,
    percent = 5,
    ci = FALSE,
    level = 0.95,
    nboot = 1000,
    average = TRUE,
    delta = 7,
    min_pboot = 0.99,
    parametric = TRUE,
    averaging_method = "arithmetic",
    keep_samples = FALSE,
    control = NULL,
    ...) {
  chk_vector(percent)
  chk_numeric(percent)
  chk_range(percent, c(0, 100))
  chk_flag(ci)
  chk_number(level)
  chk_range(level)
  chk_whole_number(nboot)
  chk_gt(nboot)
  chk_flag(average)
  chk_number(delta)
  chk_gte(delta)
  chk_number(min_pboot)
  chk_range(min_pboot)
  chk_flag(parametric)
  # chk_flag(averaging_method)
  chk_null_or(control, vld = vld_list)
  chk_unused(...)

  x <- subset(x, delta = delta)
  hc <- .ssd_hc_fitdists(
    x,
    percent,
    ci = ci,
    level = level,
    nboot = nboot,
    average = average,
    min_pboot = min_pboot,
    parametric = parametric,
    averaging_method = averaging_method,
    keep_samples = keep_samples,
    control = control
  )
  warn_min_pboot(hc, min_pboot)
}

#' @describeIn ssd_hc Hazard Concentrations for fitburrlioz Object
#' @export
#' @examples
#' fit <- ssd_fit_burrlioz(ssddata::ccme_boron)
#' ssd_hc(fit)
#'
#' @export
ssd_hc.fitburrlioz <- function(x, percent = 5, ci = FALSE, level = 0.95, nboot = 1000,
                               min_pboot = 0.99, parametric = FALSE, ...) {
  check_dim(x, values = 1L)
  chk_named(x)
  chk_subset(names(x), c("burrIII3", "invpareto", "llogis", "lgumbel"))
  chk_vector(percent)
  chk_numeric(percent)
  chk_range(percent, c(0, 100))
  chk_flag(ci)
  chk_number(level)
  chk_range(level)
  chk_whole_number(nboot)
  chk_gt(nboot)
  chk_number(min_pboot)
  chk_range(min_pboot)
  chk_flag(parametric)
  chk_unused(...)

  if (names(x) != "burrIII3" || !ci || !length(percent)) {
    class(x) <- class(x)[-1]
    return(ssd_hc(x,
      percent = percent, ci = ci, level = level,
      nboot = nboot, min_pboot = min_pboot,
      average = FALSE, parametric = parametric
    ))
  }
  hc <- .ssd_hc_burrlioz_fitdists(x,
    percent = percent, level = level, nboot = nboot,
    min_pboot = min_pboot, parametric = parametric
  )
  warn_min_pboot(hc, min_pboot)
}
