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

weighted_mean <- function(x, wt, geometric) {
  if(geometric) {
    return(exp(weighted.mean(log(x), w = wt)))
  }
  weighted.mean(x, w = wt)
}

## TODO: becky to confirm she agrees with this transform
## https://stats.stackexchange.com/questions/123514/calculating-standard-error-after-a-log-transform
exp_se <- function(se, est) {
  se * est
}

log_se <- function(se, est) {
  se / est
}

ma_est <- function(est, wt, est_method) {
  if(est_method %in% c("arithmetic", "geometric")) {
    return(weighted_mean(est, wt, geometric = est_method == "geometric"))
  }
  NA_real_
}

maw1 <- function(se, est, est_ma, wt, adj) {
  sum(wt * sqrt(((se * adj)^2 + (est - est_ma)^2)))
}

maw2 <- function(se, est, est_ma, wt, adj) {
  sqrt(sum(wt * ((se * adj)^2 + (est - est_ma)^2)))
}

ma_se <- function(se, log_se, est, wt, adj, ci_method) {
  if(ci_method %in% "MACL") {
    return(weighted_mean(se, wt, geometric = FALSE))
  } 
  if(ci_method %in% "GMACL") {
    log_se <- weighted_mean(log_se, wt, geometric = FALSE)
    est_ma <- ma_est(est, wt = wt, est_method = "geometric")
    return(exp_se(se = log_se, est = est_ma))
  } ## FIXME: Turek and Fletcher 2021 have additional adjuster on se^2
  ## but not seeing in Burnham and Anderson
  if(ci_method %in% c("MAW1", "MAW2")) {
    est_ma <- ma_est(est, wt = wt, est_method = "arithmetic")
    if(ci_method == "MAW1") {
      return(maw1(se, est = est, est_ma = est_ma, wt = wt, adj = adj))
    } else {
      return(maw2(se, est = est, est_ma = est_ma, wt = wt, adj = adj))
    }
  }
  if(ci_method %in% c("GMAW1", "GMAW2")) {
    est_ma <- ma_est(est, wt = wt, est_method = "geometric")
    if(ci_method == "GMAW1") {
      maw_se <- maw1(log_se, est = log(est), est_ma = log(est_ma), wt = wt, adj = adj)
    } else {
      maw_se <- maw2(log_se, est = log(est), est_ma = log(est_ma), wt = wt, adj = adj)
    }
    return(exp_se(maw_se, est_ma))
  }
  NA_real_
}

ma_ci <- function(est, se, log_se, wt, df, level, ci_method) {
  tail <- 1-(1-level)/2 
  adj <- stats::qt(tail, df = df)/stats::qnorm(tail)
  se_adj <- ma_se(se = se, log_se = log_se, est = est, wt = wt,
                  adj = adj, ci_method = ci_method)
  quantiles <- stats::qnorm(c(1-tail, tail))
        
  if(ci_method %in% c("MAW1", "MAW2")) {
    est_ma <- ma_est(est, wt = wt, est_method = "arithmetic")
    ci <- est_ma + se_adj * quantiles
    ci[ci < 0] <- 0
    return(ci)
  }
  if(ci_method %in% c("GMAW1", "GMAW2")) {
    est_ma <- ma_est(est, wt = wt, est_method = "geometric")
    ci <- exp(log(est_ma) + log_se(se_adj, est_ma) * quantiles)
    return(ci)
  }
  NA_real_ 
}


ma_lcl <- function(est, lcl, se, log_se, wt, df, level, ci_method) {
  if(ci_method %in% c("MACL", "GMACL")) {
    return(weighted_mean(lcl, wt, geometric = ci_method == "GMACL"))
  }
  if(ci_method %in% c("MAW1", "MAW2", "GMAW1", "GMAW2")) {
    ci <- ma_ci(est = est, se = se, log_se = log_se, wt = wt, df = df, level = level, ci_method = ci_method)
    return(ci[1])
  }
  NA_real_ 
}

ma_ucl <- function(ucl, est, se, log_se, wt, df, level, ci_method) {
  if(ci_method %in% c("MACL", "GMACL")) {
    return(weighted_mean(ucl, wt, geometric = ci_method == "GMACL"))
  }
  if(ci_method %in% c("MAW1", "MAW2", "GMAW1", "GMAW2")) {
    ci <- ma_ci(est = est, se = se, log_se = log_se, wt = wt, df = df, level = level, ci_method = ci_method)
    return(ci[2])
  }
  NA_real_ 
}

hcp_ma2 <- function(hcp, weight, ndata, level, est_method, ci_method) {
  data_dists <- ssdtools::dist_data |>
    dplyr::select("dist", "npars")
  
  hcp <- hcp |>
    dplyr::bind_rows() |>
    dplyr::inner_join(data_dists, by = "dist") |>
    dplyr::mutate(df = ndata - .data$npars) |>
    dplyr::group_by(.data$value) |>
    dplyr::mutate(weight = weight) |>
    dplyr::summarise(
      est_ma = ma_est(.data$est, .data$weight, est_method = est_method),
      se_ma = ma_se(se = .data$se, log_se = .data$log_se, est = .data$est, wt = .data$weight, 
                    adj = 1, ci_method = ci_method),
      lcl_ma = ma_lcl(lcl = .data$lcl, est = .data$est, se = .data$se, log_se = .data$log_se, wt = .data$weight, df = .data$df, level = level, ci_method = ci_method),
      ucl_ma = ma_ucl(ucl = .data$ucl, est = .data$est, se = .data$se, log_se = .data$log_se, wt = .data$weight, df = .data$df, level = level, ci_method = ci_method),
      pboot = min(.data$pboot),
      samples = list(unlist(.data$samples))) |>
    dplyr::ungroup() |>
    dplyr::rename(
      est = "est_ma", lcl = "lcl_ma", ucl = "ucl_ma", se = "se_ma")
}

hcp_ma <- function(x, value, ci, level, nboot, est_method, min_pboot,
                   data, rescale, weighted, censoring, min_pmix,
                   range_shape1, range_shape2, parametric, control,
                   save_to, samples, ci_method, hc, fun) {
  
  hcp <- purrr::map(
    x, hcp_tmbfit, nboot = nboot, value = value, ci = ci, level = level,
    min_pboot = min_pboot, data = data, rescale = rescale, weighted = weighted, censoring = censoring,
    min_pmix = min_pmix, range_shape1 = range_shape1, range_shape2 = range_shape2,
    parametric = parametric, est_method = est_method, ci_method = ci_method, average = TRUE, control = control,
    hc = hc, save_to = save_to, samples = samples, fun = fun
  )
  weight <- glance(x, wt = TRUE)$wt
  ndata <- ndata(data)
  
  hcp_ma2(hcp, weight, ndata, level = level, est_method = est_method, ci_method = ci_method)
}
