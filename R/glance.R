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

#' @export
generics::glance

.glance_tmbfit <- function(x, nobs) {
  dist <- .dist_tmbfit(x)
  npars <- npars(x)
  log_lik <- logLik(x)
  aic <- 2 * npars - 2 * log_lik
  aicc <- aic + 2 * npars * (npars + 1) / (nobs - npars - 1)

  tibble(
    dist = dist,
    npars = npars,
    nobs = nobs,
    log_lik = log_lik,
    aic = aic,
    aicc = aicc
  )
}

#' Get a tibble summarizing each distribution
#'
#' Gets a tibble with a single row for each distribution.
#'
#' @inheritParams params
#' @return A tidy tibble of the distributions.
#' @family generics
#' @seealso [`ssd_gof()`]
#' @export
#' @examples
#' fits <- ssd_fit_dists(ssddata::ccme_boron)
#' glance(fits, wt = TRUE)
glance.fitdists <- function(x, ..., wt = FALSE) {
  chk_unused(...)
  
  if(vld_flag(wt) && !wt) {
    lifecycle::deprecate_soft("2.3.1", I("glance(wt = FALSE)"), I("glance(wt = TRUE)"),
                              "Please set the `wt` argument to `glance()` to be TRUE which will rename the 'weight' column to 'wt' and then update your downstream code accordingly.")    
  }
  chk_flag(wt)
  
  nobs <- nobs(x)
  tbl <- lapply(x, .glance_tmbfit, nobs = nobs)
  tbl <- dplyr::bind_rows(tbl)
  tbl$delta <- tbl$aicc - min(tbl$aicc)
  if (is.na(tbl$delta[1]) && all(tbl$npars == tbl$npars[1])) {
    tbl$delta <- tbl$aic - min(tbl$aic)
  }
  tbl$wt <- exp(-tbl$delta / 2) / sum(exp(-tbl$delta / 2))
  if(!wt) {
    tbl <- dplyr::rename(tbl, "weight" = "wt")
  }
  tbl
}
