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

#' BCANZ Distributions
#'
#' Gets a character vector of the names of the distributions
#' adopted by BC, Canada, Australia and New Zealand for official guidelines.
#'
#' @inheritParams params
#' @return A unique, sorted character vector of the distributions.
#' @seealso [`ssd_dists()`]
#' @family dists BCANZ
#' @export
#'
#' @examples
#' ssd_dists_bcanz()
#' ssd_dists_bcanz(npars = 2)
ssd_dists_bcanz <- function(npars = c(2L, 5L)) {
  chk_whole_numeric(npars)
  chk_not_any_na(npars)
  chk_unique(npars)
  check_dim(npars, values = 1:2)
  chk_subset(npars, c(2L, 5L))
  
  ssd_dists(bcanz = TRUE, npars = npars)
}

#' Fit BCANZ Distributions
#'
#' Fits distributions using settings adopted by
#' BC, Canada, Australia and New Zealand for official guidelines.
#'
#' @inheritParams params
#' @return An object of class fitdists.
#' @seealso [`ssd_fit_dists()`]
#' @family BCANZ
#' @export
#' @examples
#' ssd_fit_bcanz(ssddata::ccme_boron)
ssd_fit_bcanz <- function(data, left = "Conc", ..., dists = ssd_dists_bcanz()) {
  chk_data(data)
  chk_unused(...)
  chk_subset(dists, ssd_dists_bcanz())
  
  ## all arguments manually specified to ensure robust to 
  ## changes in default values in ssd_fit_dists()
  ssd_fit_dists(data,
                left = left,
                right = left, 
                weight = NULL,
                dists = dists,
                nrow = 6L,
                rescale = FALSE,
                reweight = FALSE,
                computable = FALSE,
                at_boundary_ok = TRUE,
                all_dists = FALSE,
                min_pmix = ssd_min_pmix(nrow(data)),
                range_shape1 = c(0.05, 20),
                range_shape2 = c(0.05, 20),
                control = list(),
                silent = FALSE
  )
}

#' BCANZ Hazard Concentrations
#'
#' Gets hazard concentrations with confidence intervals that protect
#' 1, 5, 10 and 20% of species using settings adopted by
#' BC, Canada, Australia and New Zealand for official guidelines.
#' This function can take several minutes to run with recommended 10,000 iterations.
#'
#' @inheritParams params
#' @return A tibble of corresponding hazard concentrations.
#' @seealso [`ssd_hc()`].
#' @family BCANZ
#' @export
#' @examples
#' fits <- ssd_fit_bcanz(ssddata::ccme_boron)
#' ssd_hc_bcanz(fits, nboot = 100)
ssd_hc_bcanz <- function(x, ..., nboot = 10000, min_pboot = 0.95) {
  chk_unused(...)
  ssd_hc(x,
         proportion = c(0.01, 0.05, 0.1, 0.2),
         average = TRUE,
         ci = TRUE,
         level = 0.95,
         nboot = nboot,
         min_pboot = min_pboot,
         est_method = "multi",
         ci_method = "weighted_samples",
         parametric = TRUE,
         delta = 9.21,
         samples = FALSE,
         save_to = NULL,
         control = NULL
  )
}

#' BCANZ Hazard Proportion
#'
#' Gets  proportion of species affected at specified concentration(s)
#' using settings adopted by BC, Canada, Australia and New Zealand for official guidelines.
#' This function can take several minutes to run with recommended 10,000 iterations.
#'
#' @inheritParams params
#' @inheritParams ssd_hp
#' @return A tibble of corresponding hazard concentrations.
#' @seealso [`ssd_hp()`].
#' @family BCANZ
#' @export
#' @examples
#' fits <- ssd_fit_bcanz(ssddata::ccme_boron)
#' ssd_hp_bcanz(fits, nboot = 100)
ssd_hp_bcanz <- function(x, conc = 1, ..., nboot = 10000, min_pboot = 0.95, proportion = FALSE) {
  chk_unused(...)
  
  if (missing(proportion) || isFALSE(proportion)) {
    lifecycle::deprecate_soft("2.3.1", I("ssd_hp(proportion = FALSE)"), I("ssd_hp(proportion = TRUE)"), 
                              "Please set the `proportion` argument to `ssd_hp_bcanz()` to be TRUE which will cause it to return hazard proportions instead of percentages then update your downstream code accordingly.")    
  }
  chk_flag(proportion)
  
  ssd_hp(x,
         conc = conc,
         average = TRUE,
         ci = TRUE,
         level = 0.95,
         nboot = nboot,
         min_pboot = min_pboot,
         est_method = "multi",
         ci_method = "weighted_samples",
         parametric = TRUE,
         delta = 9.21,
         proportion = proportion,
         samples = FALSE,
         save_to = NULL,
         control = NULL
  )
}
