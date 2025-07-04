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

#' Subset fitdists Object
#'
#' Select a subset of distributions from a fitdists object.
#' The Akaike Information-theoretic Criterion differences are calculated after
#' selecting the distributions named in select.
#'
#' @inheritParams params
#' @export
#' @examples
#' fits <- ssd_fit_dists(ssddata::ccme_boron)
#' subset(fits, c("gamma", "lnorm"))
subset.fitdists <- function(x, select = names(x),  ..., delta = Inf, strict = TRUE) {
  if (!length(x)) {
    return(x)
  }
  chk_unused(...)
  chk_s3_class(select, "character")
  chk_vector(select)
  chk_unique(select)
  chk_number(delta)
  chk_gte(delta)
  chk_named(x)
  chk_flag(strict)
  if(strict) {
    chk_superset(names(x), select)
  }

  attrs <- .attrs_fitdists(x)

  class <- class(x)
  x <- x[names(x) %in% select]

  class(x) <- class
  .attrs_fitdists(x) <- attrs

  if (!length(x)) {
    return(x)
  }

  d <- glance(x, wt = TRUE)$delta
  x <- x[is.na(d) | abs(d) <= delta]

  class(x) <- class
  .attrs_fitdists(x) <- attrs

  x
}
