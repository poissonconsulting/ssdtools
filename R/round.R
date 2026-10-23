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

#' Round Hazard Concentrations
#'
#' Rounds hazard concentrations for reporting according to the rule
#' in the Australian and New Zealand Guidelines for Fresh and Marine Water
#' Quality (ANZG):
#'
#' - The significant figure to which a value is reported corresponds to the
#'   place of the first significant figure in its standard error (SE).
#' - Where the SE is greater than the value, the value is reported to one
#'   significant figure.
#' - Where the SE is missing, zero, or not finite (for example when
#'   `ssd_hc()` is called with `ci = FALSE`), the value is reported to three
#'   significant figures.
#' - A value whose dropped digits are exactly 5 (with no non-zero digits to
#'   the right) is rounded to the nearest even digit
#'   (Australian Standard AS 2706-2003).
#'   The rule is applied to the decimal representation of the value
#'   so that, for example, 0.15 rounds to 0.2 and 0.25 rounds to 0.2.
#'
#' Rounding should only be performed when reporting.
#' Calculations should be performed on the unrounded values.
#' The returned values are numeric and so do not retain trailing zeros.
#'
#' @param x A numeric vector of values to round.
#' @param se A numeric vector of the standard errors of `x`
#'   of length 1 or the same length as `x`.
#' @return A numeric vector of the rounded values.
#' @seealso [ssd_hc()]
#' @export
#' @examples
#' ssd_round(1.2649, se = 0.825)
#' ssd_round(1.2649, se = 2.1)
#' ssd_round(1.2649)
#' ssd_round(c(0.15, 0.25), se = 0.1)
ssd_round <- function(x, se = NA_real_) {
  chk_numeric(x)
  chk_numeric(se)
  if (!length(se) %in% c(1L, length(x))) {
    chk::abort_chk(
      "`se` must be length 1 or the same length as `x` (",
      length(x),
      ")."
    )
  }
  se <- rep_len(se, length(x))
  round_half_even(x, round_digits(x, se))
}

# The SE of a bootstrap SE is approximately se / sqrt(2 * (nboot - 1)).
# For hazard proportions significant figures are counted from the distance
# to the nearest bound so that 0.99987 is not rounded to 1.00.
round_hcp <- function(x, bound = Inf) {
  digits <- round_digits(x$est, x$se, bound)
  for (cl in c("lcl", "ucl")) {
    cl_digits <- pmax(digits, one_signif_digits(x[[cl]], bound))
    x[[cl]] <- round_half_even(x[[cl]], cl_digits)
  }
  se_se <- rep(NA_real_, nrow(x))
  boot <- !is.na(x$nboot) & x$nboot > 1
  se_se[boot] <- x$se[boot] / sqrt(2 * (x$nboot[boot] - 1))
  x$se <- round_half_even(x$se, round_digits(x$se, se_se))
  x$est <- round_half_even(x$est, digits)
  x
}

first_place <- function(x) {
  floor(log10(abs(x)))
}

bound_distance <- function(x, bound) {
  pmin(abs(x), bound - x)
}

one_signif_digits <- function(x, bound = Inf) {
  dist <- bound_distance(x, bound)
  digits <- rep(NA_real_, length(x))
  ok <- is.finite(dist) & dist > 0
  digits[ok] <- -first_place(dist[ok])
  digits
}

round_digits <- function(x, se, bound = Inf) {
  dist <- bound_distance(x, bound)
  ok <- is.finite(dist) & dist > 0
  ok_se <- ok & is.finite(se) & se > 0
  digits <- rep(NA_real_, length(x))
  digits[ok] <- 2 - first_place(dist[ok])
  digits[ok_se] <- -first_place(pmin(se[ok_se], dist[ok_se]))
  digits
}

round_half_even <- function(x, digits) {
  ok <- !is.na(digits)
  x[ok] <- vapply(
    which(ok),
    function(i) round_half_even_digits(x[i], digits[i]),
    numeric(1)
  )
  x
}

# signif() to 15 digits removes binary representation error so that
# the half-even rule applies to the decimal value (e.g. 0.15 to 0.2).
round_half_even_digits <- function(x, digits) {
  if (digits >= 0) {
    return(round(signif(x * 10^digits, 15)) / 10^digits)
  }
  round(signif(x / 10^-digits, 15)) * 10^-digits
}
