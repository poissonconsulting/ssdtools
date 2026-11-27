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
#' Rounds values for reporting using the rule in the Australian and New Zealand
#' Guidelines for Fresh and Marine Water Quality (ANZG) page
#' "Toxicant default guideline values for water quality in aquatic ecosystems"
#' (section "Significant figures and rounding"):
#'
#' - A value is reported to the decimal place of the first significant figure
#'   of its standard error (SE) as reported, that is, after the SE is rounded
#'   to one significant figure.
#' - Where the reported SE is greater than the value, the value is reported to
#'   one significant figure.
#' - A value whose dropped digits are exactly 5 is rounded to the nearest even
#'   digit (Australian Standard AS 2706-2003).
#'   The rule is applied to the decimal representation of the value
#'   so that, for example, 0.15 rounds to 0.2 and 0.25 rounds to 0.2.
#'
#' Values whose SE is missing, zero, or not finite are not rounded and a
#' warning is given.
#' Calculations should be performed on the unrounded values.
#' The returned values are numeric and so do not retain trailing zeros.
#'
#' @param x A numeric vector of values to round.
#' @param se A non-negative numeric vector of the standard errors of `x`
#'   of length 1 or the same length as `x`.
#' @return A numeric vector of the rounded values.
#' @seealso [ssd_hc()]
#' @export
#' @examples
#' ssd_round(5.6738, se = 0.4536)
#' ssd_round(1713, se = 1200)
#' ssd_round(c(7.45, 7.55), se = 0.1)
ssd_round <- function(x, se) {
  chk_numeric(x)
  chk_numeric(se)
  chk_gte(se, 0)
  if (!length(se) %in% c(1L, length(x))) {
    abort_chk(
      "`se` must be length 1 or the same length as `x` (",
      length(x),
      ")."
    )
  }
  se <- rep_len(se, length(x))
  digits <- round_digits(x, se)
  warn_unrounded(x, digits)
  round_half_even(x, digits)
}

chk_round <- function(round, ci) {
  chk_flag(ci)
  chk_flag(round)
  if (round && !ci) {
    abort_chk(
      "`round` must be FALSE when `ci` is FALSE ",
      "because the rounding is based on the standard error."
    )
  }
}

round_hc <- function(x) {
  digits <- round_digits(x$est, x$se)
  warn_unrounded(x$est, digits)
  for (cl in c("lcl", "ucl")) {
    cl_digits <- pmax(digits, one_signif_digits(x[[cl]]), na.rm = TRUE)
    cl_digits[is.na(digits)] <- NA
    x[[cl]] <- round_half_even(x[[cl]], cl_digits)
  }
  x$se <- round_half_even(x$se, one_signif_digits(x$se))
  x$est <- round_half_even(x$est, digits)
  x
}

warn_unrounded <- function(x, digits) {
  n <- sum(is.finite(x) & x != 0 & is.na(digits))
  if (n) {
    wrn(
      n,
      " value",
      if (n > 1) "s were" else " was",
      " not rounded because the standard error is missing, zero, or not finite."
    )
  }
  invisible(n)
}

first_place <- function(x) {
  floor(log10(abs(x)))
}

one_signif_digits <- function(x) {
  digits <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x != 0
  digits[ok] <- -first_place(x[ok])
  digits
}

round_digits <- function(x, se) {
  ok <- is.finite(x) & x != 0 & is.finite(se) & se > 0
  se <- round_half_even(se, one_signif_digits(se))
  digits <- rep(NA_real_, length(x))
  digits[ok] <- -first_place(pmin(se[ok], abs(x[ok])))
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
