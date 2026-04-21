# hazard-concentrations Specification

## Purpose

Estimate concentrations that protect a specified proportion of species
(HCx, e.g. HC5) and, inversely, the proportion of species affected at a
given concentration (HPx), from fitted species sensitivity distributions —
either from a single distribution or model-averaged across a `fitdists`
object.

## Requirements

### Requirement: Hazard concentration from a fitted model
The package SHALL provide `ssd_hc()` which, given a `fitdists` object and a
vector of proportions, returns a tibble of hazard concentration estimates
with one row per distribution-by-proportion combination (or one row per
proportion when `average = TRUE`).

#### Scenario: Default call
- **WHEN** `ssd_hc(fits)` is called on a `fitdists` object with default
  arguments
- **THEN** the function returns a single row for `proportion = 0.05` with
  `average = TRUE`, `ci = FALSE`, `est_method = "multi"`, and the returned
  tibble contains the columns `dist`, `proportion`, `est`, `se`, `lcl`,
  `ucl`, `wt`, `nboot`, `pboot`

#### Scenario: Multiple proportions
- **WHEN** `ssd_hc(fits, proportion = c(0.01, 0.05, 0.10))` is called with
  `average = TRUE`
- **THEN** the returned tibble has three rows, one per proportion, ordered
  as supplied

#### Scenario: Per-distribution output
- **WHEN** `ssd_hc(fits, average = FALSE)` is called
- **THEN** the returned tibble has one row per retained distribution per
  proportion, with `dist` identifying each distribution and `wt` giving
  that distribution's AIC weight

### Requirement: Valid proportion range
`ssd_hc()` SHALL accept only proportions in the closed interval `[0, 1]`
and SHALL raise an error for values outside that range.

#### Scenario: Proportion outside range
- **WHEN** `ssd_hc(fits, proportion = 1.5)` is called
- **THEN** the function raises an error without computing any estimates

#### Scenario: Deprecated percent argument
- **WHEN** the deprecated `percent` argument is supplied
- **THEN** the function emits a soft-deprecation message, converts `percent`
  to `proportion = percent / 100`, and otherwise behaves as if `proportion`
  had been supplied

### Requirement: Model averaging estimator options
`ssd_hc()` SHALL support at least the estimator methods named by
`ssd_est_methods()` via its `est_method` argument, defaulting to `"multi"`
so that `ssd_hc()` is the numerical inverse of `ssd_hp()`.

#### Scenario: Multi estimator is the default
- **WHEN** `ssd_hc(fits, average = TRUE)` is called with no explicit
  `est_method`
- **THEN** the estimator is `"multi"` — the distributions are treated as a
  single mixture distribution, and the returned estimate at proportion `p`
  is the value `c` such that the mixture's CDF at `c` equals `p`

#### Scenario: Arithmetic mean estimator
- **WHEN** `ssd_hc(fits, average = TRUE, est_method = "arithmetic")` is
  called
- **THEN** the returned estimate is the AIC-weighted arithmetic mean of the
  per-distribution HCx estimates

#### Scenario: Invalid est_method
- **WHEN** `est_method` is a string that is not returned by
  `ssd_est_methods()`
- **THEN** the function raises an error

### Requirement: Delta-AIC trimming prior to averaging
Before computing model-averaged hazard concentrations, `ssd_hc()` SHALL
drop distributions whose delta AIC exceeds the `delta` argument (default
`9.21`) so that distributions with negligible support do not contribute to
the result or to bootstrap sampling.

#### Scenario: Default delta filter
- **WHEN** `ssd_hc(fits)` is called with default arguments and one
  distribution has `delta > 9.21`
- **THEN** that distribution is excluded from the averaging and from any
  bootstrap resampling, and its weight is not reported in the output

#### Scenario: Delta filter disabled
- **WHEN** the caller supplies `delta = Inf`
- **THEN** all distributions (subject to boundary and computability
  filtering) participate in averaging

### Requirement: Inverse hazard proportion from a fitted model
The package SHALL provide `ssd_hp()` which, given a `fitdists` object and
a vector of concentrations, returns a tibble of proportions affected with
the same structure as `ssd_hc()`.

#### Scenario: Default call
- **WHEN** `ssd_hp(fits, conc = 1)` is called on a `fitdists` object with
  default arguments
- **THEN** the function returns a single row of the proportion affected at
  concentration `1` with `average = TRUE`, `ci = FALSE`,
  `est_method = "multi"` and a column structure matching `ssd_hc()` but
  with `conc` in place of `proportion`

#### Scenario: ssd_hc and ssd_hp are mutual inverses under est_method "multi"
- **WHEN** `ssd_hp(fits, conc = ssd_hc(fits, proportion = p)$est)` is called
  using default `est_method = "multi"`
- **THEN** the returned proportion equals `p` up to numerical tolerance

### Requirement: Output columns are self-describing
The tibbles returned by `ssd_hc()` and `ssd_hp()` SHALL carry the metadata
columns `est_method`, `ci_method`, `level`, `nboot` and `pboot` so callers
can distinguish results produced under different settings without
inspecting the call site.

#### Scenario: Metadata columns present even without CI
- **WHEN** `ssd_hc(fits)` is called with `ci = FALSE`
- **THEN** the returned tibble still includes the `est_method`, `ci_method`,
  `level` and `nboot` columns, with `se`, `lcl` and `ucl` set to `NA` and
  `nboot = 0`
