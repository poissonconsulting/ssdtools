# goodness-of-fit Specification

## Purpose

Summarise how well each fitted distribution describes the data and rank the
distributions by information criteria so callers can select, compare, and
weight them in model-averaged inference.

## Requirements

### Requirement: Goodness-of-fit summary table
The package SHALL provide `ssd_gof()` which returns a tibble with one row
per fitted distribution in a `fitdists` object reporting, at minimum, the
columns `dist`, `ad`, `ks`, `cvm`, `aic`, `aicc`, `bic`, `delta`, `wt`,
`nobs`, `npars`, `at_boundary` and `computable`.

#### Scenario: Default call on a fitdists object
- **WHEN** `ssd_gof(fits)` is called on a `fitdists` object returned by
  `ssd_fit_dists()`
- **THEN** the function returns a tibble with one row per fitted
  distribution and columns reporting Anderson-Darling (`ad`),
  Kolmogorov-Smirnov (`ks`) and Cramér-von Mises (`cvm`) statistics, plus
  the information criteria `aic`, `aicc`, `bic`, the delta and model weight
  columns `delta` and `wt`, and the bookkeeping columns `nobs`, `npars`,
  `at_boundary` and `computable`

### Requirement: Information criteria basis for model weights
The `delta` and `wt` columns returned by `ssd_gof()` SHALL be derived from
AICc for uncensored data and from AIC for censored data, because AICc is
not defined when the likelihood is a censored likelihood.

#### Scenario: Uncensored data
- **WHEN** the fit is on uncensored data
- **THEN** `delta` is the difference between each row's `aicc` and the
  minimum `aicc`, and `wt` is the corresponding AICc weight (summing to 1
  across distributions that are not `at_boundary` and are `computable`)

#### Scenario: Censored data
- **WHEN** the fit is on censored data
- **THEN** `aicc` is `NA` for every row, `delta` and `wt` are derived from
  `aic` instead, and the sum of `wt` across retained distributions is 1

### Requirement: Filtering of non-computable or boundary fits
`ssd_gof()` SHALL flag distributions whose parameters sit on a boundary
(`at_boundary = TRUE`) or for which standard errors could not be computed
(`computable = FALSE`), and the weights in the `wt` column SHALL be zero
for such distributions so that they do not contribute to model averaging.

#### Scenario: Boundary distribution contributes zero weight
- **WHEN** a distribution in the `fitdists` object has `at_boundary = TRUE`
- **THEN** its `wt` value in the `ssd_gof()` result is `0` and the
  remaining weights sum to `1`

### Requirement: Optional goodness-of-fit p-values
`ssd_gof()` SHALL accept a `pvalue` flag. When `pvalue = TRUE` it SHALL
return p-values for the goodness-of-fit statistics in place of the raw test
statistic values; when `pvalue = FALSE` (the default) it SHALL return the
raw statistics.

#### Scenario: pvalue = TRUE
- **WHEN** `ssd_gof(fits, pvalue = TRUE)` is called
- **THEN** the `ad`, `ks` and `cvm` columns contain p-values in `[0, 1]`
  rather than raw test statistics

### Requirement: Equal-parameter constraint for model averaging on censored data
`ssd_gof()` and the model-averaging code paths SHALL require, when the
fitted data are censored, that the retained distributions all have the
same number of parameters; otherwise the software SHALL raise an error
rather than produce a weighted average across distributions of different
parameter counts.

#### Scenario: Mixed parameter counts, censored data
- **WHEN** a `fitdists` object fitted to censored data contains distributions
  with different `npars`, and the caller invokes a function that would
  average across them
- **THEN** the function raises an error identifying the parameter-count
  mismatch

### Requirement: Weighted fits omit standard errors
`ssd_gof()` SHALL, when called on a `fitdists` object produced with
per-row weights, report only point parameter estimates and omit standard
errors and any statistics that rely on them.

#### Scenario: Weighted fit
- **WHEN** `ssd_gof()` is called on a `fitdists` object produced with a
  non-null `weight` argument
- **THEN** the returned tibble's `computable` column is `FALSE` for every
  row, and standard-error-dependent columns are `NA`
