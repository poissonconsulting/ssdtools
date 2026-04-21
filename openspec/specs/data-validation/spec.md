# data-validation Specification

## Purpose

Validate the shape and contents of data passed into the fitting and
plotting functions and expose inspector functions that let callers query
the fitted object's parameters, observation counts and diagnostics.

## Requirements

### Requirement: Data frame structural requirements
All public functions that consume a concentration data frame SHALL verify
that the input is a data frame with the named `left` column (and, if
supplied, the named `right` and `weight` columns); they SHALL raise an
error when the input is not a data frame or when a named column is
missing.

#### Scenario: Input is not a data frame
- **WHEN** a function such as `ssd_fit_dists()`, `ssd_plot()` or
  `ssd_is_censored()` is called with a non-data-frame argument
- **THEN** the function raises an error identifying the expected type

#### Scenario: Named left column missing
- **WHEN** `ssd_fit_dists(data, left = "Conc")` is called with a data frame
  lacking a `Conc` column
- **THEN** the function raises an error identifying the missing column

### Requirement: Concentration values must be finite and non-negative
Concentration columns SHALL be numeric; every `left` value SHALL be
finite and non-negative, and every `right` value SHALL be greater than or
equal to the corresponding `left` value (with `Inf` permitted to denote a
right-censored upper bound).

#### Scenario: Negative concentration
- **WHEN** the `left` column contains a negative value
- **THEN** the function raises an error

#### Scenario: right less than left
- **WHEN** the `right` column contains a value strictly less than the
  corresponding `left` value
- **THEN** the function raises an error

### Requirement: Optional weight column
When a `weight` argument names a column, that column SHALL contain
non-negative finite numeric values; rows with `weight == 0` SHALL be
treated as if they had been omitted.

#### Scenario: Negative weights rejected
- **WHEN** a `weight` column contains a negative value
- **THEN** the function raises an error

#### Scenario: Zero-weight rows
- **WHEN** the `weight` column contains zeros for some rows
- **THEN** those rows do not contribute to the likelihood or to bootstrap
  resampling

### Requirement: Fit inspection helpers
The package SHALL provide inspector functions on `fitdists` objects for
programmatic access to parameter estimates, observation counts and
parameter counts: `coef()`, `estimates()`, `tidy()`, `glance()`,
`nobs()` and `npars()`.

#### Scenario: coef returns named parameter vector
- **WHEN** `coef(fits)` is called on a `fitdists` object with a single
  fitted distribution
- **THEN** the function returns a named numeric vector of that
  distribution's maximum-likelihood parameter estimates

#### Scenario: estimates includes standard errors when available
- **WHEN** `estimates(fits)` is called on a `fitdists` object for which
  standard errors were computed
- **THEN** the returned object includes per-parameter standard errors, and
  for distributions where standard errors could not be computed the
  corresponding `se` values are `NA`

#### Scenario: tidy returns a tidy per-parameter tibble
- **WHEN** `tidy(fits)` is called
- **THEN** the function returns a tibble with one row per
  distribution-parameter combination and columns at least `dist`, `term`,
  `est`, `se`

#### Scenario: glance returns a one-row-per-distribution summary
- **WHEN** `glance(fits)` is called
- **THEN** the function returns a tibble with one row per distribution
  summarising nobs, npars and information criteria

#### Scenario: nobs and npars
- **WHEN** `nobs(fits)` and `npars(fits)` are called
- **THEN** `nobs()` returns the number of observations used to fit the
  model and `npars()` returns a named integer vector of parameter counts
  per fitted distribution

### Requirement: augment returns predictions aligned to the original data
The package SHALL provide an `augment()` method for `fitdists` objects
that returns a tibble with the original observations plus a fitted
proportion (empirical distribution value) column suitable for diagnostic
plotting.

#### Scenario: augment output shape
- **WHEN** `augment(fits)` is called
- **THEN** the returned tibble has at least one row per observation and
  contains the original `left` and `right` columns plus a fitted
  proportion column
