# distribution-fitting Specification

## Purpose

Fit one or more parametric probability distributions to species toxicity
concentration data by maximum likelihood, returning a `fitdists` object that
downstream capabilities (goodness-of-fit, model averaging, hazard
concentrations, plotting) consume.

## Requirements

### Requirement: Multi-distribution maximum-likelihood fitting
The package SHALL provide `ssd_fit_dists()` which fits each of the requested
distributions to a data frame of concentrations by maximum likelihood and
returns a `fitdists` object containing one fit per successfully converged
distribution.

#### Scenario: Fit default distributions to uncensored data
- **WHEN** `ssd_fit_dists(data)` is called with `data` containing a `Conc`
  column and at least 6 rows and no other arguments
- **THEN** the function fits the BCANZ default distribution set returned by
  `ssd_dists_bcanz()` (gamma, lgumbel, llogis, lnorm, lnorm_lnorm, weibull),
  reads concentrations from the `Conc` column, and returns an object of class
  `fitdists`

#### Scenario: Caller selects a subset of distributions
- **WHEN** `ssd_fit_dists(data, dists = c("lnorm", "llogis"))` is called
- **THEN** only the log-normal and log-logistic distributions are fit and
  returned, in that order

#### Scenario: Requested distribution is not recognised
- **WHEN** `dists` contains a name that is not in `ssd_dists_all()`
- **THEN** the function raises an error identifying the unknown distribution
  without attempting to fit any distribution

### Requirement: Input column and data requirements
The fitting functions SHALL validate the input data frame and its column
arguments, and SHALL require enough rows to identify the requested
distributions.

#### Scenario: Default column name
- **WHEN** the caller does not supply `left`
- **THEN** the function reads concentrations from the `Conc` column, and if no
  `Conc` column exists an error is raised

#### Scenario: Custom column name
- **WHEN** the caller supplies `left = "my_conc"`
- **THEN** the function reads concentrations from the `my_conc` column

#### Scenario: Fewer rows than required
- **WHEN** `data` has fewer than `nrow` rows (default `6L`)
- **THEN** the function raises an error and does not fit any distribution

### Requirement: Handling of non-converging or boundary fits
`ssd_fit_dists()` SHALL report distributions that fail to converge or whose
maximum-likelihood estimates sit on a parameter boundary, and SHALL allow the
caller to control whether such fits are kept, via the `at_boundary_ok`,
`computable`, `all_dists` and `silent` arguments.

#### Scenario: Failed distributions are dropped by default
- **WHEN** a distribution fails to converge and `all_dists = FALSE` (default)
  and `silent = FALSE` (default)
- **THEN** the returned `fitdists` object omits that distribution and a
  warning is emitted naming the distribution and the cause

#### Scenario: All distributions failed with all_dists = TRUE
- **WHEN** every distribution fails and `all_dists = TRUE`
- **THEN** the function raises an error rather than returning an empty
  `fitdists` object

#### Scenario: Silent suppresses warnings for dropped fits
- **WHEN** a distribution fails and `silent = TRUE`
- **THEN** no warning is emitted for the failure, and the distribution is
  still absent from the returned object

#### Scenario: Standard errors required but not computable
- **WHEN** `computable = TRUE` and standard errors cannot be computed for a
  given distribution
- **THEN** that distribution is dropped from the returned object (with a
  warning unless `silent = TRUE`)

### Requirement: Optional rescaling and reweighting
`ssd_fit_dists()` SHALL support numerical rescaling of concentrations prior
to fitting via the `rescale` argument and optional reweighting via the
`reweight` argument without changing the reported parameter estimates from
the caller's perspective.

#### Scenario: Rescale disabled by default
- **WHEN** the caller does not set `rescale`
- **THEN** fitting is performed on the concentration values as supplied
  (`rescale = FALSE`)

#### Scenario: Rescale enabled
- **WHEN** `rescale = TRUE` is supplied
- **THEN** the function divides concentrations by the geometric mean of the
  minimum and maximum positive finite values before fitting, fits on the
  rescaled scale, and back-transforms parameters so the returned estimates
  describe the original scale

### Requirement: Jurisdiction-specific wrapper for BCANZ guideline fitting
The package SHALL provide `ssd_fit_bcanz()` that configures
`ssd_fit_dists()` with the settings used for the BC, Canada, Australia and
New Zealand guidelines, restricting the caller to BCANZ-approved
distributions.

#### Scenario: Default BCANZ fit
- **WHEN** `ssd_fit_bcanz(data)` is called
- **THEN** the call is equivalent to `ssd_fit_dists()` with
  `dists = ssd_dists_bcanz()`, `nrow = 6L`, `reweight = FALSE`,
  `computable = FALSE`, `at_boundary_ok = TRUE`, `all_dists = FALSE` and the
  shape-parameter ranges `(0.05, 20)`

#### Scenario: Distribution outside the BCANZ set
- **WHEN** `ssd_fit_bcanz(data, dists = "burrIII3")` is called with a
  distribution not in `ssd_dists_bcanz()`
- **THEN** the function raises an error without fitting

### Requirement: Object interface of fitted results
The `fitdists` object returned by the fitting functions SHALL behave as a
named list of fits that exposes the original data via `ssd_data()` and
supports the `[` and `[[` subset operators for selecting individual
distributions.

#### Scenario: Accessing the original data
- **WHEN** `ssd_data(fit)` is called on a `fitdists` object
- **THEN** the original data frame passed to the fitting function is returned

#### Scenario: Subsetting by distribution name
- **WHEN** `fit["lnorm"]` is called on a `fitdists` object that contains a
  log-normal fit
- **THEN** a `fitdists` object containing only the `lnorm` fit is returned
