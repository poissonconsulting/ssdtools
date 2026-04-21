# censored-data Specification

## Purpose

Represent toxicity observations as intervals — point, left-, right- or
interval-censored — and provide helpers to detect and construct censored
data sets so downstream fitting and inference functions honour the
likelihood contributions of non-detects and greater-thans.

## Requirements

### Requirement: Interval representation via left and right columns
Censored data SHALL be represented in a data frame by two numeric columns,
`left` and `right`, where an observation is a point value when
`left == right`, left-censored when `left == 0`, right-censored when
`right == Inf`, and interval-censored otherwise.

#### Scenario: Point observation
- **WHEN** a row has `left == right == 10`
- **THEN** the observation is treated as the exact value `10`

#### Scenario: Left-censored (below detection limit)
- **WHEN** a row has `left = 0` and `right = 1`
- **THEN** the observation is treated as a non-detect with true value in
  `(0, 1]`

#### Scenario: Right-censored (above reporting limit)
- **WHEN** a row has `left = 100` and `right = Inf`
- **THEN** the observation is treated as a value greater than `100`

### Requirement: Detecting censoring in a data frame
The package SHALL provide `ssd_is_censored()` to test whether a data frame
contains any censored observations, defaulting the right-hand column name to
the same column as `left` (i.e. treating data as uncensored unless a
separate `right` column is named).

#### Scenario: Data frame with single concentration column
- **WHEN** `ssd_is_censored(data.frame(Conc = c(1, 2, 3)))` is called
- **THEN** the function returns `FALSE` because `left` defaults to `"Conc"`
  and `right` defaults to `left`, so every row has `left == right`

#### Scenario: Explicit right column naming a distinct upper bound
- **WHEN** `ssd_is_censored(data.frame(Conc = 1, right = 2), right = "right")`
  is called
- **THEN** the function returns `TRUE`

#### Scenario: Empty data frame
- **WHEN** `ssd_is_censored(data.frame(Conc = numeric()))` is called
- **THEN** the function returns `NA`

### Requirement: Constructing censored data from detection bounds
The package SHALL provide `ssd_censor_data()` that takes a data frame of
point concentrations and returns a data frame in which values below the
supplied lower bound become left-censored and values above the upper bound
become right-censored.

#### Scenario: No censoring when default bounds are used
- **WHEN** `ssd_censor_data(data)` is called with `censoring = c(0, Inf)`
  (the default)
- **THEN** the returned data frame has `right` equal to `left` for every row

#### Scenario: Lower bound censoring
- **WHEN** `ssd_censor_data(data, censoring = c(5, Inf))` is called on data
  with concentrations below 5
- **THEN** rows with concentration below 5 have `left = 0` and `right = 5`
  in the returned data frame, while rows at or above 5 are unchanged

#### Scenario: Invalid censoring vector
- **WHEN** `censoring` is not a length-2 numeric vector, or contains `NA`
- **THEN** `ssd_censor_data()` raises an error

### Requirement: Constraints imposed by censored data on downstream analysis
When the fitted data are censored, the package SHALL restrict downstream
operations that are not defined for censored likelihoods: model averaging
SHALL be limited to distribution sets whose members have the same number of
parameters, and confidence intervals SHALL be available only via
non-parametric bootstrap (`parametric = FALSE`).

#### Scenario: Model averaging with mixed parameter counts on censored data
- **WHEN** the caller requests model-averaged output from `ssd_hc()` on a
  `fitdists` object fitted to censored data, and the retained distributions
  do not all have the same number of parameters
- **THEN** the function raises an error rather than producing misleading
  averaged estimates

#### Scenario: Parametric bootstrap on censored data
- **WHEN** `ssd_hc()` or `ssd_hp()` is called on censored data with
  `ci = TRUE` and `parametric = TRUE`
- **THEN** the function raises an error instructing the caller to use
  `parametric = FALSE`
