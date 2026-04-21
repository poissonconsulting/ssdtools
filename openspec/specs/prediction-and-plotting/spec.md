# prediction-and-plotting Specification

## Purpose

Generate predicted species sensitivity distribution curves from fitted
models and render them as ggplot2 objects, together with the observed
data, optional confidence-interval ribbons, and hazard-concentration
annotations.

## Requirements

### Requirement: Prediction across a proportion grid
The package SHALL provide a `predict()` method for `fitdists` objects that
returns a tibble of model-averaged concentration estimates across a grid
of proportions, with the same column structure as `ssd_hc()`.

#### Scenario: Default proportion grid
- **WHEN** `predict(fits)` is called on a `fitdists` object with default
  arguments
- **THEN** the returned tibble contains 99 rows at proportions
  `1:99 / 100`, with `average = TRUE`, `est_method = "multi"` and `ci`
  defaulting to `FALSE`

#### Scenario: Custom proportion grid
- **WHEN** `predict(fits, proportion = c(0.05, 0.5, 0.95))` is called
- **THEN** the returned tibble has exactly three rows at the requested
  proportions

#### Scenario: Predict with confidence intervals
- **WHEN** `predict(fits, ci = TRUE)` is called
- **THEN** the returned tibble contains non-`NA` `lcl` and `ucl` columns
  and inherits all bootstrap defaults from the bootstrap-confidence-
  intervals capability (`nboot = 1000`, `level = 0.95`, `parametric = TRUE`)

### Requirement: Predict honours the same deprecations as ssd_hc
`predict.fitdists()` SHALL accept the deprecated `percent` argument with
the same soft-deprecation behaviour as `ssd_hc()`: it emits a
soft-deprecation message and converts the value to `proportion`.

#### Scenario: Deprecated percent argument
- **WHEN** `predict(fits, percent = c(5, 50))` is called
- **THEN** a soft-deprecation message is emitted and the function behaves
  as if `proportion = c(0.05, 0.50)` had been supplied

### Requirement: SSD plot of data plus model-averaged curve
The package SHALL provide `ssd_plot()` which takes a data frame and a
predictions tibble (as returned by `predict()`) and returns a ggplot2
object showing the empirical distribution of the data and the
model-averaged predicted curve.

#### Scenario: Default ssd_plot call
- **WHEN** `ssd_plot(data, pred)` is called where `data` has a `Conc`
  column and `pred` is the output of `predict(fits, ci = TRUE)`
- **THEN** the returned object is a ggplot with a log10-transformed x axis
  labelled "Concentration", a y axis labelled "Species Affected" formatted
  as a percentage, a line for the predicted model, and a CI ribbon

#### Scenario: HC annotation
- **WHEN** `ssd_plot(data, pred, hc = 0.05)` is called
- **THEN** the plot includes a vertical reference line at the HC5 estimate

#### Scenario: Disabling the CI ribbon
- **WHEN** `ssd_plot(data, pred, ribbon = FALSE)` is called with a
  predictions tibble that contains `lcl`/`ucl`
- **THEN** the returned plot contains the predicted curve but no CI ribbon

#### Scenario: Disabling the HC line
- **WHEN** `ssd_plot(data, pred, hc = NULL)` is called
- **THEN** no HC reference line is drawn

### Requirement: CDF plot from a fitted object
The package SHALL provide `ssd_plot_cdf()` and an `autoplot()` method for
`fitdists` objects that render the cumulative distribution function of
each fitted distribution as a ggplot2 object.

#### Scenario: CDF plot
- **WHEN** `ssd_plot_cdf(fits)` is called
- **THEN** the returned ggplot has one line per distribution in `fits`
  plotted against concentration

#### Scenario: autoplot is a convenience wrapper
- **WHEN** `autoplot(fits)` is called
- **THEN** the call is equivalent to `ssd_plot_cdf(fits)` with default
  arguments

### Requirement: Empirical cumulative distribution helper
The package SHALL provide `ssd_ecd()` that, given a numeric vector of
concentrations, returns the empirical cumulative distribution values, and
`ssd_sort_data()` that returns the input data sorted by concentration so
these helpers support custom plots and diagnostics.

#### Scenario: ssd_ecd on a vector
- **WHEN** `ssd_ecd(c(1, 2, 3, 4))` is called
- **THEN** the function returns a numeric vector of the same length with
  values strictly increasing in `(0, 1)`

#### Scenario: ssd_sort_data orders by concentration
- **WHEN** `ssd_sort_data(data)` is called with `data` containing a `Conc`
  column
- **THEN** the returned data frame is sorted ascending by `Conc`
