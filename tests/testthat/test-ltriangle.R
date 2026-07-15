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

test_that("ltriangle", {
  test_dist("ltriangle")
  expect_equal(ssd_pltriangle(1), 0.5)
  expect_equal(ssd_qltriangle(0.5), 1)
  expect_snapshot_value(ssd_pltriangle(2), style = "deparse")
  expect_snapshot_value(ssd_qltriangle(0.75), style = "deparse")
  withr::with_seed(50, {
    expect_snapshot_value(ssd_rltriangle(2), style = "deparse")
  })
})

test_that("sltriangle returns finite starting values with no spread on log scale", {
  data <- data.frame(left = rep(5, 8), right = rep(5, 8), weight = rep(1, 8))
  start <- sltriangle(data)
  expect_true(is.finite(start$locationlog))
  expect_true(is.finite(start$log_scalelog))
})

test_that("ltriangle hazard concentrations are finite and monotonic across the full range", {
  fit <- ssd_fit_dists(ssddata::ccme_boron, dists = "ltriangle")

  hc <- ssd_hc(fit, proportion = 1:99 / 100)
  expect_identical(nrow(hc), 99L)
  expect_true(all(is.finite(hc$est)))
  expect_false(is.unsorted(hc$est))
})

test_that("ltriangle predict returns finite estimates across the full range", {
  fit <- ssd_fit_dists(ssddata::ccme_boron, dists = "ltriangle")

  pred <- predict(fit)
  expect_identical(nrow(pred), 99L)
  expect_true(all(is.finite(pred$est)))
})

test_that("ltriangle fit is invariant to scaling the concentrations", {
  data <- ssddata::ccme_boron
  fit <- ssd_fit_dists(data, dists = "ltriangle")

  data_scaled <- data
  data_scaled$Conc <- data_scaled$Conc * 1000
  fit_scaled <- ssd_fit_dists(data_scaled, dists = "ltriangle")

  est <- estimates(fit)
  est_scaled <- estimates(fit_scaled)

  # scaling concentrations by a constant only shifts locationlog by its log and
  # leaves scalelog unchanged, so hazard concentrations scale by the constant.
  # The scale-relative density floor keeps this invariance across magnitudes.
  expect_equal(
    est_scaled$ltriangle.scalelog,
    est$ltriangle.scalelog,
    tolerance = 1e-5
  )
  expect_equal(
    est_scaled$ltriangle.locationlog,
    est$ltriangle.locationlog + log(1000),
    tolerance = 1e-5
  )
  expect_equal(
    ssd_hc(fit_scaled, proportion = 0.05)$est,
    ssd_hc(fit, proportion = 0.05)$est * 1000,
    tolerance = 1e-4
  )
})

test_that("ltriangle fits the endosulfan species sensitivity data", {
  # endosulfan acute toxicity values (Hose and Van den Brink 2004), a public
  # species sensitivity dataset from the fitdistrplus package. The triangular
  # distribution is fitted to log-transformed toxicity data by the US EPA SSD
  # Toolbox, so this exercises ltriangle on data it is used for in practice.
  skip_if_not_installed("fitdistrplus")

  env <- new.env()
  utils::data("endosulfan", package = "fitdistrplus", envir = env)
  data <- data.frame(Conc = env$endosulfan$ATV)

  fit <- ssd_fit_dists(data, dists = "ltriangle")
  expect_s3_class(fit, "fitdists")

  est <- estimates(fit)
  expect_equal(est$ltriangle.locationlog, 2.972006, tolerance = 1e-3)
  expect_equal(est$ltriangle.scalelog, 8.020607, tolerance = 1e-3)

  glance <- glance(fit, wt = TRUE)
  expect_identical(glance$nobs, 104L)
  expect_identical(glance$npars, 2L)

  hc <- ssd_hc(fit, proportion = 0.05)
  expect_equal(hc$est, 0.08108457, tolerance = 1e-3)
})
