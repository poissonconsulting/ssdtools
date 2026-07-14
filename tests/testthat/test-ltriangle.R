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
