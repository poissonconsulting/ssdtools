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

test_that("weibull", {
  expect_snapshot_value(ssd_pweibull(1), style = "deparse")
  expect_snapshot_value(ssd_qweibull(0.75), style = "deparse")
  withr::with_seed(50, {
    expect_snapshot_value(ssd_rweibull(2), style = "deparse")
  })
  skip_on_cran() # ssd_fit_dists() failing on CRAN M1mac likely due to much older version of the OS and toolchain
  test_dist("weibull")
})

test_that("weibull", {
  withr::with_seed(50, {
    data <- data.frame(Conc = ssd_rweibull(1000, shape = 0.5, scale = 2))
  })
  fit <- ssd_fit_dists(data, dists = c("weibull", "lgumbel"))
  tidy <- tidy(fit)
  expect_snapshot_data(tidy, "tidy")
})

test_that("weibull works anona", {
  withr::with_seed(50, {
    fit <- ssd_fit_dists(ssddata::anon_a, dists = c("weibull", "lgumbel"))
  })
  tidy <- tidy(fit)
  expect_snapshot_data(tidy, "tidy_anona")
})

test_that("weibull works anon_e", {
  withr::with_seed(50, {
    fit <- ssd_fit_dists(ssddata::anon_e, dists = c("weibull", "lgumbel"))
  })
  tidy <- tidy(fit)
  expect_snapshot_data(tidy, "tidy_anon_e")
})

test_that("weibull bootstraps anona", {
  fit <- ssd_fit_dists(ssddata::anon_a, dists = "weibull")
  withr::with_seed(50, {
    hc <- ssd_hc(fit, nboot = 1000, ci = TRUE, ci_method = "weighted_samples", samples = TRUE)
  })
  expect_snapshot_data(hc, "hc_anona")
})
