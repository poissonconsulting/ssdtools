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

test_that("ssd_round rounds to first significant figure of se", {
  expect_identical(ssd_round(1.2649, se = 0.825), 1.3)
  expect_identical(ssd_round(1.2649, se = 0.0825), 1.26)
  expect_identical(ssd_round(12345, se = 210), 12300)
  expect_identical(
    ssd_round(c(1.2649, 12345), se = c(0.0825, 210)),
    c(1.26, 12300)
  )
})

test_that("ssd_round uses one significant figure when se greater than x", {
  expect_identical(ssd_round(1.2649, se = 2.1), 1)
  expect_identical(ssd_round(0.0372, se = 0.05), 0.04)
  expect_identical(ssd_round(1.2649, se = 1.2649), 1)
})

test_that("ssd_round uses three significant figures when se missing", {
  expect_identical(ssd_round(1.2649), 1.26)
  expect_identical(ssd_round(1.2649, se = NA_real_), 1.26)
  expect_identical(ssd_round(1.2649, se = 0), 1.26)
  expect_identical(ssd_round(1.2649, se = Inf), 1.26)
  expect_identical(
    ssd_round(c(1.2649, 1.2649), se = c(0.825, NA)),
    c(1.3, 1.26)
  )
})

test_that("ssd_round rounds half to even on decimal value", {
  expect_identical(ssd_round(0.15, se = 0.1), 0.2)
  expect_identical(ssd_round(0.25, se = 0.1), 0.2)
  expect_identical(ssd_round(0.35, se = 0.1), 0.4)
  expect_identical(ssd_round(0.125, se = 0.01), 0.12)
  expect_identical(ssd_round(0.1251, se = 0.01), 0.13)
  expect_identical(ssd_round(25, se = 10), 20)
  expect_identical(ssd_round(35, se = 10), 40)
  expect_identical(ssd_round(2.5, se = 3), 2)
})

test_that("ssd_round handles special values", {
  expect_identical(ssd_round(numeric()), numeric())
  expect_identical(ssd_round(NA_real_, se = 1), NA_real_)
  expect_identical(ssd_round(0, se = 1), 0)
  expect_identical(ssd_round(Inf, se = 1), Inf)
  expect_identical(ssd_round(-1.2649, se = 0.0825), -1.26)
})

test_that("ssd_round errors on incompatible se", {
  expect_error(ssd_round(1:3, se = 1:2), "`se` must be length 1")
  expect_chk_error(ssd_round("1"))
  expect_chk_error(ssd_round(1, se = "1"))
})

test_that("ssd_round_hc rounds est and cls", {
  hc <- data.frame(
    dist = c("a", "b", "c"),
    est = c(1.2649, 1.2649, 1.2649),
    se = c(0.825, 2.1, NA),
    lcl = c(0.2987, 0.0372, NA),
    ucl = c(3.3812, 8.62, NA)
  )
  rounded <- ssd_round_hc(hc)
  expect_identical(rounded$dist, hc$dist)
  expect_identical(rounded$se, hc$se)
  expect_identical(rounded$est, c(1.3, 1, 1.26))
  expect_identical(rounded$lcl, c(0.3, 0.04, NA))
  expect_identical(rounded$ucl, c(3.4, 9, NA))
})

test_that("ssd_round_hc does not require cls", {
  hc <- data.frame(est = 1.2649, se = 0.825)
  expect_identical(ssd_round_hc(hc), data.frame(est = 1.3, se = 0.825))
})

test_that("ssd_round_hc errors without est and se", {
  expect_chk_error(ssd_round_hc(data.frame(est = 1)))
  expect_chk_error(ssd_round_hc(1))
})

test_that("ssd_round_hc works on ssd_hc output", {
  fit <- ssd_fit_dists(ssddata::ccme_boron, dists = "lnorm")
  hc <- ssd_hc(fit)
  rounded <- ssd_round_hc(hc)
  expect_identical(colnames(rounded), colnames(hc))
  expect_identical(rounded$est, signif(hc$est, 3))
})
