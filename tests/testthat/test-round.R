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

test_that("ssd_round reproduces ANZG worked examples", {
  expect_identical(
    ssd_round(
      c(5.6738, 1713, 0.335, 336.6, 938.6, 1503, 2525),
      se = c(0.4536, 1200, 3.53, 400, 500, 600, 800)
    ),
    c(5.7, 2000, 0.3, 300, 900, 1500, 2500)
  )
  expect_identical(
    ssd_round(c(7.450, 7.45, 7.550, 7.55, 7.4501, 7.5501), se = 0.1),
    c(7.4, 7.4, 7.6, 7.6, 7.5, 7.6)
  )
})

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

test_that("ssd_round warns and does not round when se missing", {
  expect_warning(
    expect_identical(ssd_round(1.2649, se = NA_real_), 1.2649),
    "^1 value was not rounded"
  )
  expect_warning(expect_identical(ssd_round(1.2649, se = 0), 1.2649))
  expect_warning(expect_identical(ssd_round(1.2649, se = Inf), 1.2649))
  expect_warning(
    expect_identical(
      ssd_round(c(1.2649, 1.2649, 1.2649), se = c(0.825, NA, NA)),
      c(1.3, 1.2649, 1.2649)
    ),
    "^2 values were not rounded"
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
  expect_identical(ssd_round(numeric(), se = 1), numeric())
  expect_identical(ssd_round(NA_real_, se = 1), NA_real_)
  expect_identical(ssd_round(0, se = 1), 0)
  expect_identical(ssd_round(0, se = NA_real_), 0)
  expect_identical(ssd_round(Inf, se = 1), Inf)
  expect_identical(ssd_round(-1.2649, se = 0.0825), -1.26)
})

test_that("ssd_round errors on incompatible se", {
  expect_error(ssd_round(1), "argument \"se\" is missing")
  expect_error(ssd_round(1:3, se = 1:2), "`se` must be length 1")
  expect_chk_error(ssd_round("1", se = 1))
  expect_chk_error(ssd_round(1, se = "1"))
})

test_that("round_hc rounds est, cls and se", {
  hc <- data.frame(
    dist = c("a", "b", "c"),
    est = c(1.2649, 1.2649, 0.0123),
    se = c(0.8254, 2.1, 0.00961),
    lcl = c(0.2987, 0.0372, 0.0000412),
    ucl = c(3.3812, 8.62, 12345.678),
    nboot = c(1000, 1000, 1000)
  )
  rounded <- round_hc(hc)
  expect_identical(rounded$dist, hc$dist)
  expect_identical(rounded$nboot, hc$nboot)
  expect_identical(rounded$est, c(1.3, 1, 0.012))
  expect_identical(rounded$lcl, c(0.3, 0.04, 0.00004))
  expect_identical(rounded$ucl, c(3.4, 9, 12345.678))
  expect_identical(rounded$se, c(0.8, 2, 0.01))
})

test_that("round_hc warns and does not round est without se", {
  hc <- data.frame(
    est = c(1.2649, 1.2649),
    se = c(0.8254, NA),
    lcl = c(0.2987, NA),
    ucl = c(3.3812, NA),
    nboot = c(1000, 1000)
  )
  expect_warning(rounded <- round_hc(hc), "^1 value was not rounded")
  expect_identical(rounded$est, c(1.3, 1.2649))
  expect_identical(rounded$se, c(0.8, NA))
  expect_identical(rounded$lcl, c(0.3, NA))
})

test_that("ssd_hc round errors when ci is FALSE", {
  fit <- ssd_fit_dists(ssddata::ccme_boron, dists = "lnorm")
  expect_chk_error(
    ssd_hc(fit, round = TRUE),
    "`round` must be FALSE when `ci` is FALSE"
  )
  expect_chk_error(ssd_hc(fit, round = NA))
  expect_identical(ssd_hc(fit, round = FALSE), ssd_hc(fit))
})

test_that("ssd_hc round with ci", {
  fit <- ssd_fit_dists(ssddata::ccme_boron, dists = "lnorm")
  withr::with_seed(10, {
    hc <- ssd_hc(fit, proportion = c(0.05, 0.2), ci = TRUE, nboot = 100)
  })
  withr::with_seed(10, {
    rounded <- ssd_hc(
      fit,
      proportion = c(0.05, 0.2),
      ci = TRUE,
      nboot = 100,
      round = TRUE
    )
  })
  expect_identical(colnames(rounded), colnames(hc))
  expect_identical(rounded$est, ssd_round(hc$est, hc$se))
  expect_identical(rounded$se, signif(hc$se, 1))
})

test_that("ssd_hc_bcanz round defaults to ci", {
  fit <- ssd_fit_dists(ssddata::ccme_boron, dists = "lnorm")
  expect_identical(
    ssd_hc_bcanz(fit),
    ssd_hc_bcanz(fit, round = FALSE)
  )
  expect_chk_error(ssd_hc_bcanz(fit, round = TRUE))
  withr::with_seed(10, {
    hc <- ssd_hc_bcanz(fit, ci = TRUE, nboot = 100, round = FALSE)
  })
  withr::with_seed(10, {
    rounded <- ssd_hc_bcanz(fit, ci = TRUE, nboot = 100)
  })
  expect_identical(rounded$est, ssd_round(hc$est, hc$se))
})
