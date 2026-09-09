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

test_that("ltriangle3", {
  # bounded support: the defaults locationlog = 0, scalelog = 3, skew = 0 give
  # exp(-3), exp(3)
  test_dist("ltriangle3", lower = exp(-3), upper = exp(3))
  expect_equal(ssd_pltriangle3(1), 0.5)
  expect_equal(ssd_qltriangle3(0.5), 1)
  expect_snapshot_value(ssd_pltriangle3(2, skew = 0.5), style = "deparse")
  expect_snapshot_value(ssd_qltriangle3(0.75, skew = 0.5), style = "deparse")
  withr::with_seed(50, {
    expect_snapshot_value(ssd_rltriangle3(2, skew = 0.5), style = "deparse")
  })
})

test_that("ltriangle3 with skew = 0 is ltriangle", {
  q <- c(0.01, 0.1, 0.5, 1, 2, 10, 100)
  p <- c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1)
  expect_equal(
    ssd_pltriangle3(q, locationlog = 1, scalelog = 2),
    ssd_pltriangle(q, locationlog = 1, scalelog = 2)
  )
  expect_equal(
    ssd_qltriangle3(p, locationlog = 1, scalelog = 2),
    ssd_qltriangle(p, locationlog = 1, scalelog = 2)
  )
})

test_that("ltriangle3 support limits follow the skew", {
  # a = locationlog - (1 - skew) scalelog, b = locationlog + (1 + skew) scalelog
  expect_identical(
    ssd_qltriangle3(c(0, 1), locationlog = 1, scalelog = 2, skew = 0.5),
    exp(c(0, 4))
  )
  expect_identical(
    ssd_qltriangle3(c(0, 1), locationlog = 1, scalelog = 2, skew = -0.5),
    exp(c(-2, 2))
  )
  # right-triangular limits: the mode sits at a support limit
  expect_identical(
    ssd_qltriangle3(c(0, 1), locationlog = 1, scalelog = 2, skew = 1),
    exp(c(1, 5))
  )
  expect_identical(
    ssd_qltriangle3(c(0, 1), locationlog = 1, scalelog = 2, skew = -1),
    exp(c(-3, 1))
  )
  expect_identical(
    ssd_pltriangle3(exp(c(0, 4)), locationlog = 1, scalelog = 2, skew = 0.5),
    c(0, 1)
  )
  # the cumulative probability at the mode is the lower limb's share of the width
  expect_equal(
    ssd_pltriangle3(exp(1), locationlog = 1, scalelog = 2, skew = 0.5),
    0.25
  )
})

test_that("ltriangle3 quantile function inverts the distribution function", {
  for (skew in c(-1, -0.5, 0, 0.5, 1)) {
    p <- c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)
    q <- ssd_qltriangle3(p, locationlog = 1, scalelog = 2, skew = skew)
    expect_equal(
      ssd_pltriangle3(q, locationlog = 1, scalelog = 2, skew = skew),
      p
    )
  }
})

test_that("ltriangle3 is NaN outside the parameter space", {
  expect_identical(ssd_pltriangle3(1, skew = 1.5), NaN)
  expect_identical(ssd_pltriangle3(1, skew = -1.5), NaN)
  expect_identical(ssd_qltriangle3(0.5, skew = 1.5), NaN)
  expect_identical(ssd_pltriangle3(1, scalelog = 0), NaN)
  expect_identical(ssd_rltriangle3(2, skew = 1.5), c(NaN, NaN))
})

test_that("ltriangle3 random values respect the support and the skew", {
  withr::with_seed(42, {
    x <- ssd_rltriangle3(10000, locationlog = 1, scalelog = 2, skew = 0.5)
  })
  expect_true(all(x >= exp(0)))
  expect_true(all(x <= exp(4)))
  expect_equal(mean(log(x) <= 1), 0.25, tolerance = 0.05)
})

test_that("sltriangle3 returns finite starting values with no spread on log scale", {
  data <- data.frame(left = rep(5, 8), right = rep(5, 8), weight = rep(1, 8))
  start <- sltriangle3(data)
  expect_identical(names(start), c("locationlog", "log_scalelog", "skew"))
  expect_true(is.finite(start$locationlog))
  expect_true(is.finite(start$log_scalelog))
  expect_identical(start$skew, 0)
})

test_that("sltriangle3 starting support covers the data with the skew inside its bounds", {
  withr::with_seed(1, {
    x <- ssd_rltriangle3(50, locationlog = 1, scalelog = 2, skew = 0.7)
  })
  data <- data.frame(left = x, right = x, weight = 1)
  start <- sltriangle3(data)
  scale <- exp(start$log_scalelog)
  expect_lt(start$locationlog - (1 - start$skew) * scale, min(log(x)))
  expect_gt(start$locationlog + (1 + start$skew) * scale, max(log(x)))
  expect_gt(start$skew, -1)
  expect_lt(start$skew, 1)
  expect_gt(start$skew, 0)
})

test_that("ltriangle3 fits the boron data", {
  # the maximum likelihood fit is a right triangle with the mode at the largest
  # observation, so skew is at its bound
  fit <- ssd_fit_dists(
    ssddata::ccme_boron,
    dists = "ltriangle3",
    at_boundary_ok = TRUE
  )
  expect_s3_class(fit, "fitdists")

  glance <- glance(fit, wt = TRUE)
  expect_identical(glance$npars, 3L)
  expect_true(is.finite(glance$log_lik))

  est <- estimates(fit)
  expect_identical(
    names(est),
    c(
      "ltriangle3.weight",
      "ltriangle3.locationlog",
      "ltriangle3.scalelog",
      "ltriangle3.skew"
    )
  )
  expect_gte(est$ltriangle3.skew, -1)
  expect_lte(est$ltriangle3.skew, 1)
  expect_snapshot(est)

  expect_true(ssd_at_boundary(fit$ltriangle3))

  # the asymmetric distribution nests the symmetric one
  fit2 <- ssd_fit_dists(ssddata::ccme_boron, dists = "ltriangle")
  expect_gte(glance$log_lik, glance(fit2, wt = TRUE)$log_lik)

  # fitted support covers the data
  lower <- exp(
    est$ltriangle3.locationlog -
      (1 - est$ltriangle3.skew) * est$ltriangle3.scalelog
  )
  upper <- exp(
    est$ltriangle3.locationlog +
      (1 + est$ltriangle3.skew) * est$ltriangle3.scalelog
  )
  expect_lte(lower, min(ssddata::ccme_boron$Conc))
  expect_gte(upper, max(ssddata::ccme_boron$Conc))
})

test_that("ltriangle3 hazard concentrations are finite and monotonic across the full range", {
  fit <- ssd_fit_dists(
    ssddata::ccme_boron,
    dists = "ltriangle3",
    at_boundary_ok = TRUE
  )

  hc <- ssd_hc(fit, proportion = 1:99 / 100)
  expect_identical(nrow(hc), 99L)
  expect_true(all(is.finite(hc$est)))
  expect_false(is.unsorted(hc$est))

  pred <- predict(fit)
  expect_identical(nrow(pred), 99L)
  expect_true(all(is.finite(pred$est)))
})

test_that("ltriangle3 fit is invariant to scaling the concentrations", {
  data <- ssddata::ccme_boron
  fit <- ssd_fit_dists(data, dists = "ltriangle3", at_boundary_ok = TRUE)

  data_scaled <- data
  data_scaled$Conc <- data_scaled$Conc * 1000
  fit_scaled <- ssd_fit_dists(
    data_scaled,
    dists = "ltriangle3",
    at_boundary_ok = TRUE
  )

  est <- estimates(fit)
  est_scaled <- estimates(fit_scaled)

  # the log-likelihood has a kink at every observation so the optimizer stops
  # within a looser tolerance of the optimum than for the smooth distributions
  expect_equal(
    est_scaled$ltriangle3.scalelog,
    est$ltriangle3.scalelog,
    tolerance = 1e-2
  )
  expect_equal(
    est_scaled$ltriangle3.skew,
    est$ltriangle3.skew,
    tolerance = 1e-2
  )
  expect_equal(
    est_scaled$ltriangle3.locationlog,
    est$ltriangle3.locationlog + log(1000),
    tolerance = 1e-2
  )
  expect_equal(
    ssd_hc(fit_scaled, proportion = 0.05)$est,
    ssd_hc(fit, proportion = 0.05)$est * 1000,
    tolerance = 1e-2
  )
})

test_that("ltriangle3 recovers the skew of simulated data", {
  withr::with_seed(99, {
    conc <- ssd_rltriangle3(500, locationlog = 1, scalelog = 2, skew = 0.6)
  })
  data <- data.frame(Conc = conc, Species = paste0("sp", seq_along(conc)))
  fit <- ssd_fit_dists(data, dists = "ltriangle3")
  est <- estimates(fit)
  expect_equal(est$ltriangle3.skew, 0.6, tolerance = 0.1)
  expect_equal(est$ltriangle3.scalelog, 2, tolerance = 0.1)
  expect_equal(est$ltriangle3.locationlog, 1, tolerance = 0.2)
})

test_that("ltriangle3 flags a right-triangular fit as at the boundary", {
  # data from a right triangle with the mode at the upper limit
  withr::with_seed(7, {
    conc <- ssd_rltriangle3(200, locationlog = 2, scalelog = 1, skew = -1)
  })
  data <- data.frame(Conc = conc, Species = paste0("sp", seq_along(conc)))
  fit <- ssd_fit_dists(data, dists = "ltriangle3", at_boundary_ok = TRUE)
  est <- estimates(fit)
  expect_equal(est$ltriangle3.skew, -1, tolerance = 0.05)
  expect_true(ssd_at_boundary(fit$ltriangle3))
})

test_that("ltriangle3 fitted support covers an extreme low outlier", {
  withr::with_seed(42, {
    conc <- c(exp(stats::rnorm(30, log(10), 0.1)), 10 * 1e-7)
  })
  data <- data.frame(Conc = conc, Species = paste0("sp", seq_along(conc)))

  fit <- ssd_fit_dists(data, dists = "ltriangle3", at_boundary_ok = TRUE)
  est <- estimates(fit)
  lower <- exp(
    est$ltriangle3.locationlog -
      (1 - est$ltriangle3.skew) * est$ltriangle3.scalelog
  )
  upper <- exp(
    est$ltriangle3.locationlog +
      (1 + est$ltriangle3.skew) * est$ltriangle3.scalelog
  )

  expect_lt(lower, min(data$Conc))
  expect_gt(upper, max(data$Conc))
  expect_true(is.finite(ssd_hc(fit, proportion = 0.05, ci = FALSE)$est))
  expect_gt(
    ssd_hp(fit, conc = min(data$Conc), ci = FALSE, proportion = TRUE)$est,
    0
  )
})

test_that("ltriangle3 fits interval censored data lying outside the initial support", {
  withr::with_seed(42, {
    conc <- exp(stats::rnorm(20, log(10), 0.1))
  })
  data <- data.frame(
    Conc = c(conc, 1e-8),
    Right = c(conc, 1e-6),
    Species = paste0("sp", seq_len(21))
  )

  fit <- ssd_fit_dists(
    data,
    left = "Conc",
    right = "Right",
    dists = "ltriangle3",
    at_boundary_ok = TRUE
  )
  est <- estimates(fit)
  expect_true(is.finite(est$ltriangle3.locationlog))
  expect_true(is.finite(est$ltriangle3.scalelog))
  expect_true(is.finite(est$ltriangle3.skew))
  expect_true(is.finite(glance(fit, wt = TRUE)$log_lik))
})

test_that("ltriangle3 fitted support covers a left-censored non-detect below the data", {
  withr::with_seed(42, {
    conc <- exp(stats::rnorm(30, log(10), 0.1))
  })
  data <- data.frame(
    Conc = c(conc, 0),
    Right = c(conc, 1),
    Species = paste0("sp", seq_len(31))
  )

  fit <- ssd_fit_dists(
    data,
    left = "Conc",
    right = "Right",
    dists = "ltriangle3",
    at_boundary_ok = TRUE
  )
  est <- estimates(fit)
  lower <- exp(
    est$ltriangle3.locationlog -
      (1 - est$ltriangle3.skew) * est$ltriangle3.scalelog
  )

  expect_lt(lower, 1)
  expect_gt(ssd_hp(fit, conc = 1, ci = FALSE, proportion = TRUE)$est, 0)
})

test_that("a mixture weighted on ltriangle3 has its finite support limits", {
  expect_identical(
    ssd_qmulti(
      c(0, 1),
      ltriangle3.weight = 1,
      ltriangle3.locationlog = 1,
      ltriangle3.scalelog = 2,
      ltriangle3.skew = 0.5
    ),
    exp(c(0, 4))
  )
  expect_identical(
    ssd_qmulti(c(0, 1), ltriangle3.weight = 0.5, lnorm.weight = 0.5),
    c(0, Inf)
  )
})

test_that("ltriangle3 hazard concentrations at proportion 0 and 1 are the fitted support", {
  fit <- ssd_fit_dists(
    ssddata::ccme_boron,
    dists = "ltriangle3",
    at_boundary_ok = TRUE
  )
  est <- estimates(fit)
  support <- exp(
    est$ltriangle3.locationlog +
      c(-(1 - est$ltriangle3.skew), 1 + est$ltriangle3.skew) *
        est$ltriangle3.scalelog
  )

  hc <- ssd_hc(fit, proportion = c(0, 1), ci = FALSE)
  expect_equal(hc$est, support)
  expect_equal(
    ssd_hp(fit, conc = support, ci = FALSE, proportion = TRUE)$est,
    c(0, 1)
  )
})
