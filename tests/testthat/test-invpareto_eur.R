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

test_that("invpareto_eur", {
  test_dist("invpareto_eur")
  expect_snapshot_value(ssd_pinvpareto_eur(0.5), style = "deparse")
  expect_snapshot_value(ssd_qinvpareto_eur(0.125), style = "deparse")
  withr::with_seed(50, {
    expect_snapshot_value(ssd_rinvpareto_eur(2), style = "deparse")
  })
})

test_that("invpareto_eur cdf is unbounded above unlike north american invpareto", {
  # the north american invpareto is bounded above at scale (its cdf reaches 1
  # at q = scale); the european form is unbounded above so its cdf stays
  # strictly below 1 for any finite q.
  expect_equal(ssd_pinvpareto(1, shape = 3, scale = 1), 1)
  expect_lt(ssd_pinvpareto_eur(1, shape = 3, scale = 1), 1)
  expect_lt(ssd_pinvpareto_eur(1e6, shape = 3, scale = 1), 1)
  expect_gt(ssd_pinvpareto_eur(1e6, shape = 3, scale = 1), 0.999)
})

test_that("invpareto_eur matches actuar parameterisation", {
  skip_if_not_installed("actuar")
  q <- c(0.1, 0.5, 1, 2, 10)
  expect_equal(
    ssd_pinvpareto_eur(q, shape = 2.5, scale = 1.3),
    actuar::pinvpareto(q, shape = 2.5, scale = 1.3)
  )
  expect_equal(
    ssd_qinvpareto_eur(c(0.1, 0.5, 0.9), shape = 2.5, scale = 1.3),
    actuar::qinvpareto(c(0.1, 0.5, 0.9), shape = 2.5, scale = 1.3)
  )
})

test_that("invpareto_eur fits with anon_a", {
  fit <- ssd_fit_dists(ssddata::anon_a, dists = "invpareto_eur")
  expect_s3_class(fit, "fitdists")
  tidy <- tidy(fit)
  expect_snapshot_data(tidy, "anon_a")
})

test_that("invpareto_eur mle matches actuar mle with ccme_boron", {
  skip_if_not_installed("actuar")
  fit <- ssd_fit_dists(ssddata::ccme_boron, dists = "invpareto_eur")
  est <- estimates(fit)

  x <- ssddata::ccme_boron$Conc
  nll <- function(pars) {
    -sum(actuar::dinvpareto(
      x,
      shape = exp(pars[1]),
      scale = exp(pars[2]),
      log = TRUE
    ))
  }
  opt <- stats::optim(c(log(2), log(5)), nll, method = "BFGS")
  expect_equal(est$invpareto_eur.shape, exp(opt$par[1]), tolerance = 1e-4)
  expect_equal(est$invpareto_eur.scale, exp(opt$par[2]), tolerance = 1e-4)
})

test_that("invpareto_eur gives cis with ccme_boron", {
  fit <- ssd_fit_dists(ssddata::ccme_boron, dists = "invpareto_eur")
  expect_s3_class(fit, "fitdists")
  withr::with_seed(50, {
    hc <- ssd_hc(
      fit,
      nboot = 100,
      ci = TRUE,
      ci_method = "multi_fixed",
      samples = TRUE
    )
  })
  expect_snapshot_data(hc, "hc_boron")
})

test_that("invpareto_eur ssd_hp gives cis with ccme_boron", {
  fit <- ssd_fit_dists(ssddata::ccme_boron, dists = "invpareto_eur")
  expect_s3_class(fit, "fitdists")
  withr::with_seed(50, {
    hp <- ssd_hp(
      fit,
      nboot = 100,
      ci = TRUE,
      ci_method = "multi_fixed",
      samples = TRUE,
      proportion = FALSE
    )
  })
  expect_snapshot_data(hp, "hp_boron")
})

test_that("invpareto_eur can be fitted alongside other distributions", {
  fit <- ssd_fit_dists(
    ssddata::ccme_boron,
    dists = c("lnorm", "invpareto_eur")
  )
  expect_s3_class(fit, "fitdists")
  expect_identical(names(fit), c("lnorm", "invpareto_eur"))
})
