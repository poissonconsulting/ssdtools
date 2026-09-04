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
  # bounded support: the defaults locationlog = 0, scalelog = 3 give exp(-3), exp(3)
  test_dist("ltriangle", lower = exp(-3), upper = exp(3))
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

test_that("ltriangle reproduces the US EPA SSD Toolbox triangular fit to the permethrin acute data", {
  # Permethrin acute LC50 data (PermethrinAcuteData.xlsx), reproduced from the US
  # EPA SSD Toolbox v1.0 User's Manual (EPA/600/R-18/116, Table 2), originally
  # from Fojut et al. (2012, Table 10). The SSD Toolbox is a work of the US
  # federal government and is in the public domain (17 U.S.C. sec. 105);
  data <- tibble::tribble(
    ~Genus         , ~species          , ~LC50  ,
    "Ceriodaphnia" , "dubia"           , 0.25   ,
    "Ceriodaphnia" , "dubia"           , 0.652  ,
    "Ceriodaphnia" , "dubia"           , 0.788  ,
    "Ceriodaphnia" , "dubia"           , 0.622  ,
    "Ceriodaphnia" , "dubia"           , 0.772  ,
    "Ceriodaphnia" , "dubia"           , 0.745  ,
    "Ceriodaphnia" , "dubia"           , 0.858  ,
    "Ceriodaphnia" , "dubia"           , 0.571  ,
    "Ceriodaphnia" , "dubia"           , 0.58   ,
    "Ceriodaphnia" , "dubia"           , 0.609  ,
    "Ceriodaphnia" , "dubia"           , 0.57   ,
    "Ceriodaphnia" , "dubia"           , 0.827  ,
    "Ceriodaphnia" , "dubia"           , 0.585  ,
    "Ceriodaphnia" , "dubia"           , 0.849  ,
    "Ceriodaphnia" , "dubia"           , 0.889  ,
    "Ceriodaphnia" , "dubia"           , 0.865  ,
    "Chironomus"   , "dilutus"         , 0.189  ,
    "Danio"        , "rerio"           , 2.5    ,
    "Daphnia"      , "magna"           , 0.32   ,
    "Erimonax"     , "monachus"        , 1.7    ,
    "Etheostoma"   , "fonticola"       , 3.34   ,
    "Etheostoma"   , "lepidum"         , 2.71   ,
    "Hyalella"     , "azteca"          , 0.0211 ,
    "Ictalurus"    , "punctatus"       , 5.4    ,
    "Notropis"     , "mekistocholas"   , 4.16   ,
    "Oncorhynchus" , "apache"          , 1.71   ,
    "Oncorhynchus" , "clarki henshawi" , 1.58   ,
    "Oncorhynchus" , "mykiss"          , 7      ,
    "Oreonectes"   , "immunis"         , 0.21   ,
    "Pimephales"   , "promelas"        , 9.38   ,
    "Procambarus"  , "blandingi"       , 0.21   ,
    "Procloeon"    , "Sp."             , 0.0896 ,
    "Salmo"        , "salar"           , 1.5    ,
    "Xyrauchen"    , "texanus"         , 5.95
  )

  gm <- tapply(
    data$LC50,
    paste(data$Genus, data$species),
    function(x) exp(mean(log(x)))
  )
  conc <- data.frame(Conc = as.numeric(gm))

  fit <- ssd_fit_dists(conc, dists = "ltriangle")
  expect_s3_class(fit, "fitdists")

  glance <- glance(fit, wt = TRUE)
  expect_identical(glance$nobs, 19L)
  expect_identical(glance$npars, 2L)

  est <- estimates(fit)
  expect_equal(est$ltriangle.locationlog, -0.1815457, tolerance = 1e-4)
  expect_equal(est$ltriangle.scalelog, 4.092351, tolerance = 1e-4)

  a <- (est$ltriangle.locationlog - est$ltriangle.scalelog) / log(10)
  b <- (est$ltriangle.locationlog + est$ltriangle.scalelog) / log(10)
  expect_equal(a, -1.85613, tolerance = 1e-4)
  expect_equal(b, 1.698441, tolerance = 1e-4)

  hc <- ssd_hc(fit, proportion = 0.05)
  expect_equal(hc$est, 0.05080394, tolerance = 1e-4)
})

test_that("ltriangle fitted support covers an extreme low outlier", {
  # The log-triangular distribution has bounded support, so a candidate support
  # that fails to cover an observation assigns it zero density. The likelihood
  # barrier in ll_ltriangle.hpp must make that unprofitable, otherwise the fit
  # silently discards the most sensitive species and reports it as having zero
  # probability of occurring.
  withr::with_seed(42, {
    conc <- c(exp(stats::rnorm(30, log(10), 0.1)), 10 * 1e-7)
  })
  data <- data.frame(Conc = conc, Species = paste0("sp", seq_along(conc)))

  fit <- ssd_fit_dists(data, dists = "ltriangle")
  est <- estimates(fit)
  lower <- exp(est$ltriangle.locationlog - est$ltriangle.scalelog)
  upper <- exp(est$ltriangle.locationlog + est$ltriangle.scalelog)

  expect_lt(lower, min(data$Conc))
  expect_gt(upper, max(data$Conc))
  expect_true(is.finite(ssd_hc(fit, proportion = 0.05, ci = FALSE)$est))
  # the excluded observation previously had an estimated hazard of exactly 0
  expect_gt(
    ssd_hp(fit, conc = min(data$Conc), ci = FALSE, proportion = TRUE)$est,
    0
  )
})

test_that("ltriangle fits interval censored data lying outside the initial support", {
  # `ptri_ltriangle()` returns exactly 0 and exactly 1 outside the support, so
  # the censored interval mass reaches exactly 0 and log() would be non-finite.
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
    dists = "ltriangle"
  )
  est <- estimates(fit)
  expect_true(is.finite(est$ltriangle.locationlog))
  expect_true(is.finite(est$ltriangle.scalelog))
  expect_true(is.finite(glance(fit, wt = TRUE)$log_lik))
})

test_that("ltriangle fitted support covers a left-censored non-detect below the data", {
  # A non-detect (left = 0, right = detection limit) below the candidate
  # support has interval mass exactly 0 with zero gradient, so without a
  # distance barrier on the censored branch the fit sits on a flat plateau and
  # silently excludes the most sensitive observation.
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
    dists = "ltriangle"
  )
  est <- estimates(fit)
  lower <- exp(est$ltriangle.locationlog - est$ltriangle.scalelog)

  expect_lt(lower, 1)
  expect_gt(ssd_hp(fit, conc = 1, ci = FALSE, proportion = TRUE)$est, 0)
})

test_that("ssd_qltriangle returns the support limits at p = 0 and p = 1", {
  # Bounded on both sides, so the extreme quantiles are exp(locationlog -+
  # scalelog) rather than 0 and Inf (poissonconsulting/ssdtools#195).
  expect_identical(ssd_qltriangle(c(0, 1)), exp(c(-3, 3)))
  expect_identical(
    ssd_qltriangle(c(0, 1), locationlog = 1, scalelog = 2),
    exp(c(-1, 3))
  )
})

test_that("ssd_pltriangle is exactly 0 below and 1 above the support", {
  expect_identical(ssd_pltriangle(c(0, 0.01, 100, Inf)), c(0, 0, 1, 1))
  expect_identical(ssd_pltriangle(exp(c(-3, 3))), c(0, 1))
  expect_identical(pltriangle_ssd(c(0.01, 100), 0, 3), c(0, 1))
})

test_that("ssd_qltriangle round trips through ssd_pltriangle at the support limits", {
  limits <- exp(c(-3, 3))
  expect_equal(ssd_qltriangle(ssd_pltriangle(limits)), limits)
})

test_that("a mixture weighted on ltriangle has its finite support limits", {
  expect_identical(
    ssd_qmulti(
      c(0, 1),
      ltriangle.weight = 1,
      ltriangle.locationlog = 1,
      ltriangle.scalelog = 2
    ),
    exp(c(-1, 3))
  )
  # any unbounded component makes the mixture unbounded again
  expect_identical(
    ssd_qmulti(c(0, 1), ltriangle.weight = 0.5, lnorm.weight = 0.5),
    c(0, Inf)
  )
})

test_that("ltriangle hazard concentrations at proportion 0 and 1 are the fitted support", {
  fit <- ssd_fit_dists(ssddata::ccme_boron, dists = "ltriangle")
  est <- estimates(fit)
  support <- exp(est$ltriangle.locationlog + c(-1, 1) * est$ltriangle.scalelog)

  hc <- ssd_hc(fit, proportion = c(0, 1), ci = FALSE)
  expect_equal(hc$est, support)
  expect_identical(
    ssd_hp(fit, conc = support, ci = FALSE, proportion = TRUE)$est,
    c(0, 1)
  )
})
