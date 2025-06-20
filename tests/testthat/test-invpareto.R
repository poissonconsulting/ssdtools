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

test_that("invpareto", {
  test_dist("invpareto", upadj = 1e-03)
  expect_snapshot_value(ssd_pinvpareto(0.5), style = "deparse")
  expect_snapshot_value(ssd_qinvpareto(0.125), style = "deparse")
  withr::with_seed(50, {
    expect_snapshot_value(ssd_rinvpareto(2), style = "deparse")
  })
})

test_that("invpareto fits with anon_a", {
  fit <- ssd_fit_dists(ssddata::anon_a, dists = "invpareto")
  expect_s3_class(fit, "fitdists")
  tidy <- tidy(fit)
  expect_snapshot_data(tidy, "anon_a")
})

test_that("invpareto gives cis with ccme_boron", {
  fit <- ssd_fit_dists(ssddata::ccme_boron, dists = "invpareto")
  expect_s3_class(fit, "fitdists")
  withr::with_seed(50, {
    hc <- ssd_hc(fit, nboot = 100, ci = TRUE, ci_method = "multi_fixed", samples = TRUE)
  })
  expect_snapshot_data(hc, "hc_boron")
})

test_that("invpareto ssd_hp gives cis with ccme_boron", {
  fit <- ssd_fit_dists(ssddata::ccme_boron, dists = "invpareto")
  expect_s3_class(fit, "fitdists")
  withr::with_seed(50, {
    hp <- ssd_hp(fit, nboot = 100, ci = TRUE, ci_method = "multi_fixed", samples = TRUE, proportion = FALSE)
  })
  expect_snapshot_data(hp, "hp_boron")
})

test_that("invpareto initial shape is MLEs", {
  withr::with_seed(50, {
    data <- data.frame(Conc = ssd_rinvpareto(6), weight = 1)
  })
  data$left <- data$Conc
  data$right <- data$left
  initial <- ssdtools:::sinvpareto(data)
  expect_snapshot_value(lapply(initial, exp), style = "deparse")
  fit <- ssd_fit_dists(data, dists = "invpareto")
  expect_snapshot_value(estimates(fit), style = "deparse")
})

test_that("invpareto unbiased scale estimator small n", {
  fun <- function(n) {
    exp(ssdtools:::sinvpareto(data.frame(right = ssd_rinvpareto(n)))$log_scale)
  }
  withr::with_seed(50, {
    expect_snapshot_value(mean(vapply(rep(6, 1000), fun, 1)), style = "deparse")
  })
})

test_that("invpareto biased shape estimator small n", {
  fun <- function(n) {
    exp(ssdtools:::sinvpareto(data.frame(right = ssd_rinvpareto(n)))$log_shape)
  }
  withr::with_seed(50, {
    expect_snapshot_value(mean(vapply(rep(6, 1000), fun, 1)), style = "deparse")
  })
})

test_that("invpareto unbiased scale estimator large n", {
  fun <- function(n) {
    exp(ssdtools:::sinvpareto(data.frame(right = ssd_rinvpareto(n)))$log_scale)
  }
  withr::with_seed(50, {
    expect_snapshot_value(mean(vapply(rep(1000, 1000), fun, 1)), style = "deparse")
  })
})

test_that("invpareto unbiased shape estimator large n", {
  fun <- function(n) {
    exp(ssdtools:::sinvpareto(data.frame(right = ssd_rinvpareto(n)))$log_shape)
  }
  withr::with_seed(50, {
    expect_snapshot_value(mean(vapply(rep(1000, 1000), fun, 1)), style = "deparse")
  })
})

test_that("invpareto with extreme data", {
  data <- data.frame(Conc = c(
    2.48892649039671, 2.5258371156749, 2.51281264491458,
    2.49866046657748, 2.56572740160664, 2.49440006912093, 2.4817062813665,
    2.47546618759501, 2.53571697416386, 2.50242492575677, 2.50112253589808,
    2.5287786019635, 2.57780684900776, 2.53608336578284, 2.58101156958599,
    2.47461770234486, 2.49063194551244, 2.5856619890231, 2.48695693688166,
    2.57378026021983, 2.51235308389976, 2.48522032692049, 2.49973051106759,
    2.53625648406357, 2.51192819101941, 2.48564121012588, 2.47989185141965,
    2.47104478254847, 2.53704987914894, 2.48182203478124, 2.51943279158882,
    2.47875248023764, 2.52955571948405, 2.53413505298479, 2.4857126516631,
    2.55015093854307, 2.50566701101757, 2.5134323318284, 2.49793441210188,
    2.49424215906085, 2.48960347486455, 2.55358332496617, 2.55446292958609,
    2.48210193691792, 2.46945069890001, 2.48557684661491, 2.56460608968987,
    2.53708962699444, 2.48214951933889, 2.54412439394134, 2.59518068845417,
    2.55975671870397, 2.493434223589, 2.53455956396635, 2.49737837236316,
    2.54900643026637, 2.50513718347292, 2.54882879624245, 2.51814393193009,
    2.46420777049251, 2.46410824439861, 2.52375449633473, 2.50472480352834,
    2.47468853687034, 2.49903375287477, 2.51052484516152, 2.52440831022558,
    2.48241564711347, 2.57274003332032, 2.48966764017043, 2.5690823103684,
    2.50354051434315, 2.57783696959855, 2.55278129417344, 2.49091327122561,
    2.4858726676362, 2.50704022976757, 2.60120582374815, 2.48030852436464,
    2.58234455069583, 2.54629314447072, 2.52650700793897, 2.4871602238994,
    2.50569757079671, 2.49183442063104, 2.50165889380711, 2.47934668379978,
    2.47510756679179, 2.53369127110563, 2.46868451852079, 2.61321699644183,
    2.52987952199996, 2.58987810707128, 2.46777896999791, 2.51447342615507,
    2.48618482994608, 2.51794970929166, 2.49716394702713, 2.49218587262049
  ))

  fit98 <- ssd_fit_dists(data[1:98, , drop = FALSE], dists = "invpareto")
  expect_equal(
    estimates(fit98),
    list(invpareto.weight = 1, invpareto.scale = 2.61422908501617, invpareto.shape = 26.0909009531098)
  )

  fit99r <- ssd_fit_dists(data, dists = "invpareto", rescale = TRUE)
  expect_equal(
    estimates(fit99r),
    list(invpareto.weight = 1, invpareto.scale = 1.03020756694085, invpareto.shape = 26.0278618888664)
  )
})
