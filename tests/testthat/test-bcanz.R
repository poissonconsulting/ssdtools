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

test_that("ssd_dists_bcanz works", {
  expect_identical(
    ssd_dists_bcanz(),
    c(
      "gamma",
      "lgumbel",
      "llogis",
      "lnorm",
      "lnorm_lnorm",
      "weibull"
    )
  )
})

test_that("ssd_dists_bcanz works", {
  fit <- ssd_fit_bcanz(data = ssddata::ccme_boron)
  withr::with_seed(50, {
    hc <- ssd_hc_bcanz(fit)
    lifecycle::expect_deprecated(
      hp <- ssd_hp_bcanz(fit),
      "ssd_hp\\(proportion = FALSE\\) was deprecated"
    )
  })
  expect_snapshot_data(hc, "hc_chloride")
  expect_snapshot_data(hp, "hp_chloride")
})

test_that("ssd_dists_bcanz proportion = TRUE", {
  fit <- ssd_fit_bcanz(data = ssddata::ccme_boron)
  withr::with_seed(50, {
    hp <- ssd_hp_bcanz(fit, proportion = TRUE)
  })
  expect_snapshot_data(hp, "hp_chloride_proportion")
})

test_that("ssd_fit_bcanz drops lnorm_lnorm collapsed onto tied values", {
  # 26 of 78 values tied at the maximum; full precision required to reproduce
  data <- data.frame(Conc = c(
    882, 882, 314.152694603729, 99.934369602017597, 882, 31.410745692261401,
    217.23125610933201, 548.68470187722301, 41.240569143720698, 882,
    26.1934042762534, 882, 286.30579376186699, 882, 882, 143.44931088188901,
    94.392197608371106, 882, 56.070642028172102, 235.241033063988, 882,
    273.66899429740198, 7.5035334574408399, 67.382640559928504, 882, 882, 882,
    125.40087037586601, 882, 425.54118271294499, 262.16736467427302,
    48.664968593418102, 373.62520761574098, 200.21792509352201, 882,
    101.719906682825, 882, 87.289545820365305, 415.03393363020302, 882,
    217.90223182246601, 882, 2.1923546854067602, 882, 40.251820114145701,
    1.45440646634085, 322.115289288285, 882, 43.719013886286596, 882,
    32.775135005641097, 882, 457.80476903359499, 159.462758794525,
    45.988450304559599, 39.010971942049302, 882, 33.816123703167101,
    34.242819592093802, 3.8805474452034701, 37.929171644778599,
    146.89570409266099, 56.532739819203897, 34.8395083190237,
    261.30717414838398, 882, 235.34486470400699, 124.019428939394,
    97.889697348463102, 22.764786082144301, 259.42516653573699, 882,
    102.010527771756, 882, 1.75410407168703, 378.01409319284602, 882,
    130.50168925218799
  ))
  expect_warning(
    fit <- ssd_fit_bcanz(data),
    "Distribution 'lnorm_lnorm' failed to fit.*collapsed onto tied values"
  )
  expect_identical(names(fit), c("gamma", "lgumbel", "llogis", "lnorm", "weibull"))

  fit_no_mix <- ssd_fit_bcanz(data, dists = ssd_dists_bcanz(npars = 2))
  expect_equal(
    ssd_hc(fit, ci = FALSE)$est,
    ssd_hc(fit_no_mix, ci = FALSE)$est
  )
  expect_silent(ssd_fit_bcanz(data, silent = TRUE))
})

test_that("ssd_fit_bcanz retains lnorm_lnorm without ties", {
  fit <- ssd_fit_bcanz(ssddata::ccme_boron)
  expect_true("lnorm_lnorm" %in% names(fit))
})
