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

test_that("ssd_plot_cdf", {
  fits <- ssd_fit_dists(ssddata::ccme_boron)

  expect_snapshot_plot(ssd_plot_cdf(fits), "fits")
  expect_snapshot_plot(ssd_plot_cdf(fits, average = TRUE), "fits_average")
  expect_snapshot_plot(ssd_plot_cdf(fits, average = NA), "fits_average_na")
})

test_that("ssd_plot_cdf deals with rescaled data", {
  fits <- ssd_fit_dists(ssddata::ccme_boron, rescale = TRUE)
  expect_snapshot_plot(ssd_plot_cdf(fits), "fits_rescale")
})

test_that("ssd_plot_cdf deals with named list", {
  expect_snapshot_plot(
    ssd_plot_cdf(list(
      llogis = c(locationlog = 2, scalelog = 1),
      lnorm = c(meanlog = 2, sdlog = 2)
    )),
    "list"
  )
})

test_that("autoplot deals with delta", {
  dists <- ssd_dists_all()
  withr::with_seed(50, {
    fits <- ssd_fit_dists(ssddata::ccme_boron, dists = dists, at_boundary_ok = TRUE, computable = FALSE)
  })
  expect_snapshot_plot(ssd_plot_cdf(fits, delta = Inf), "fits_delta")
})

test_that("ssd_plot_cdf language", {
  data <- ssddata::ccme_boron
  data$Conc <- data$Conc * 100
  fits <- ssd_fit_dists(data)

  expect_snapshot_plot(ssd_plot_cdf(fits, suffix = "="), "suffix")
  expect_snapshot_plot(ssd_plot_cdf(fits, big.mark = "_"), "fits_bigmark")
})
