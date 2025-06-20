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

test_that("gof paper also", {
  fits <- ssd_fit_dists(ssddata::ccme_boron)
  
  gof_statistic <- ssd_gof(fits, wt = TRUE)
  expect_snapshot_data(gof_statistic, "gof_statistic")
  
  gof <- ssd_gof(fits, pvalue = TRUE, wt = TRUE)
  expect_snapshot_data(gof, "gof")
})

test_that("gof wt deprecated", {
  fits <- ssd_fit_dists(ssddata::ccme_boron)
  
  expect_deprecated(
    gof_statistic <- ssd_gof(fits)
  )
  expect_snapshot_data(gof_statistic, "gof_statisticdep")
})

test_that("gof wt = FALSE still deprecated message", {
  fits <- ssd_fit_dists(ssddata::ccme_boron)
  
  expect_deprecated( 
    gof_statistic <- ssd_gof(fits, wt = FALSE)
  )
  expect_snapshot_data(gof_statistic, "gof_statisticwtFALSE")
})

test_that("gof censored same parameters2", {
  data <- ssddata::ccme_boron
  data$right <- data$Conc
  data$Conc[c(3, 6, 8)] <- NA
  
  fits <- ssd_fit_dists(data, right = "right", dists = c("llogis", "lnorm"))
  
  gof_statistic <- ssd_gof(fits, wt = TRUE)
  expect_snapshot_data(gof_statistic, "gof_statistic2")
  
  gof <- ssd_gof(fits, pvalue = TRUE, wt = TRUE)
  expect_snapshot_data(gof, "gof2")
})

test_that("gof censored same parameters5", {
  data <- ssddata::ccme_boron
  data$right <- data$Conc
  data$Conc[c(3, 6, 8)] <- NA
  
  fits <- ssd_fit_dists(data, right = "right", dists = c("llogis_llogis", "lnorm_lnorm"))
  
  gof_statistic <- ssd_gof(fits, wt = TRUE)
  expect_snapshot_data(gof_statistic, "gof_statistic5")
  
  gof <- ssd_gof(fits, pvalue = TRUE, wt = TRUE)
  expect_snapshot_data(gof, "gof5")
})

test_that("gof censored same diff parameters", {
  data <- ssddata::ccme_boron
  data$right <- data$Conc
  data$Conc[c(3, 6, 8)] <- NA
  
  fits <- ssd_fit_dists(data, right = "right", dists = c("llogis", "lnorm_lnorm"))
  
  gof_statistic <- ssd_gof(fits, wt = TRUE)
  expect_snapshot_data(gof_statistic, "gof_statisticn")
  
  gof <- ssd_gof(fits, pvalue = TRUE, wt = TRUE)
  expect_snapshot_data(gof, "gofn")
})

test_that("gof fits2.3", {
  fits <- ssdtools:::fits2.3
  gof <- ssd_gof(fits, wt = TRUE)
  expect_snapshot_data(gof, "gof23")
})
