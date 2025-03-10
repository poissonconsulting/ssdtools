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

test_that("ssd_plot_data ccme_boron", {
  expect_snapshot_plot(ssd_plot_data(ssddata::ccme_boron), "ccme_boron")
})

test_that("ssd_plot_data ccme_boron color", {
  expect_snapshot_plot(ssd_plot_data(ssddata::ccme_boron,
    color = "Group", label = "Species", trans = "identity",
    shift_x = 1, add_x = 10,
  ), "ccme_boron2")
})

test_that("ssd_plot_data ccme_boron language", {
  data <- ssddata::ccme_boron
  data$Conc <- data$Conc * 100
  expect_snapshot_plot(ssd_plot_data(data, suffix = " %"), "suffix")
  expect_snapshot_plot(ssd_plot_data(data), "big_mark_comma")
  expect_snapshot_plot(ssd_plot_data(data, big.mark = " "), "big_mark_space")
})
