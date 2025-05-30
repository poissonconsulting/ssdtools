# Copyright 2015-2023 Province of British Columbia
# Copyright 2021 Environment and Climate Change Canada
# Copyright 2023-2025 Australian Government Department of Climate Change,
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

library(ssdtools)
library(ssddata)
library(tibble)
library(usethis)

dist_data <- tibble::tribble(
  ~dist, ~bcanz, ~tails, ~npars, ~valid, ~bound,
  "burrIII3", FALSE, TRUE, 3L, TRUE, TRUE,
  "gamma", TRUE, TRUE, 2L, TRUE, FALSE,
  "gompertz", FALSE, TRUE, 2L, TRUE, FALSE,
  "invpareto", FALSE, FALSE, 2L, FALSE, FALSE,
  "lgumbel", TRUE, TRUE, 2L, TRUE, FALSE,
  "llogis", TRUE, TRUE, 2L, TRUE, FALSE,
  "llogis_llogis", FALSE, TRUE, 5L, TRUE, TRUE,
  "lnorm", TRUE, TRUE, 2L, TRUE, FALSE,
  "lnorm_lnorm", TRUE, TRUE, 5L, TRUE, TRUE,
  "weibull", TRUE, TRUE, 2L, TRUE, FALSE
)

use_data(dist_data, overwrite = TRUE)

withr::with_seed(
  99,
  fits <- ssd_fit_dists(ssddata::ccme_boron)
)

withr::with_seed(
  99,
  boron_pred <- predict(fits, ci = TRUE)
)

use_data(boron_pred, overwrite = TRUE)

## preserves old version
fits2.3 <- ssdtools:::fits2.3

use_data(fits, fits2.3, overwrite = TRUE, internal = TRUE)
