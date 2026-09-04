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

# Define functions of the given names on the search path that error if called
# and remove them when the calling test finishes.
mask_globals <- function(names, frame = parent.frame()) {
  for (name in names) {
    local({
      nm <- name
      assign(
        nm,
        function(...) stop("search-path ", nm, " should not be called"),
        envir = .GlobalEnv
      )
    })
  }
  withr::defer(rm(list = names, envir = .GlobalEnv), envir = frame)
}

test_that("ns_fun finds exported and internal ssdtools functions", {
  expect_identical(ssdtools:::ns_fun("ssd_plnorm"), ssdtools::ssd_plnorm)
  expect_identical(ssdtools:::ns_fun("plnorm_ssd"), ssdtools:::plnorm_ssd)
  expect_identical(ssdtools:::ns_fun("slnorm"), ssdtools:::slnorm)
})

test_that("ns_fun returns NULL for functions outside the ssdtools namespace", {
  expect_null(ssdtools:::ns_fun("blnorm"))
  expect_null(ssdtools:::ns_fun("rlgumbel_ssd"))
  # imported, not defined, in ssdtools
  expect_null(ssdtools:::ns_fun("list_assign"))
  expect_null(ssdtools:::ns_fun("ks.test"))
})

test_that("ns_fun ignores functions of the same name on the search path", {
  mask_globals("plnorm_ssd")

  expect_identical(ssdtools:::ns_fun("plnorm_ssd"), ssdtools:::plnorm_ssd)
})

test_that("optional distribution hooks are not satisfied by the search path (#187)", {
  # lnorm has no b*, r*_ssd or m* hook in ssdtools; lgumbel has no r*_ssd hook.
  mask_globals(c("blnorm", "rlnorm_ssd", "rlgumbel_ssd", "mlnorm"))

  expect_false(ssdtools:::is_bounds("lnorm"))
  expect_identical(
    ssdtools:::bdist(
      "lnorm",
      data = NULL,
      min_pmix = 0,
      range_shape1 = 0,
      range_shape2 = 0
    ),
    list(lower = -Inf, upper = Inf)
  )
  expect_identical(ssdtools:::mdist("lnorm"), list())

  set.seed(10)
  r <- ssdtools::ssd_rlgumbel(3)
  set.seed(10)
  expect_identical(r, ssdtools::ssd_qlgumbel(runif(3)))
  expect_identical(ssdtools::ssd_rlnorm(0), numeric(0))
})

test_that("internal dispatch by name resolves to ssdtools with ssd_* masked (#187)", {
  mask_globals(c(
    "ssd_plnorm",
    "ssd_qlnorm",
    "ssd_rlnorm",
    "ssd_elnorm",
    "plnorm_ssd",
    "qlnorm_ssd",
    "rlnorm_ssd",
    "slnorm",
    "blnorm"
  ))

  fit <- ssd_fit_dists(ssddata::ccme_boron, dists = c("lnorm", "gamma"))
  expect_s3_class(fit, "fitdists")
  expect_identical(names(fit), c("lnorm", "gamma"))

  expect_no_error(ssd_gof(fit, wt = TRUE))
  set.seed(1)
  hc <- ssd_hc(fit, ci = TRUE, nboot = 10, average = FALSE, parametric = TRUE)
  expect_true(all(!is.na(hc$lcl)))
  set.seed(1)
  hc <- ssd_hc(fit, ci = TRUE, nboot = 10, average = TRUE, parametric = FALSE)
  expect_true(all(!is.na(hc$lcl)))
  expect_no_error(ssd_hp(fit, 1, ci = TRUE, nboot = 10, proportion = TRUE))

  expect_no_error(ssd_pmulti_fitdists(1, fit))
  expect_no_error(ssd_qmulti_fitdists(0.5, fit))
  expect_no_error(ssd_rmulti_fitdists(2, fit))
  expect_identical(ssdtools::ssd_elnorm(), list(meanlog = 0, sdlog = 1))
  expect_identical(
    ssdtools:::emulti_ssd()$lnorm,
    c(list(weight = 1 / length(ssd_dists_bcanz())), ssdtools::ssd_elnorm())
  )
})
