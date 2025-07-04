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

local_multisession <- function(.local_envir = parent.frame(), workers = 2) {
  oldDoPar <- doFuture::registerDoFuture()
  withr::defer_parent(with(oldDoPar, foreach::setDoPar(fun = fun, data = data, info = info)))
  oldPlan <- future::plan("future::multisession", workers = workers)
  withr::defer_parent(future::plan(oldPlan))
  invisible(oldDoPar)
}

save_png <- function(x, width = 400, height = 400) {
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = width, height = height)
  on.exit(grDevices::dev.off())
  print(x)
  
  path
}

save_csv <- function(x) {
  path <- tempfile(fileext = ".csv")
  readr::write_csv(x, path)
  path
}

expect_snapshot_plot <- function(x, name) {
  testthat::skip_on_os("windows")
  testthat::skip_on_os("linux")
  
  path <- save_png(x)
  testthat::expect_snapshot_file(path, paste0(name, ".png"))
}

expect_snapshot_boot_data <- function(x, name, digits = 6, min_pboot = 0.9, max_pboot = 1) {
  if (!is.na(min_pboot) && min_pboot > 0) {
    testthat::expect_true(all(x$pboot >= min_pboot))
  }
  if (!is.na(min_pboot) && max_pboot < 1) {
    testthat::expect_true(all(x$pboot <= max_pboot))
  }
  x$pboot <- NULL
  expect_snapshot_data(x, name, digits = digits)
}

expect_snapshot_data <- function(x, name, digits = 6) {
  fun <- function(x) if(is.numeric(x)) signif(x, digits = digits) else x
  lapply_fun <- function(x) I(lapply(x, fun))
  x <- dplyr::mutate(x, dplyr::across(where(is.numeric), fun))
  x <- dplyr::mutate(x, dplyr::across(where(is.list), lapply_fun))
  path <- save_csv(x)
  testthat::expect_snapshot_file(
    path,
    paste0(name, ".csv"),
    compare = testthat::compare_file_text
  )
}

ep <- function(text) {
  invisible(eval(parse(text = text)))
}

test_dist <- function(dist, qroottolerance = 1.490116e-08, upadj = 0, multi = FALSE) {
  if (!multi) {
    ep(glue::glue("expect_identical(ssd_p{dist}(numeric(0)), numeric(0))"))
    ep(glue::glue("expect_identical(ssd_p{dist}(NA), NA_real_)"))
    ep(glue::glue("expect_identical(ssd_p{dist}(NaN), NaN)"))
    ep(glue::glue("expect_identical(ssd_p{dist}(0), 0)"))
    ep(glue::glue("expect_identical(ssd_p{dist}(-Inf), 0)"))
    ep(glue::glue("expect_identical(ssd_p{dist}(Inf), 1)"))
    ep(glue::glue("expect_gt(ssd_p{dist}(1.000001), ssd_p{dist}(1))"))
    
    ep(glue::glue("expect_equal(ssd_p{dist}(1, log.p = TRUE), log(ssd_p{dist}(1)))"))
    ep(glue::glue("expect_equal(ssd_p{dist}(1, lower.tail = FALSE), 1- ssd_p{dist}(1))"))
    ep(glue::glue("expect_equal(ssd_p{dist}(1, lower.tail = FALSE, log.p = TRUE), log(1 - ssd_p{dist}(1)))"))
    
    ep(glue::glue("expect_identical(p{}(c(NA, NaN, 0, Inf, -Inf)),
                   c(NA, NaN, 0, Inf, -Inf))"))
    ep(glue::glue("expect_equal(ssd_p{dist}(1:2, 1:2, 3:4),
               c(ssd_p{dist}(1, 1, 3), ssd_p{dist}(2, 2, 4)))"))
    ep(glue::glue("expect_equal(ssd_p{dist}(1:2, c(1, NA), 3:4),
               c(ssd_p{dist}(1, 1, 3), NA_real_))"))
    ep(glue::glue("expect_gt(ssd_q{dist}(0.5000001), ssd_q{dist}(0.5))"))
    ep(glue::glue("expect_identical(ssd_q{dist}(log(0.75), log.p = TRUE), ssd_q{dist}(0.75))"))
    ep(glue::glue("expect_identical(ssd_q{dist}(0.75, lower.tail = FALSE), ssd_q{dist}(0.25))"))
    ep(glue::glue("expect_identical(ssd_q{dist}(log(0.75), lower.tail = FALSE, log.p = TRUE), ssd_q{dist}(0.25))"))
  } else {
    ep(glue::glue("expect_gt(ssd_q{dist}(0.5000001, lnorm.weight = 1), ssd_q{dist}(0.5, lnorm.weight = 1))"))
    ep(glue::glue("expect_identical(ssd_q{dist}(log(0.75), log.p = TRUE, lnorm.weight = 1), ssd_q{dist}(0.75, lnorm.weight = 1))"))
    ep(glue::glue("expect_identical(ssd_q{dist}(0.75, lower.tail = FALSE, lnorm.weight = 1), ssd_q{dist}(0.25, lnorm.weight = 1))"))
    ep(glue::glue("expect_identical(ssd_q{dist}(log(0.75), lower.tail = FALSE, log.p = TRUE, lnorm.weight = 1), ssd_q{dist}(0.25, lnorm.weight = 1))"))
  }
  
  ep(glue::glue("expect_identical(ssd_q{dist}(numeric(0)), numeric(0))"))
  ep(glue::glue("expect_identical(ssd_q{dist}(NA), NA_real_)"))
  ep(glue::glue("expect_identical(ssd_q{dist}(NaN), NaN)"))
  ep(glue::glue("expect_identical(ssd_q{dist}(0), 0)"))
  ep(glue::glue("expect_identical(ssd_q{dist}(1), Inf)"))
  ep(glue::glue("expect_identical(ssd_q{dist}(-1), NaN)"))
  ep(glue::glue("expect_identical(ssd_q{dist}(2), NaN)"))
  ep(glue::glue("expect_identical(ssd_q{dist}(-Inf), NaN)"))
  ep(glue::glue("expect_identical(ssd_q{dist}(Inf), NaN)"))
  ep(glue::glue("expect_identical(ssd_q{dist}(0.75, log.p = TRUE), NaN)"))
  ep(glue::glue("expect_identical(ssd_q{dist}(c(NA, NaN, 0, Inf, -Inf)), c(NA, NaN, 0, NaN, NaN))"))
  
  if (!multi) {
    ep(glue::glue("expect_identical(ssd_q{dist}(c(0.25, 0.75), 1:2, 3:4), c(ssd_q{dist}(0.25, 1, 3), ssd_q{dist}(0.75, 2, 4)))"))
    ep(glue::glue("expect_identical(ssd_q{dist}(c(0.25, 0.75), c(1,NA), 3:4), c(ssd_q{dist}(0.25, 1, 3), NA_real_))"))
    ep(glue::glue("expect_equal(ssd_q{dist}(ssd_p{dist}(c(0, 0.1, 0.5, 0.9, 0.99))), c(0, 0.1, 0.5, 0.9, 0.99), tolerance = {qroottolerance})"))
    ep(glue::glue("expect_identical(ssd_r{dist}(1, NA), NA_real_)"))
    ep(glue::glue("expect_identical(ssd_r{dist}(2, NA), c(NA, NA_real_))"))
    ep(glue::glue("expect_error(ssd_r{dist}(1, 1:2))"))
  }
  
  ep(glue::glue("expect_identical(ssd_r{dist}(numeric(0)), numeric(0))"))
  ep(glue::glue("expect_identical(ssd_r{dist}(0), numeric(0))"))
  ep(glue::glue("expect_error(ssd_r{dist}(NA))"))
  ep(glue::glue("expect_error(ssd_r{dist}(-1))"))
  
  if (!multi) {
    ep(glue::glue("expect_identical(length(ssd_r{dist}(1)), 1L)"))
    ep(glue::glue("expect_identical(length(ssd_r{dist}(2)), 2L)"))
    ep(glue::glue("expect_identical(length(ssd_r{dist}(3:4)), 2L)"))
    ep(glue::glue("expect_identical(length(ssd_r{dist}(c(NA, 1))), 2L)"))
    
    ests <- ep(glue::glue("ssd_e{dist}()"))
    testthat::expect_true(vld_list(ests))
    testthat::expect_true(vld_all(ests, vld_number))
    testthat::expect_true(vld_length(ests, length = 2L, upper = 5L))
    testthat::expect_true(vld_named(ests))
    
    withr::with_seed(97, {
      data <- data.frame(Conc = ep(glue::glue("ssd_r{dist}(1000)")))
    })
    fits <- ssd_fit_dists(data = data, dists = dist)
    tidy <- tidy(fits)
    testthat::expect_s3_class(tidy, "tbl_df")
    testthat::expect_identical(names(tidy), c("dist", "term", "est", "se"))
    testthat::expect_identical(tidy$dist[1], dist)
    tidy$lower <- tidy$est - tidy$se * 3
    tidy$upper <- tidy$est + tidy$se * 3
    
    default <- ep(glue::glue("formals(ssd_r{dist})"))
    default$n <- NULL
    default$chk <- NULL
    default <- data.frame(term = names(default), default = unlist(default))
    
    tidy <- merge(tidy, default, by = "term", all = "TRUE")
    testthat::expect_true(all(tidy$default > tidy$lower - upadj))
    testthat::expect_true(all(tidy$default < tidy$upper + upadj))
  } else {
    ep(glue::glue("expect_identical(length(ssd_r{dist}(1, lnorm.weight = 1)), 1L)"))
    ep(glue::glue("expect_identical(length(ssd_r{dist}(2, lnorm.weight = 1)), 2L)"))
    ep(glue::glue("expect_identical(length(ssd_r{dist}(3:4, lnorm.weight = 1)), 2L)"))
    ep(glue::glue("expect_identical(length(ssd_r{dist}(c(NA, 1), lnorm.weight = 1)), 2L)"))
  }
}
